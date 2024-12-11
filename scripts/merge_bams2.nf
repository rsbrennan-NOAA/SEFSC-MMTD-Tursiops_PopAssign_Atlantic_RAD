#!/usr/bin/env nextflow

// Define input parameters
params.input_dir = null  // Parent directory containing subdirectories
params.output_dir = null // where to output merged files
params.stats_file = null // summary stats file name, txt file
params.merge_log = "merge_summary.csv"

// merge, sort, and index BAM files for each sample
process MERGE_BAMS {
    publishDir "${params.output_dir}", mode: 'move'
    
    input:
    tuple val(sample_id), path(bam_files, stageAs: '*.bam')
    
    output:
    path "${sample_id}.merged.sorted.bam", emit: bam
    path "${sample_id}.merged.sorted.bam.bai", emit: bai
    stdout emit: merge_log
    
    script:
    def bam_files_str = bam_files.join(' ')
    """
    module load bio/samtools/1.19
    
    # Count input files
    FILE_COUNT=\$(echo "${bam_files_str}" | wc -w)
    
    # Get input directories then summarize with ; delimiter for summary file output
    INPUT_DIRS=\$(for file in ${bam_files_str}; do dirname \$file; done | sort -u | tr '\\n' ';' | sed 's/;\$//')
    
    if [ \$FILE_COUNT -eq 1 ]; then
        # If only one BAM file, copy and sort it
        samtools sort -o ${sample_id}.merged.sorted.bam ${bam_files_str}
    else
        # Merge BAM files and sort
        samtools merge -f - ${bam_files_str} | samtools sort -o ${sample_id}.merged.sorted.bam
    fi
    
    # Index the sorted BAM file
    samtools index ${sample_id}.merged.sorted.bam
    
    # add to summary file
    echo "${sample_id},\$FILE_COUNT,\$INPUT_DIRS"
    """
}

// calculate alignment statistics
process BAMSTATS {
    input:
    path bam
    
    output:
    stdout
    
    script:
    """
    BAM_PATH="${params.output_dir}/${bam}"
    
    samtools flagstat \$BAM_PATH > tmp.${bam.baseName}.txt
    TOTAL=\$(grep "in total" tmp.${bam.baseName}.txt | cut -f 1 -d " ")
    MAPPED=\$(grep "mapped (" tmp.${bam.baseName}.txt | cut -f 1 -d " " | head -n 1)
    PERCENT=\$(grep "mapped (" tmp.${bam.baseName}.txt | cut -f 2 -d "(" | cut -f 1 -d "%" | head -n 1)
    MAPPED_q=\$(samtools view -F 4 -q 20 \$BAM_PATH | wc -l)
    PERCENT_q=\$(echo "scale=2 ; \$MAPPED_q / \$TOTAL" | bc)
    echo "${bam.baseName},\$TOTAL,\$MAPPED,\$PERCENT,\$MAPPED_q,\$PERCENT_q"
    rm tmp.${bam.baseName}.txt
    """
}

workflow {
    // Create a channel that recursively finds all BAM files in subdirectories
    bam_files = Channel
        .fromPath("${params.input_dir}/**/*.bam")
        .map { file ->
            // get sample ID from file path/name
            def sample_id = file.name.toString().split('\\.')[0]
            return tuple(sample_id, file.toAbsolutePath())
        }
        .groupTuple(sort: true)

    // Log the samples found
    bam_files.count().view { count ->
        println "Number of unique samples found: $count"
    }

    // Merge, sort, and index BAM files for each sample
    merged_results = MERGE_BAMS(bam_files)

    // Collect merge logs
    merged_results.merge_log.collectFile(
        name: "${params.output_dir}/${params.merge_log}",
        newLine: true
    )

    // Run BAMSTATS on merged and sorted files
    alignment_stats = BAMSTATS(merged_results.bam)

    // Collect and process results
    alignment_stats.collectFile(
        name: "${params.output_dir}/${params.stats_file}",
        seed: "sample,total_reads,total_mapped,map_percent,mapped_q20,map_percent_q20\n",
        newLine: true
    )
}
