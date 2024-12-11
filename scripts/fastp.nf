#!/usr/bin/env nextflow

params.input_dir = null
params.outdir = null
params.pattern = '*.{fastq,fq,fastq.gz,fq.gz}'

if (!params.input_dir) {
    error "Specify input directories with --input_dir '/path/to/dir1,/path/to/dir2'"
}

if (!params.outdir) {
    error "Specify output directory with --outdir"
}

// Verify input directories exist
params.input_dir.split(',').each { dir ->
    if (!file(dir.trim()).exists()) {
        error "Input directory '${dir.trim()}' does not exist"
    }
}

process fastp {
    publishDir "${params.outdir}/${dir_name}", mode: 'move', pattern: '*.clip.fastq.gz'
    tag "${dir_name}/${sample_id}"
    cpus 4

    input:
        tuple val(dir_name), val(sample_id), path(read)

    output:
        tuple val(dir_name), val(sample_id), path("${sample_id}.clip.fastq.gz"), emit: trimmed_reads

    script:
    """
    module load bio/fastp/0.23.2

    fastp -i ${read} \
          -o ${sample_id}.clip.fastq.gz \
          -w ${task.cpus} \
          --cut_right \
          -l 35 \
          -h /dev/null \
          -j /dev/null \
	  --dont_eval_duplication
    """
}

workflow {
    Channel
        .fromPath(params.input_dir.split(',').collect { dir ->
            file(dir.trim()) + "/${params.pattern}"
        })
        .map { file -> 
            def dir_name = file.parent.name
            def sample_id = file.simpleName
            tuple(dir_name, sample_id, file)
        }
        .set { reads_ch }

    fastp(reads_ch)
}

