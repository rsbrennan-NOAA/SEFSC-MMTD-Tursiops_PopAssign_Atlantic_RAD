#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_dir = null
params.outdir = null
params.pattern = '*.{fastq,fq,fastq.gz,fq.gz}'
params.sampleMaps = null  // Changed to sampleMaps for multiple map files
params.refDir = "/home/rbrennan/reference_genomes/bottlenose_dolphin"
params.refGenome = "${params.refDir}/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna.gz"

// Input validation
if (!params.input_dir) {
    error "Specify input directories with --input_dir '/path/to/dir1,/path/to/dir2'"
}
if (!params.outdir) {
    error "Specify output directory with --outdir"
}
if (!params.sampleMaps) {
    error "Sample maps must be specified with --sampleMaps '/path/to/map1,/path/to/map2'"
}

// Verify input directories and sample maps match
def inputDirs = params.input_dir.split(',').collect { it.trim() }
def sampleMapFiles = params.sampleMaps.split(',').collect { it.trim() }

if (inputDirs.size() != sampleMapFiles.size()) {
    error "Number of input directories (${inputDirs.size()}) must match number of sample maps (${sampleMapFiles.size()})"
}

// Verify input directories exist
inputDirs.each { dir ->
    if (!file(dir).exists()) {
        error "Input directory '${dir}' does not exist"
    }
}

// Verify sample maps exist
sampleMapFiles.each { mapFile ->
    if (!file(mapFile).exists()) {
        error "Sample map file '${mapFile}' does not exist"
    }
}

// Create output directory
file(params.outdir).mkdirs()

def writeToCustomLog(message) {
    def logFile = file("${params.outdir}/sample_processing.log")
    logFile.append("""${new Date().format("yyyy-MM-dd HH:mm:ss")} - ${message}\n""")
}

// Create a mapping from directory to its sample map
def createSampleMaps(dirs, mapFiles) {
    def dirToSampleMap = [:]
    
    // Loop through input directories with their index
    dirs.eachWithIndex { dir, idx ->
	def sampleMap = [:] // Creates an empty HashMap
        def mapFile = mapFiles[idx] // get the corresponding sample map file using the same index
        def dirName = file(dir).getName() // get just the directory name without the full path
        
        println "\nLoading sample map for directory ${dirName} from: ${mapFile}"
        // Read the sample map file line by line
        file(mapFile).eachLine { line ->
            def (barcode, sampleName) = line.split('\t')
            sampleMap[barcode] = sampleName
            println "Dir ${dirName}: Mapped barcode: '$barcode' to sample name: '$sampleName'"
        }
        dirToSampleMap[dirName] = sampleMap // add this directory's completed sample map to the main map
        println "Total mappings loaded for ${dirName}: ${sampleMap.size()}\n"
    }
    return dirToSampleMap
}

process bwamem2 {
    publishDir "${params.outdir}/${dir_name}", mode: 'move'
    tag "${dir_name}/${sample_name}"
    cpus 4

    input:
        tuple val(dir_name), val(sample_id), val(sample_name), path(reads)
        path(refGenome)
        path(ref_0123)
        path(ref_amb)
        path(ref_ann)
        path(ref_bwt)
        path(ref_pac)

    output:
        tuple val(dir_name), val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai")

    script:
    """
    # Load modules
    module load aligners/bwa-mem2/2.2.1
    module load bio/samtools/1.19

    # Set read group
    RG="@RG\\tID:${sample_name}-${dir_name}\\tSM:${sample_name}"

    # Align reads
    bwa-mem2 mem -t ${task.cpus} -R "\$RG" $refGenome $reads | \
    samtools view -S -h -u - | \
    samtools sort - > ${sample_name}.bam

    # Index BAM file
    samtools index ${sample_name}.bam
    """
}

workflow {
    // Load sample mappings for each directory
    dirToSampleMap = createSampleMaps(inputDirs, sampleMapFiles)

    // Create channel for input files from multiple directories
    reads_ch = Channel
        .fromPath(params.input_dir.split(',').collect { dir ->
            file(dir.trim()) + "/${params.pattern}"
        })
        .filter { !it.name.startsWith('Undetermined') }
        .map { file ->
            def dir_name = file.parent.name
            def parts = file.name.split('_')
            def barcode = "${parts[1]}"
            
            // Get the sample map for this directory
            def sampleMap = dirToSampleMap[dir_name]
            if (!sampleMap) {
                error "No sample map found for directory: ${dir_name}"
            }
            
            def sample_name = sampleMap[barcode]
            if (sample_name == null) {
                error """
                ERROR: No sample name mapping found for barcode: ${barcode}
                Directory: ${dir_name}
                File: ${file.name}
                Check the sample map file for this directory
                """
            }

            def message = """
                Processing file:
                - Directory: ${dir_name}
                - Original filename: ${file.name}
                - Extracted barcode: ${barcode}
                - Mapped sample name: ${sample_name}"""
            println message
            writeToCustomLog(message)

            tuple(dir_name, file.simpleName, sample_name, file)
        }

    // Reference genome and associated files
    refGenome = file(params.refGenome)
    ref_0123 = file("${params.refDir}/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna.gz.0123")
    ref_amb = file("${params.refDir}/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna.gz.amb")
    ref_ann = file("${params.refDir}/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna.gz.ann")
    ref_bwt = file("${params.refDir}/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna.gz.bwt.2bit.64")
    ref_pac = file("${params.refDir}/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna.gz.pac")

    // Run alignment
    bwamem2(reads_ch, refGenome, ref_0123, ref_amb, ref_ann, ref_bwt, ref_pac)
}
