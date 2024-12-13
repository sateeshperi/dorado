#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reads = "reads.fastq"
params.outdir = "results"

process GENERATE_OVERLAPS {
    container "aharish/dorado-public:0.8.3"
    
    input:
    path reads

    output:
    path "*_overlaps.paf", emit: overlaps

    script:
    def prefix = reads.baseName()
    """
    dorado correct ${reads} --to-paf > ${prefix}_overlaps.paf
    """
}

// process CORRECT_READS {
//     conda "bioconda::dorado=0.3.1"
//     publishDir "${params.outdir}", mode: 'copy'
//     input:
//     path reads
//     path overlaps
//     output:
//     path "corrected_reads.fasta"
//     script:
//     """
//     dorado correct ${reads} --from-paf ${overlaps} > corrected_reads.fasta
//     """
// }

workflow {
    reads_ch = channel.fromPath(params.reads)

    GENERATE_OVERLAPS(reads_ch)
    //CORRECT_READS(reads_ch, GENERATE_OVERLAPS.out.overlaps)
}