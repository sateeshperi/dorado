#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process BAM2FQ {
    container "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"

    publishDir "${params.outdir}/fastq/", mode: 'copy'

    input:
    path bam

    output:
    path "*.fq", emit: fastq

    script:
    def prefix = bam.baseName
    """
    samtools fastq --threads 28 ${bam} > ${prefix}.fq
    """
}

process FQIDX {
    container "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"

    publishDir "${params.outdir}/fastq/", mode: 'copy'

    input:
    path fastq

    output:
    path "*.fq.fai", emit: index

    script:
    def prefix = fastq.baseName
    """
    samtools fqidx ${fastq} > ${prefix}.fq.fai
    """
}

process GENERATE_OVERLAPS {
    container "aharish/dorado-public:latest"

    publishDir "${params.outdir}/corrected_reads/", mode: 'copy'
    
    input:
    path reads

    output:
    path "*.paf", emit: overlaps

    script:
    def prefix = reads.baseName
    """
    dorado correct ${reads} --to-paf > ${prefix}_overlaps.paf
    """
}

process CORRECT_READS {
    container "aharish/dorado-public:latest"

    publishDir "${params.outdir}/corrected_reads/", mode: 'copy'

    input:
    path reads
    path index
    path overlaps

    output:
    path "*_corrected_reads.fasta", emit: corrected_reads

    script:
    def prefix = reads.baseName
    """
    dorado correct ${reads} --from-paf ${overlaps} > ${prefix}_corrected_reads.fasta
    """
}

workflow {
    bam_ch      = Channel.fromPath(params.bam + '*.bam')
    reads_ch    = Channel.fromPath(params.reads)
    index_ch    = Channel.fromPath(params.index)
    overlaps_ch = Channel.fromPath(params.overlaps)

    //BAM2FQ(bam_ch)
    //FQIDX(BAM2FQ.out.fastq)
    GENERATE_OVERLAPS(BAM2FQ.out.fastq)
    //CORRECT_READS(reads_ch, index_ch, overlaps_ch)
}
