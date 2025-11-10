#! /usr/bin/env nextflow
/*
========================================================================================
<PIPELINE NAME>
========================================================================================
Description:
<DESCRIPTION>

Author:
Max Walk
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

def helpMessage() {
    log.info""" 
    Usage:
    nextflow run main.nf <params>

    --help                                Print this help message and exit.

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

log.info """\
    outdir -> ${params.outdir}
    """
    .stripIndent()

include { SAMPLESHEET_CHECK } from './modules.nf'
include { BOWTIE_INDEX } from './modules.nf'
include { FASTQC as FASTQC_RAW } from './modules.nf'
include { FASTQC as FASTQC_PROCESSED } from './modules.nf'
include { CUTADAPT } from './modules.nf'
include { BOWTIE_ALIGNMENT } from './modules.nf'
include { BAMQC } from './modules.nf'
include { MULTIQC } from './modules.nf'
include { RF_COUNT } from './modules.nf'
include { RF_COUNT_MUTATION_MAP } from './modules.nf'
include { RF_NORM } from './modules.nf'
include { DRACO } from './modules.nf'

workflow {

    raw_reads_ch        = SAMPLESHEET_CHECK(params.samplesheet)
                            .samples_csv
                            .splitCsv(header:true, sep:',')
                            .map {
                                row -> [row.sample, row.treatment, [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]]
                            }
    
    bowtie_index_ch     = BOWTIE_INDEX(params.reference_transcriptome)

    raw_fastqc_ch       = FASTQC_RAW(raw_reads_ch)

    trimmed_reads_ch    = CUTADAPT(raw_reads_ch)

    processed_fastqc_ch = FASTQC_PROCESSED(trimmed_reads_ch.reads)

    bam_ch              = BOWTIE_ALIGNMENT(bowtie_index_ch, trimmed_reads_ch.reads)

    bamqc_ch            = BAMQC(bam_ch.indexed_bam)
    
    MULTIQC(raw_fastqc_ch.zip.mix(processed_fastqc_ch.zip, bam_ch.log, bamqc_ch.report).collect())

    rf_counts_ch = RF_COUNT(bam_ch.indexed_bam, params.reference_transcriptome)

    rf_counts_mm_ch = RF_COUNT_MUTATION_MAP(bam_ch.indexed_bam, params.reference_transcriptome)

    draco_ch = DRACO(rf_counts_mm_ch.mm_file)

    sample_rf_counts_ch    = rf_counts_ch.counts_file
                    .map { sample_id, treatment, counts_file -> [sample_id, [(treatment):counts_file]] }
                    .groupTuple(by: 0)
                    .map{sample_id, labelled_count_files -> [sample_id, labelled_count_files.collectEntries()]}

    reactivity_ch = RF_NORM(sample_rf_counts_ch)
}