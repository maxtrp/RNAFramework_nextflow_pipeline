#! /usr/bin/env nextflow
/*
========================================================================================
RNAFramework_nextflow_pipeline
========================================================================================
Description:
Pipeline for processing RNA probing data (DMS/SHAPE-MaP) from raw fastq files to reactivities in .shape format.

Author:
Max Walk
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

def helpMessage() {
    log.info""" 
    Usage:
    nextflow run main.nf \
    --samplesheet </path/to/samplesheet.csv> \
    --reference_transcriptome </path/to/reference.fasta> \
    --probing_reagent <'dms' or 'shape'> \
    --sequencing_protocol <'paired_end' or 'single_end'>

    Pipeline parameters:                    (if not specified in the nextflow.config file)

      --samplesheet [str]                   Path to sample information file in csv format.
      --reference_transcriptome [str]       Path to reference transcriptome in fasta format.
      --probing_reagent [str]               Reagent used to probe RNA (either 'dms' or 'shape').
      --sequencing_protocol [str]           Either 'paired_end' or 'single_end'.

      --outdir [str]                        Path to output directory where the results will be saved (default: results/YYYYMMDD_HHMM/).
      --read1_adapter [str]                 Sequence of R1 3' adapter (default: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA").
      --read2_adapter [str]                 Sequence of R2 3' adapter (default: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT").
      --quality_cutoff [int]                Base quality threshold applied during adapter trimming (default: 17).
      --min_read_length [int]               Lower limit for read length (default: 40).
      --run_draco [bool]                    Whether to run ensemble deconvolution (default: false).
      --draco_subsampling [int]             Number of reads to sample for running draco (default: 10000).

      --help                                Print this help message and exit.

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

if (!params.probing_reagent) {
    exit 1, "--probing_reagent not specified (must be 'dms' or 'shape')"
}

log.info """\
    samplesheet   -> ${params.samplesheet}
    transcriptome -> ${params.reference_transcriptome}
    reagent       -> ${params.probing_reagent}
    protocol      -> ${params.sequencing_protocol}
    outdir        -> ${params.outdir}
    """
    .stripIndent()

include { SAMPLESHEET_CHECK } from './modules.nf'
include { BOWTIE_INDEX } from './modules.nf'
include { PEAR } from './modules.nf'
include { FASTQC as FASTQC_RAW } from './modules.nf'
include { FASTQC as FASTQC_PROCESSED } from './modules.nf'
include { FASTP_DEDUPLICATION } from './modules.nf'
include { CUTADAPT } from './modules.nf'
include { BOWTIE_ALIGNMENT } from './modules.nf'
include { BAMQC } from './modules.nf'
include { MULTIQC } from './modules.nf'
include { RF_COUNT } from './modules.nf'
include { RF_NORM } from './modules.nf'
include { RAW_COUNTS } from './modules.nf'
include { RF_COUNT_MUTATION_MAP_SUBSAMPLED } from './modules.nf'
include { DRACO } from './modules.nf'

workflow {

    if (params.sequencing_protocol == 'paired_end') {

        raw_reads_ch        = SAMPLESHEET_CHECK(params.samplesheet)
                                .samples_csv
                                .splitCsv(header:true, sep:',')
                                .map {
                                    row -> [row.sample, row.treatment, [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]]
                                }

        raw_reads_ch = PEAR(raw_reads_ch).fastq

    } else if (params.sequencing_protocol == 'single_end') {

        raw_reads_ch        = SAMPLESHEET_CHECK(params.samplesheet)
                                .samples_csv
                                .splitCsv(header:true, sep:',')
                                .map {
                                    row -> [row.sample, row.treatment, file(row.fastq_1, checkIfExists: true)]
                                }

    } else {
        exit 1, "must specify --sequencing_protocol (either 'single_end' or 'paired_end')"
    }
    
    bowtie_index_ch     = BOWTIE_INDEX(params.reference_transcriptome)

    raw_fastqc_ch       = FASTQC_RAW(raw_reads_ch)

    deduplicated_reads_ch = FASTP_DEDUPLICATION(raw_reads_ch)
    
    trimmed_reads_ch    = CUTADAPT(deduplicated_reads_ch.reads)

    processed_fastqc_ch = FASTQC_PROCESSED(trimmed_reads_ch.reads)

    bam_ch              = BOWTIE_ALIGNMENT(bowtie_index_ch, trimmed_reads_ch.reads)

    bamqc_ch            = BAMQC(bam_ch.indexed_bam)
    
    MULTIQC(raw_fastqc_ch.zip.mix(processed_fastqc_ch.zip, bam_ch.log, bamqc_ch.report).collect())

    rf_counts_ch = RF_COUNT(bam_ch.indexed_bam, params.reference_transcriptome)

    sample_rf_counts_ch    = rf_counts_ch.counts_files
                    .map { sample_id, treatment, rc_file, txt_file -> [sample_id, [(treatment):[rc_file, txt_file]]] }
                    .groupTuple(by: 0)
                    .map{sample_id, labelled_count_files -> [sample_id, labelled_count_files.collectEntries()]}

    reactivity_ch = RF_NORM(sample_rf_counts_ch)

    raw_counts_ch = RAW_COUNTS(sample_rf_counts_ch)

    if (params.run_draco) {
        rf_counts_mm_ch = RF_COUNT_MUTATION_MAP_SUBSAMPLED(bam_ch.indexed_bam, params.reference_transcriptome)
        draco_ch = DRACO(rf_counts_mm_ch.mm_file)
    }

}