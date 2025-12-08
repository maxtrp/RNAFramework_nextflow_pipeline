# RNAFramework_nextflow_pipeline

Pipeline for processing RNA probing data (DMS/SHAPE-MaP) from raw fastq files to reactivities in .shape format.

## Requirements

A Linux/MacOS machine is needed with the following installed and added to `PATH`:

- [Nextflow](https://www.nextflow.io/docs/latest/index.html)
- [Docker](https://www.docker.com/)
- python3

The pipeline runs most third-party tools (such as cutadapt, bowtie2, etc.) in docker containers meaning no need to install them manually.

However, if you want to perform ensemble deconvolution, [DRACO](https://github.com/dincarnato/draco) needs to be installed in a conda environment (see [DRACO configuration](#draco-configuration)).

## Basic usage

There are four arguments that need to be passed to run the pipeline:

- Path to samplesheet containing absolute paths to input FASTQ files (see `samplesheet_example.csv`)
- Path to reference transcriptome in FASTA format (DNA sequence: T not U)
- Identity of reagent used to probe RNA (either `'dms'` or `'shape'`)
- Type of sequencing protocol performed (`'paired_end'` or `'single_end'`)

```
nextflow run main.nf \
--samplesheet </path/to/samplesheet.csv> \
--reference_transcriptome </path/to/reference.fasta> \
--probing_reagent <'dms' or 'shape'> \
--sequencing_protocol <'paired_end' or 'single_end'>
```

N.B. leave the `_fastq_2` columns of the samplesheet blank when specifying `--sequencing_protocol 'single_end'`.

Information on the optional arguments that can be passed to the pipeline can be obtained from the help statement:

```
nextflow run main.nf --help
```

## DRACO configuration

Due to the difficulty getting [DRACO](https://github.com/dincarnato/draco) up and running, it needs to be installed and manually compiled in its own conda environment. Use `DRACO_environment.yml` to install the dependencies and follow the instructions in the [DRACO documentation](https://draco-docs.readthedocs.io/en/latest/#installation). 

Once set up, add the path of the environment to the `nextflow.config` file:

```
...
if (params.run_draco) {

    conda.enabled = true
    
    process {
        withName: 'DRACO' {
            conda = '/path/to/conda/envs/DRACO' <- CHANGE THIS LINE
        }
    }
}
...
```

## Running the pipeline with DRACO

Ensemble deconvolution with DRACO can be activated by passing `--run_draco` along with the required arguments specified in [basic usage](#basic-usage).

DRACO takes a long time to run on high-depth samples, so the pipeline downsamples aligned reads (default of 10,000 reads) to speed up the analysis. This behaviour can be modified with the `--draco_subsampling` argument to specify a different number of reads to sample. Use  `--draco_subsampling -1` to use all reads available (i.e. no downsampling).
