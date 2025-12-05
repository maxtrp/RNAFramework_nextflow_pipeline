# RNAFramework_nextflow_pipeline

Pipeline for processing RNA probing data (DMS/SHAPE-MaP) from raw fastq files to reactivities in .shape format.

## Requirements

A Linux/MacOS machine is needed with the following installed and added to `PATH`:

- [Nextflow](https://www.nextflow.io/docs/latest/index.html)
- [Docker](https://www.docker.com/)
- python3

The pipeline runs most third-party tools (such as cutadapt, bowtie2, etc.) in docker containers meaning no need to install them manually.

However, if you want to perform ensemble deconvolution, [DRACO](https://github.com/dincarnato/draco) needs to be installed in a conda environment (see [DRACO configuration](#draco-configuration)).

## DRACO configuration

Due to the difficulty getting [DRACO](https://github.com/dincarnato/draco) up and running, it needs to be installed and compiled in its own conda environment. Use `DRACO_environment.yml` to install the dependencies and follow the instructions in the [DRACO documentation](https://draco-docs.readthedocs.io/en/latest/#installation). 

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

