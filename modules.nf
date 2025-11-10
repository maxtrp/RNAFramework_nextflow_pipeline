process SAMPLESHEET_CHECK {
    tag "${sample_sheet}"

    input:
    path sample_sheet

    output:
    path "*.csv", emit: samples_csv

    script:
    """
    check_samplesheet.py ${sample_sheet} valid_${sample_sheet}
    """
}

process BOWTIE_INDEX {
    tag "${fasta}"

    input:
    path fasta
    
    output:
    tuple val(ref_name), file("*_index.*")

    script:
    ref_name = fasta.baseName
    """
    bowtie2-build --threads $task.cpus $fasta ${ref_name}_index
    """
}

process FASTQC {
    tag "${sample_id}_${treatment}"

    input: 
    tuple val(sample_id), val(treatment), path(reads)

    output:
    path "*.zip", emit: zip

    script:
    """
    fastqc $reads
    """
}

process CUTADAPT {
    tag "${sample_id}_${treatment}"
    publishDir "${params.outdir}/logs/${sample_id}/", mode: "copy", pattern: "*.log"

    input:
    tuple val(sample_id), val(treatment), path(reads)

    output:
    tuple val(sample_id), val(treatment), path("*.fastq.gz"), emit: reads
    path "*.log"                                            , emit: log

    script:
    """
    cutadapt --cores $task.cpus --adapter ${params.read1_adapter} -A ${params.read2_adapter} \
    --quality-cutoff ${params.quality_cutoff} --minimum-length ${params.min_read_length} \
    --output trimmed_${sample_id}_${treatment}_R1.fastq.gz --paired-output trimmed_${sample_id}_${treatment}_R2.fastq.gz \
    ${reads[0]} ${reads[1]} \
    > ${sample_id}_${treatment}_cutadapt.log
    """
}

process BOWTIE_ALIGNMENT {
    tag "${sample_id}_${treatment} on ${ref_name}"
    publishDir "${params.outdir}/logs/${sample_id}/", mode: 'copy', pattern: '*.log'

    input:
    tuple val(ref_name), file(bowtie_index)
    tuple val(sample_id), val(treatment), path(reads)

    output:
    tuple val(sample_id), val(treatment), file("*.bam"), file("*.bai"), emit: indexed_bam
    path "*.log"                                                      , emit: log

    script:
    """
    bowtie2 --threads $task.cpus --norc -x ${ref_name}_index -1 ${reads[0]} -2 ${reads[1]} \
    2> ${sample_id}_${treatment}_bowtie2.log | samtools sort -o ${sample_id}_${treatment}.bam

    samtools index ${sample_id}_${treatment}.bam --output ${sample_id}_${treatment}.bam.bai
    """
}

process BAMQC {
    tag "${sample_id}_${treatment}"

    input:
    tuple val(sample_id), val(treatment), file(bam), file(bai)

    output:
    path "*_qualimap_report", emit: report

    script:
    """
    qualimap bamqc -bam ${bam} -outdir ${sample_id}_${treatment}_qualimap_report
    """
}

process MULTIQC {
    publishDir "${params.outdir}", mode: "copy", pattern: "*.html"

    input:
    path '*'

    output:
    path "*.html", emit: report

    script:
    """
    multiqc .
    """
}

process RF_COUNT {
    tag "${sample_id}_${treatment}"
    container 'dincarnato/rnaframework:latest' // run in docker container
    publishDir "${params.outdir}/logs/${sample_id}/", mode: 'copy', pattern: '*.log'

    input:
    tuple val(sample_id), val(treatment), file(bam), file(bai)
    path fasta

    output:
    tuple val(sample_id), val(treatment), file("rf_count/*.rc"), emit: counts_file
    path "*.log"                                               , emit: log

    script:
    """
    rf-count --working-threads $task.cpus \
    --count-mutations --fasta ${fasta} ${bam} \
    > ${sample_id}_${treatment}_rf_count.log
    """
}

process RF_COUNT_MUTATION_MAP {
    tag "${sample_id}_${treatment}"
    container 'dincarnato/rnaframework:latest' // run in docker container
    publishDir "${params.outdir}/logs/${sample_id}/", mode: 'copy', pattern: '*.log'

    input:
    tuple val(sample_id), val(treatment), file(bam), file(bai)
    path fasta

    output:
    tuple val(sample_id), val(treatment), file("rf_count/*.mm"), emit: mm_file
    path "*.log"                                               , emit: log

    script:
    """
    rf-count --working-threads $task.cpus \
    --max-coverage 1000 --count-mutations --mutation-map \
    --fasta ${fasta} ${bam} \
    > ${sample_id}_${treatment}_rf_count_downsampled.log
    """
}

process DRACO {
    tag "${sample_id}_${treatment}"
    // publishDir "${params.outdir}/logs/${sample_id}/", mode: 'copy', pattern: '*.log'
    
    input:
    tuple val(sample_id), val(treatment), file(mutation_map_file)

    output:
    // path "*.log"                         , emit: log


    script:
    """
    draco --mm ${mutation_map_file}
    """
}

process RF_NORM {
    tag "${sample_id}"
    container 'dincarnato/rnaframework:latest' // run in docker container
    publishDir "${params.outdir}/${sample_id}/"     , mode: 'copy', pattern: '*.shape'
    publishDir "${params.outdir}/logs/${sample_id}/", mode: 'copy', pattern: '*.log'
    
    input:
    tuple val(sample_id), val(counts_files)

    output:
    tuple val(sample_id), file('*.shape'), emit: shape_file
    path "*.log"                         , emit: log


    script:
    """
    rf-norm --processors $task.cpus --scoring-method 3 --raw \
    --treated  ${counts_files['treated']} --untreated ${counts_files['control']} \
    --output-dir rf_norm > ${sample_id}_rf_norm.log
    
    xml_to_shape.py rf_norm/*.xml ${sample_id}.shape
    """
}