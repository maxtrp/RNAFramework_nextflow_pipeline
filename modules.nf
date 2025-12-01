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

process PEAR {
    tag "${sample_id}_${treatment}"
    publishDir "${params.outdir}/logs/${sample_id}/", mode: "copy", pattern: "*.log"

    input:
    tuple val(sample_id), val(treatment), path(reads)

    output:
    tuple val(sample_id), val(treatment), path("*.assembled.fastq.gz"), emit: fastq
    path "*.log"                                                      , emit: log

    script:
    """
    pear --forward-fastq ${reads[0]} --reverse-fastq ${reads[1]} --output ${sample_id}_${treatment} \
    --max-uncalled-base 0 --threads $task.cpus > ${sample_id}_${treatment}_pear.log

    cat ${sample_id}_${treatment}.discarded.fastq >> ${sample_id}_${treatment}.assembled.fastq

    gzip ${sample_id}_${treatment}.assembled.fastq
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

process FASTP_DEDUPLICATION {

    tag "${sample_id}_${treatment}"
    publishDir "${params.outdir}/logs/${sample_id}/", mode: "copy", pattern: "*.log"

    input:
    tuple val(sample_id), val(treatment), path(reads)

    output:
    tuple val(sample_id), val(treatment), path("*.fastq.gz"), emit: reads
    path "*.log"                                            , emit: log

    script:
    """
    fastp --thread $task.cpus --dedup --disable_adapter_trimming \
    --in1 ${reads} --out1 deduplicated_${sample_id}_${treatment}.fastq.gz \
    2> ${sample_id}_${treatment}_fastp_dedup.log
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
    cutadapt --cores $task.cpus --adapter ${params.read1_adapter} --adapter ${params.read2_adapter} \
    --quality-cutoff ${params.quality_cutoff} --minimum-length ${params.min_read_length} \
    --output trimmed_${sample_id}_${treatment}.fastq.gz \
    ${reads} > ${sample_id}_${treatment}_cutadapt.log
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
    bowtie2 --threads $task.cpus --norc -x ${ref_name}_index -U ${reads} \
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
    tuple val(sample_id), val(treatment), file("rf_count/*.rc"), file("*_counts.txt"), emit: counts_files
    path "*.log"                                                                     , emit: log

    script:
    """
    rf-count --working-threads $task.cpus \
    --count-mutations --fasta ${fasta} ${bam} \
    > ${sample_id}_${treatment}_rf_count.log

    rf-rctools view rf_count/*.rc > ${sample_id}_counts.txt
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


    shell:
    if (params.probing_reagent == 'dms') {
            reactive_bases = 'AC'
        } else if(params.probing_reagent == 'shape') {
            reactive_bases = 'ACGT'
        } else {
            exit 1, "--probing_reagent must be 'dms' or 'shape'"
        }
    '''
    rf-norm --processors !{task.cpus} --scoring-method 3 --raw --reactive-bases !{reactive_bases} \
    --treated  !{counts_files['treated'][0]} --untreated !{counts_files['control'][0]} \
    --output-dir rf_norm > !{sample_id}_rf_norm.log
    
    # create reactivity file in .shape format
    for file in rf_norm/*.xml;
        do xml_to_shape.py $file $(basename --suffix=.xml $file).shape;
    done
    '''
}

process RAW_COUNTS {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/", mode: 'copy', pattern: '*.csv'
    
    input:
    tuple val(sample_id), val(counts_files)

    output:
    tuple val(sample_id), file('*.csv')


    script:
    """
    raw_counts.py ${counts_files['treated'][1]} ${counts_files['control'][1]}
    """
}

process RF_COUNT_MUTATION_MAP_SUBSAMPLED {
    tag "${sample_id}_${treatment}"
    container 'dincarnato/rnaframework:latest' // run in docker container
    publishDir "${params.outdir}/logs/${sample_id}/", mode: 'copy', pattern: '*.log'

    input:
    tuple val(sample_id), val(treatment), file(bam), file(bai)
    path fasta

    output:
    tuple val(sample_id), val(treatment), file("rf_count/*.mm"), emit: mm_file
    path "*.log"                                               , emit: log

    shell:
    target_num_pairs = params.draco_subsampling
    '''
    # calculate subsampling ratio
    TARGET_NUM_PAIRS=!{target_num_pairs}
    NUM_MAPPED_PAIRS=$(samtools view --excl-flags 132 !{bam} | wc -l)
    RATIO=$(awk "BEGIN {if ($TARGET_NUM_PAIRS / $NUM_MAPPED_PAIRS > 1) print 1; else print $TARGET_NUM_PAIRS / $NUM_MAPPED_PAIRS}")

    samtools view --excl-flags 12 --bam --subsample $RATIO --subsample-seed 1 !{bam} > downsampled_!{sample_id}_!{treatment}.bam &&
    samtools index downsampled_!{sample_id}_!{treatment}.bam --output downsampled_!{sample_id}_!{treatment}.bam.bai

    rf-count --working-threads !{task.cpus} \
    --count-mutations --mutation-map \
    --fasta !{fasta} downsampled_!{sample_id}_!{treatment}.bam \
    > !{sample_id}_!{treatment}_rf_count_downsampled.log
    '''
}

process DRACO {
    tag "${sample_id}_${treatment}"
    publishDir "${params.outdir}/${sample_id}/", mode: 'copy', pattern: '*.json'
    publishDir "${params.outdir}/logs/${sample_id}/", mode: 'copy', pattern: '*.log'
    
    input:
    tuple val(sample_id), val(treatment), file(mutation_map_file)

    output:
    path "*.json", emit: draco_json
    path "*.log" , emit: log


    script:
    if (params.probing_reagent == 'shape') {
        draco_shape = true
    } else {
        draco_shape = false
    }
    """
    draco --mm ${mutation_map_file} ${draco_shape ? '--shape' : ''} \
    --absWinLen 100 --absWinOffset 5 --minPermutations 10 --maxPermutations 50 --firstEigengapShift 0.95 \
    --lookaheadEigengaps 1 --softClusteringIters 30 --softClusteringInits 500 --softClusteringWeightModule 0.005 \
    --output ${sample_id}_${treatment}_draco.json > ${sample_id}_${treatment}_draco.log
    """
}