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
    container 'biocontainers/bowtie2:v2.4.1_cv1' // run bowtie2 in docker biocontainer

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
    container 'community.wave.seqera.io/library/pear:0.9.6--568591d27ee05b22' // run pear in seqera container
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
    container 'biocontainers/fastqc:v0.11.9_cv8' // run fastqc in docker biocontainer

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
    container 'community.wave.seqera.io/library/fastp:1.0.1--c8b87fe62dcc103c' // run fastp in seqera container
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
    container 'community.wave.seqera.io/library/cutadapt:5.0--991bbd2e184b7014' //run cutadapt in seqera container
    publishDir "${params.outdir}/logs/${sample_id}/", mode: "copy", pattern: "*.log"

    input:
    tuple val(sample_id), val(treatment), path(reads)

    output:
    tuple val(sample_id), val(treatment), path("*.fastq.gz"), emit: reads
    path "*.log"                                            , emit: log

    script:
    """
    cutadapt --cores $task.cpus --adapter ${params.read1_adapter} --front ${params.read2_adapter} \
    --quality-cutoff ${params.quality_cutoff} --minimum-length ${params.min_read_length} \
    --output trimmed_${sample_id}_${treatment}.fastq.gz \
    ${reads} > ${sample_id}_${treatment}_cutadapt.log
    """
}

process BOWTIE_ALIGNMENT {
    tag "${sample_id}_${treatment} on ${ref_name}"
    container 'zpqu/bowtie2-samtools:v2.5.2-v1.19.2' // run bowtie2/samtools in docker container
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
    container 'community.wave.seqera.io/library/qualimap:2.3--c1797c2253925b3a' // run bamqc in seqera container

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
    container 'multiqc/multiqc:pdf-v1.32' // run multiqc in docker container

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
    publishDir "${params.outdir}/${sample_id}/"     , mode: 'copy', pattern: '*.xml'
    publishDir "${params.outdir}/logs/${sample_id}/", mode: 'copy', pattern: '*.log'
    
    input:
    tuple val(sample_id), val(counts_files)

    output:
    tuple val(sample_id), file('*.shape'), emit: shape_file
    tuple val(sample_id), file('*.xml')  , emit: xml_file
    path "*.log"                         , emit: log


    shell:
    if (params.probing_reagent == 'dms') {
            reactive_bases = 'AC'
        } else if(params.probing_reagent == 'shape') {
            reactive_bases = 'ACGT'
        } else {
            exit 1, "--probing_reagent must be 'dms' or 'shape'"
        }
    if (!params.normalisation) {
            norm_option = "--raw"
        } else {
            norm_option = "--norm-method ${params.normalisation}"
        }
    '''
    rf-norm --processors !{task.cpus} --reactive-bases !{reactive_bases} \
    --scoring-method !{params.scoring_method} !{norm_option} \
    --treated  !{counts_files['treated'][0]} --untreated !{counts_files['control'][0]} \
    --output-dir rf_norm > !{sample_id}_rf_norm.log
    
    # create reactivity file(s) in .shape format
    mv rf_norm/*.xml .
    for file in *.xml;
        do xml_to_shape.py $file $(basename --suffix=.xml $file).shape;
    done
    '''
}

process RAW_COUNTS {
    tag "${sample_id}"
    container 'felixlohmeier/pandas:latest' // run in container
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

process RF_COUNT_SUBSAMPLED {
    // process for generating (downsampled) input files for DRACO
    tag "${sample_id}_${treatment}"
    container 'dincarnato/rnaframework:latest' // run in docker container
    publishDir "${params.outdir}/logs/${sample_id}/", mode: 'copy', pattern: '*.log'

    input:
    tuple val(sample_id), val(treatment), file(bam), file(bai)
    path fasta

    output:
    tuple val(sample_id), val(treatment), file("rf_count/*.mm"), emit: mm_file
    tuple val(sample_id), val(treatment), file("rf_count/*.rc"), emit: rc_file
    path "*.log"                                               , emit: log

    shell:
    target_num_reads = params.draco_subsampling
    '''
    # calculate subsampling ratio
    TARGET_NUM_READS=!{target_num_reads}
    if [ $TARGET_NUM_READS -eq -1 ]; then
        RATIO=1
    else
        NUM_MAPPED_READS=$(samtools view --excl-flags 4 !{bam} | wc -l)
        RATIO=$(awk "BEGIN {if ($TARGET_NUM_READS / $NUM_MAPPED_READS > 1) print 1; else print $TARGET_NUM_READS / $NUM_MAPPED_READS}")
    fi

    samtools view --excl-flags 4 --bam --subsample $RATIO --subsample-seed 1 !{bam} > downsampled_!{sample_id}_!{treatment}.bam &&
    samtools index downsampled_!{sample_id}_!{treatment}.bam --output downsampled_!{sample_id}_!{treatment}.bam.bai

    rf-count --working-threads !{task.cpus} \
    --count-mutations --mutation-map \
    --fasta !{fasta} downsampled_!{sample_id}_!{treatment}.bam \
    > !{sample_id}_!{treatment}_rf_count_downsampled.log
    '''
}

process DRACO {
    tag "${sample_id}_${treatment}"
    publishDir "${params.outdir}/${sample_id}/draco/${treatment}/", mode: 'copy', pattern: '*.json'
    publishDir "${params.outdir}/logs/${sample_id}/"              , mode: 'copy', pattern: '*.log'
    
    input:
    tuple val(sample_id), val(treatment), file(mutation_map_file)

    output:
    tuple val(sample_id), val(treatment), file("*.json"), emit: draco_json
    path "*.log"                                        , emit: log


    script:
    if (params.probing_reagent == 'shape') {
        draco_shape = true
    } else {
        draco_shape = false
    }
    """
    draco --processors $task.cpus --mm ${mutation_map_file} ${draco_shape ? '--shape' : ''} \
    --minWindowsOverlap 0.5 --winOffset 5 \
    --minPermutations 10 --maxPermutations 50 \
    --lookaheadEigengaps 1 --softClusteringIters 30 \
    --minBaseCoverage 10 --minFilteredReads 10 \
    --log-level trace \
    --output ${sample_id}_${treatment}_draco.json > ${sample_id}_${treatment}_draco.log
    """
}

process DRACO_JSON_FIX {
    tag "${sample_id}_${treatment}"

    input:
    tuple val(sample_id), val(treatment), file(draco_json_file)

    output:
    tuple val(sample_id), val(treatment), file("fixed_*.json"), emit: draco_json

    script:
    """
    fix_draco_json.py ${draco_json_file}
    """
}

process RF_JSON2RC {
    tag "${sample_id}_${treatment}"
    container 'dincarnato/rnaframework:latest' // run in docker container
    publishDir "${params.outdir}/logs/${sample_id}/"              , mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/${sample_id}/draco/${treatment}/", mode: 'copy', pattern: '*stoichiometries.txt'

    input:
    tuple val(sample_id), val(treatment), file(draco_json_file)
    tuple val(sample_id), val(treatment), file(downsampled_rc_file)

    output:
    tuple val(sample_id), val(treatment), file('rf_json2rc/*.rc')     , emit: draco_rc_file
    tuple val(sample_id), val(treatment), file('*stoichiometries.txt'), emit: stoichiometries
    path "*.log"                                                      , emit: log

    script:
    """
    rf-json2rc --json ${draco_json_file} --rc ${downsampled_rc_file} \
               --median-pre-cov 0 --min-confs 1 \
               --min-overlap-merge 0.35 --extend 20 \
               --cap-mut-freqs 0.1 --ignore-terminal 0.1 --min-corr-merge 0.65 \
               > ${sample_id}_${treatment}_json2rc.log
    
    mv rf_json2rc/stoichiometries.txt ./${sample_id}_${treatment}_draco_stoichiometries.txt
    """
}

process DRACO_RF_NORM {
    tag "${sample_id}_${treatment}"
    container 'dincarnato/rnaframework:latest' // run in docker container
    publishDir "${params.outdir}/${sample_id}/draco/${treatment}/shape_format/", mode: 'copy', pattern: '*.shape'
    publishDir "${params.outdir}/${sample_id}/draco/${treatment}/xml_format/"  , mode: 'copy', pattern: '*.xml'
    publishDir "${params.outdir}/logs/${sample_id}/"                           , mode: 'copy', pattern: '*.log'
    
    input:
    tuple val(sample_id), val(treatment), val(draco_rc_file)

    output:
    tuple val(sample_id), file('*.xml')  , emit: xml_file
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
    if (!params.normalisation) {
            norm_option = "--raw"
        } else {
            norm_option = "--norm-method ${params.normalisation}"
        }        
    '''
    rf-norm --processors !{task.cpus} --reactive-bases !{reactive_bases} \
    --scoring-method 4 !{norm_option} \
    --treated !{draco_rc_file} \
    --output-dir rf_norm > !{sample_id}_draco_rf_norm.log

    mv rf_norm/*.xml .
    for file in *.xml;
        do xml_to_shape.py $file $(basename --suffix=.xml $file).shape;
    done
    '''
}