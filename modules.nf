#!/usr/bin/env nextflow

process fastqc {
    cache 'lenient'
    publishDir params.resultsDir ,  mode: 'copy'
 
    input:
    tuple val(sample), path(reads)

    output:
    path("fastqc/fq_${sample}/*"), emit: fq_files

 
    script:
    """
    mkdir -p fastqc/fq_${sample}
    fastqc -o fastqc/fq_${sample} -f fastq -q ${reads[0]} ${reads[1]}

    """
}


process trim_Galore {
  cache 'lenient'
  publishDir params.resultsDir, mode: 'symlink'

  input:
  tuple val(sample), path(reads)

  output:
  tuple val(sample), path("trimming_files/tg_${sample}/*.fq.gz")      , emit: trim_fq
  tuple val(sample), path("trimming_files/tg_${sample}/*report.txt")          , emit: trim_report
  tuple val(sample), path("trimming_files/tg_${sample}/*_fastqc.html"), emit: trim_html
  tuple val(sample), path("trimming_files/tg_${sample}/*_fastqc.zip") , emit: trim_zip
  path "trimming_files/tg_${sample}/*", emit: trim_path
  
  script:
  """
    mkdir -p trimming_files
    trim_galore -o trimming_files/tg_${sample} --paired --fastqc ${reads[0]} ${reads[1]}
  """
}


process multiqc_fastqc {
  cache 'lenient'
  publishDir params.resultsDir, mode: 'symlink'

  input:
  path (fastqc)

  output:
  path("multiqc/fastqc/*")   , emit: multiqc_fq

  script:

  """
    echo -e "\n"

    mkdir -p multiqc
    multiqc ${fastqc} -o multiqc/fastqc

  """
}


process multiqc_trim {
  cache 'lenient'
  publishDir params.resultsDir, mode: 'symlink'

  input:
  path (trimg)

  output:
  path("multiqc/trimg/*")   , emit: multiqc_trimg

  script:

  """
    echo -e "\n"

    mkdir -p multiqc
    multiqc ${trimg} -o multiqc/trimg

  """
}


process align {
    cache 'lenient'
    publishDir params.resultsDir, mode:'symlink'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("aligned_reads/${sample}_aligned_reads.sam"),     emit: aligned_reads_ch
    path 'aligned_reads/*sam',	emit: aligned_path

    script:
    """
    mkdir -p aligned_reads

    FLOWCELL=\$(zcat ${reads[0]} | head -n 1 | awk -F ':' '{ print \$3 }')
    LANE=\$(zcat ${reads[0]} | head -n 1 | awk -F ':' '{ print \$4 }')
   
    RGID=\$FLOWCELL'.'\$LANE'.'${sample}
    ReadGroup="@RG\\tID:\${RGID}\\tLB:${sample}\\tPL:${params.pl}\\tPM:${params.pm}\\tSM:${sample}"
    
    bwa mem \
        -K 100000000 \
        -v 3 \
        -t ${task.cpus} \
        -Y \
        -R \"\${ReadGroup}\" \
        ${params.ref} \
        ${reads[0]} \
        ${reads[1]} \
        > aligned_reads/${sample}_aligned_reads.sam
    """
}

process mergeSam {
   
   cache 'lenient'
   publishDir params.resultsDir, mode:'symlink'
 
   input:
   tuple val(key), file(bam_files)

   output:
   tuple val(key), path("merged_sam/${key}_merged.sam"),     emit: merged_sam_ch

   script:
   """
   mkdir -p merged_sam
   ${params.picard} \
            MergeSamFiles \
            ${'INPUT='+bam_files.join(' INPUT=')} \
            SORT_ORDER='coordinate' \
	    OUTPUT=merged_sam/${key}_merged.sam
   """
}

process markDuplicatesSpark {
    cache 'lenient'
    publishDir params.resultsDir, mode:'symlink'

    input:
    tuple val(sample), path(merged_reads)

    output:
    tuple val(sample), path("dedup_sorted/${sample}_sorted_dedup.bam"),    emit: bam_for_variant_calling
    tuple val(sample), path("dedup_sorted/${sample}_dedup_metrics.txt"),   emit: dedup_qc_ch
    tuple val(sample), path("dedup_sorted/${sample}_sorted_dedup.bam.bai"),   emit: dedup_index_qc_ch

    script:
    """
    mkdir -p ${params.tmpdir}/${workflow.runName}/${sample}
    mkdir -p dedup_sorted
    ${params.gatk} --java-options "-Djava.io.tmpdir=${params.tmpdir}/${workflow.runName}/${sample}" \
         MarkDuplicatesSpark \
        -I ${sample}_merged.sam \
        -M dedup_sorted/${sample}_dedup_metrics.txt \
        -O dedup_sorted/${sample}_sorted_dedup.bam \
	--create-output-bam-index 

    rm -r ${params.tmpdir}/${workflow.runName}/${sample}
    """
}

process getMetrics {
    cache 'lenient'
    publishDir params.resultsDir, mode:'symlink'

    input:
    tuple val(sample), path(sorted_dedup_reads)

    output:
    tuple val(sample),
          path("metrics/${sample}_alignment_metrics.txt"), \
          path("metrics/${sample}_insert_metrics.txt"), \
          path("metrics/${sample}_insert_size_histogram.pdf"), \
          path("metrics/${sample}_depth_out.txt"),             emit: metrics_qc_ch

    script:
    """
    mkdir -p metrics
    picard-tools \
        CollectAlignmentSummaryMetrics \
        R=${params.ref} \
        I=${sorted_dedup_reads} \
        O= metrics/${sample}_alignment_metrics.txt
    picard-tools \
        CollectInsertSizeMetrics \
        INPUT=${sorted_dedup_reads} \
        OUTPUT= metrics/${sample}_insert_metrics.txt \
        HISTOGRAM_FILE= metrics/${sample}_insert_size_histogram.pdf
    samtools coverage -w 32 ${sorted_dedup_reads} > metrics/${sample}_depth_out.txt
    """
}

process mutServe {
    cache 'lenient'
    publishDir params.resultsDir, mode:'symlink'

    input:
    tuple val(sample), path(dedup_reads)

    output:
    tuple val(sample), path("vc/${sample}.txt"),    emit: vc_tabular
    tuple val(sample), path("vc/${sample}.vcf.gz"),   emit: vc_vcf

    script:
    """
    mkdir -p vc
    ${params.mutserve} \
        call \
	--reference ${params.ref} \
	--output "vc/${sample}.vcf.gz" \
	${sample}_sorted_dedup.bam
    """
}


process haploGrep {
    cache 'lenient'
    publishDir params.resultsDir, mode:'symlink'

    input:
    tuple val(sample), path(vc_vcf)

    output:
    tuple val(sample), path("haplogr/${sample}.hg.txt"),    emit: haplogroup_txt

    script:
    """
    mkdir -p haplogr
    ${params.haplogrep} \
        classify \
        --in ${sample}.vcf.gz \
        --format vcf \
        --out "haplogr/${sample}.hg.txt" 
    """
}


