#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "/labs/sbio/gpops/ADNmt/FASTQ/*{R1,R2}_001.fastq.gz"
params.ref = "/labs/sbio/mito-forense/ref/rCRS.fasta"
params.resultsDir = "results"
params.pl = 'Illumina'
params.pm = 'NextSeq'
params.tmpdir = "$baseDir/tmp" 
params.gatk = "gatk"
params.picard = "java -jar /labs/sbio/dpservs/programasd/progra_vc/picard.jar"
params.mutserve = 'java -jar /home/programs/mutserve.jar'
params.haplogrep = 'java -jar /home/programs/haplogrep.jar'


include { fastqc;
	 trim_Galore;
	 multiqc_fastqc;
	 multiqc_trim;
	 align;
	 mergeSam;
	 markDuplicates;
	 getMetrics;
	 mutServe;
	 haploGrep } from './modules.nf'


workflow {
    read_pairs_ch = Channel.fromFilePairs( params.reads )
    fastqc(read_pairs_ch)
    trim_Galore(read_pairs_ch)
    fastqc.out.fq_files \
	| collect \
	| multiqc_fastqc  
    trim_Galore.out.trim_path \
        | collect \
        | multiqc_trim
    align(trim_Galore.out.trim_fq)
    align.out.aligned_path \
        | map { file ->
           def key = file.name.toString().tokenize('_').get(1)
         return tuple(key, file)
          } \
        | groupTuple() \
        | mergeSam

   markDuplicates(mergeSam.out.merged_sam_ch)
   getMetrics(markDuplicates.out.bam_for_variant_calling)
   mutServe(markDuplicates.out.bam_for_variant_calling)
   haploGrep(mutServe.out.vc_vcf)
}

