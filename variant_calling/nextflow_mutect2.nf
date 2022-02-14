#!/usr/bin/env nextflow


params.outdir = "vcfs/"
params.tmpdir = "tmp/"


ref_chr="/nfs/users/nfs_p/pq1/scratch/ukbiobank/200k_exomes/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"
ref_pon="/nfs/users/nfs_p/pq1/scratch/reference/human/1000g_pon_intervals_ch.hg38.vcf.gz"
ref_germ="/nfs/users/nfs_p/pq1/scratch/reference/human/af-only-gnomad_intervals_ch.hg38.vcf.gz"
intervals="/nfs/users/nfs_p/pq1/scratch/ukbiobank/200k_exomes/resources/genes_CH_xgen.intervals"
tmp_dir="/lustre/scratch119/realdata/mdt1/team163/pq1/tmp"



Channel
	.fromPath(params.samples)
    .splitCsv(header:false)
    .map {row -> row[0]}
	.set{ rowsamples }


process variantCalling {

    tag "${sample}"
    echo true

	queue 'small'
	memory = '2G'
	time '10m' //check
	cpus 1
	maxForks 200
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	maxRetries 3

	publishDir "${params.outdir}", pattern: "*.filt.vcf.gz"

	input:
	val(sample) from rowsamples

	output:
	tuple sample, file("*.gz"), file("*tbi")

	script:
	"""
    cramfile=`cat /nfs/users/nfs_p/pq1/scratch/ukbiobank/200k_exomes/all_cram_files | grep ${sample}`

    # Mutect2
    gatk --java-options "-Xmx2g" Mutect2 -R ${ref_chr} -I \${cramfile} -O ${sample}.unfiltered.vcf.gz --tmp-dir ${tmp_dir} -L ${intervals} --panel-of-normals ${ref_pon} --germline-resource ${ref_germ} --f1r2-tar-gz ${sample}.f1r2.tar.gz

    gatk --java-options "-Xmx2g" LearnReadOrientationModel -I ${sample}.f1r2.tar.gz -O ${sample}.read-orientation-model.tar.gz

    gatk --java-options "-Xmx2g" FilterMutectCalls -R ${ref_chr} -V ${sample}.unfiltered.vcf.gz -O ${sample}.filt.vcf.gz -L ${intervals} --ob-priors ${sample}.read-orientation-model.tar.gz
	
	"""
}


workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}