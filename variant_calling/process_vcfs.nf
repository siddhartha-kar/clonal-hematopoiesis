#!/usr/bin/env nextflow

params.outdir = "results/"

ref_chr="/nfs/users/nfs_p/pq1/scratch/ukbiobank/200k_exomes/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"
vep_cache="/nfs/users/nfs_p/pq1/scratch/reference/vep/"
ch_genes="/nfs/users/nfs_p/pq1/scratch/ukbiobank/genes"

Channel
 	.from(01..50)
    .map { (it < 10 ? "0" : "") + it }
 	.set{ folders }


process processVCFs {

    tag "${folder}"
    echo true

	queue 'small'
	memory = '1G'
	cpus 1
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	maxRetries 3

	input:
	val(folder) from folders

	output:
	file("sel_vars_batch_${folder}") into selvars
    file("other_vars_batch_${folder}") into othervars

	script:
	"""
    find -L $launchDir/batch_$folder/vcfs -iname "*.filt.vcf.gz" > list_vcfs_${folder}.txt

	while IFS= read -r line;

    do
    
    filename=`echo \$line | awk -F"/" '{print \$13}' | awk -F"." '{print \$1}'`

    #Mutations that PASS FilterMutectCalls
    bcftools norm -Ov -m-any -f $ref_chr \$line | bcftools query -i 'FMT/DP >=7 & FMT/AD[:1] >= 2 & FMT/SB[:2]>0 & FMT/SB[:3]>0 & FILTER = "PASS"' -f '[%SAMPLE] %CHROM %POS %REF %ALT %FILTER [ %MBQ{1} %MMQ{1} %AD{0} %AD{1}]\n' | awk -v VAR="\$filename" '{print VAR,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$10/(\$9+\$10)}' >> sel_vars_batch_${folder}

    # Get the other variants with other flags
    bcftools norm -Ov -m-any -f $ref_chr \$line | bcftools query -i 'FMT/DP >=7 & FMT/AD[:1] >= 2 & FMT/SB[:2]>0 & FMT/SB[:3]>0 & FILTER = "germline"' -f '[%SAMPLE] %CHROM %POS %REF %ALT %FILTER [ %MBQ{1} %MMQ{1} %AD{0} %AD{1}]\n' | awk -v VAR="\$filename" '{print VAR,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$10/(\$9+\$10)}' >> other_vars_batch_${folder}

    done < list_vcfs_${folder}.txt
    
    """

}



process selectVars {

	publishDir "${params.outdir}/" , pattern: "all_sel_vars.txt"

	queue 'small'
	memory = '1G'
	cpus 1
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	maxRetries 3

    input:
	file(sel) from selvars.collect()
	file(other) from othervars.collect()

	output:
	file("all_sel_vars.txt") into selectedvars
	file("spdi_vars.txt") into spdivars

	script:
	"""
    cat $sel > sel_vars.txt 
	
	cat sel_vars.txt | cut -d' ' -f3-6 | sort | uniq -c | sort -nr | awk '{if(\$1 >= 2){print \$2,\$3,\$4,\$5}}' > uniq_vars_PASS_n2.txt
    
	#Rescue vars flagged as germline that passed in other samples
	cat $other | grep -E '(\\s)germline(\\s)' > other_vars_germ.txt
	grep -w -f uniq_vars_PASS_n2.txt other_vars_germ.txt > other_vars_germ_selected.txt

	cat sel_vars.txt other_vars_germ_selected.txt > all_sel_vars.txt

	## generate a SPDI file to use vep
	awk 'BEGIN {FS=" ";OFS=":"} {print \$3,\$4-1,\$5,\$6}' all_sel_vars.txt | sort -V | uniq | awk 'gsub(/*/,"",\$1)1' > spdi_vars.txt
	
	"""

}




process runVEP {

	publishDir "${params.outdir}/" , pattern: "all_vars_vep.txt"

	queue 'normal'
	memory = '4G'
	cpus 4
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	maxRetries 3

	input:
	file(spdi) from spdivars

	output:
	file("all_vars_vep.txt") into res

	"""
	#VEP
	vep --cache --merged --dir_cache $vep_cache --fork 4 --input_file $spdi --output_file tmp_vep.txt --everything --exclude_predicted --total_length --show_ref_allele --check_existing --exclude_null_alleles --clin_sig_allele 1 --flag_pick --force_overwrite --no_stats --tab

	#Filter info
	filter_vep -i tmp_vep.txt -f "SYMBOL in $ch_genes" | grep -v "^##" > all_vars_vep.txt
	"""

}
