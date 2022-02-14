#!/usr/bin/env bash

# Run nextflow workflow in batches

for i in {01..50}
do
	cd /nfs/users/nfs_p/pq1/scratch/ukbiobank/200k_exomes/analysis
    echo 'copying'
	# copy nextflow to each folder
	cp nextflow_mutect2.nf batch_$i
    cd batch_$i
	echo 'bsub'
	bsub -J batch${i} -q long -R'select[mem>2000] rusage[mem=2000] span[hosts=1]' -M 2000 -n 2 -o bsub.o -e bsub.e "nextflow run nextflow_mutect2.nf --samples samples_batch_${i}"
done