#!/bin/bash

# Determine haplogroup using script: filters are set to zero because bam files were already filtered and it is not logical to re-filter
individuals=("CpM0203 CpM0506 CpM11 CpM12 CpM13 CpM34 CpM4243")

for ind in ${individuals}; do
	
	# Extract mtDNA only
        samtools view -b ../mergePlatform/${ind}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam MT \
                > ${ind}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.mt.bam
        samtools index ${ind}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.mt.bam

	# Create list for input
        echo "${ind}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.mt.bam" > listIndividuals
	
	# Run haplogrep and all the other stuff
        python /Programs/MitoMaster_v5.py \
                        -R /Reference_Genomes/Human_Mitochondria_revised_Cambridge/Human_mitochondria_rCRS.fasta \
                        -n 14 \
                        -p 7 \
                        -mq 0 \
                        -bq 0 \
                        -dp 50000 \
                        -mem 50 \
                        -I listIndividuals \
                        -o mitochondrial_${ind} \
                        -minbp 0 \
                        -clip 0 \
                        -s "" \
                        -indel Y \
                        -remap Y
done
rm listIndividuals
