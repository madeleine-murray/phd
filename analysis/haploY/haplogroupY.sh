#!/bin/bash


# Create list of individuals
males=("CpM11 CpM12 CpM13")
rm -f listOfMales
for ind in $males; do
	echo "../mergePlatform/${ind}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam" >> listOfMales
done


# The database was copied from pg1: /Gustavo/Datasets/SGDP_americans/SGDP_0_preparation/
#	SGDP.database.pileup.bed
#	SGDP.database.pileup.bim
#	SGDP.database.pileup.fam
# Correct chromosome format
awk '{print "Y",$2,$3}'             SGDP.chrY.pileup.bed > SGDP.chrY.bed
awk '{print "chrY",$2,$3,$4,$5,$6}' SGDP.chrY.pileup.bim > SGDP.chrY.bim
awk '{print $0}'                    SGDP.chrY.pileup.fam > SGDP.chrY.fam


# Haplogroup information was copied from pg1:
#	/home/niallc/Y_Chromosome/ISOGG_Y_Chromosome_SNPs_filtered_no_indels_organized_no_repeats.csv
# Prepare info file
python /Programs/Y_Finder/info_file_maker.py \
	--info_file ISOGG_Y_Chromosome_SNPs_filtered_no_indels_organized_no_repeats \
	> SGDP.chrY.info

# Determine haplogroup using script: filters are set to zero because bam files were already filtered and it is not logical to re-filter
python /Programs/Y_Finder/Y_Finder_v3.py \
			-data SGDP.chrY \
			-I listOfMales \
			-R /Reference_Genomes/hs37d5/hs37d5.fa \
			-p 3 \
			-o /st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/hapY \
			-mq 0 \
			-minbp 0 \
			-maxbp -1 \
			-bq 0 \
			-mem 50 \
			-gatk /Software/GenomeAnalysisTK.jar \
			-Zcov 10
