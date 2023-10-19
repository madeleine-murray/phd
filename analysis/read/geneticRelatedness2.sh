#!/bin/bash


# The database was copied from pg1: /Gustavo/Datasets/SGDP_americans/SGDP_0_preparation/
#	SGDP.database.pileup.bed
#	SGDP.database.pileup.bim
#	SGDP.database.pileup.fam


# The bam files were copied to this directory otherwise filename did not correspond inside the script
scp ../mergePlatform/CpM0203.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam CpM0203.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
scp ../mergePlatform/CpM0506.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam CpM0506.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
scp ../mergePlatform/CpM11.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam CpM11.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
scp ../mergePlatform/CpM12.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam CpM12.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
scp ../mergePlatform/CpM13.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam CpM13.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
scp ../mergePlatform/CpM34.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam CpM34.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
scp ../mergePlatform/CpM4243.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam CpM4243.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam


# Index bam files
samtools index CpM0203.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
samtools index CpM0506.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
samtools index CpM11.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
samtools index CpM12.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
samtools index CpM13.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
samtools index CpM34.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
samtools index CpM4243.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam


# Create list of samples
rm -f listOfBamInput
for file in $(ls *trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam); do
	echo ${file} >> listOfBamInput
done


# Add ancient Maya from Copan site to the dataset
# (mapping quality and read length is set to zero because already filtered in the previous directory)
python HapSNP_v8_mad.py \
		-data /st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/READ2/SGDP.database \
		-I listOfBamInput \
		-R /Reference_Genomes/hs37d5/hs37d5.fa \
		-p 7 \
		-o Copan \
		-mq 0 \
		-minbp 0 \
		-mem 250 \
		-gatk /Software/GenomeAnalysisTK.jar

#/Programs/HapSNP_v8.py
#	[-h] -data DATASET
#	-I LIST_OF_BAMS
#	[-R REF]
#	[-p PARALLEL]
#	[-o OUTPUT_MERGE_FILE]
#	[-mq MAPPING_QUAL]
#	[-minbp BP_LENGTH]
#	[-mem JAVA_MEM]
#	[-gatk PATH_GATK]
# I changed line 360 from "chrom=col[0]" to "chrom=col[0][3]" to avoid "chr" in chromosomal number


# Merge all Copan individuals
awk '{print $1,$2}' Copan.SGDP.hs37d5.bq30.mq0.0bp.homozyg.fam | head > listOfSamples
plink --bfile Copan.SGDP.hs37d5.bq30.mq0.0bp.homozyg \
	--keep listOfSamples \
	--make-bed \
	--recode \
	--transpose \
	--out Copan-only


# Relationship Estimation from Ancient DNA
python READ.py Copan-only

echo "Simplify manually: Copan.SGDP.hs37d5.bq30.mq0.0bp.homozyg.fam"
echo "From 'CpM*_WGS_UDG_DS' to 'CpM*'"

# Remove bam files to save memory
rm *trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
rm *trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam.bai
