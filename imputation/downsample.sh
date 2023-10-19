#!/bin/bash
# Downsample bam files
# Please calculate the downsample ratio accordingly and change parameters

DOWNSAMPLING=("0X01 0X05 0X1 0X5 1X0 2X0 3X0 4X0 5X0")
downsampling[0]="0X01"
downsampling[1]="0X05"
downsampling[2]="0X1"
downsampling[3]="0X5"
downsampling[4]="1X0"
downsampling[5]="2X0"
downsampling[6]="3X0"
downsampling[7]="4X0"
downsampling[8]="5X0"
proportion[0]=0.00046
proportion[1]=0.00230
proportion[2]=0.00459
proportion[3]=0.02296
proportion[4]=0.04592
proportion[5]=0.09184
proportion[6]=0.13777
proportion[7]=0.18369
proportion[8]=0.22961


WD="/home/madeleine/imputationValidation"
sampleID="Yamnaya"
sampleDetail="sorted.grouped.mq25.rmdup.indel.softclip.rl34"
mkdir -p ${WD}/bam/
mkdir -p ${WD}/qualimap/


# 1_Domwnsample to match other samples of lower coverage

for (( i=0; i<"${#downsampling[@]}"; i++))
do
	java -jar /research/picard-tools-2.4.1/picard.jar DownsampleSam \
	I=/Gustavo/Downloaded_Data/deBarrosDamgaard_2018_Horse_Herders/Alignments/hg19_Cristina_with_rCRS/Bam_Files/${sampleID}.${sampleDetail}.bam \
	O=${WD}/bam/${sampleID}.${sampleDetail}.${downsampling[i]}.bam \
	RANDOM_SEED=1234 \
	P=${proportion[i]}
done


# 2_Check quality with qualimap

for downsample in ${DOWNSAMPLING}; do

	/home/niallc/bin/qualimap_v2.2.1/qualimap bamqc \
		-bam ${WD}/bam/${sampleID}.${sampleDetail}.${downsample}.bam -c \
		> ${WD}/qualimap/qualimap.${sampleID}.${downsample}.picard
done


# 3_Convert bam to vcf

for downsample in ${DOWNSAMPLING}; do

	# replace all identifiers
	java -Xmx5g -jar /The_Luggage/Home_Directory_Backup_May2016/Programming/AddOrReplaceReadGroups.jar \
		I=${WD}/bam/${sampleID}.${sampleDetail}.${downsample}.bam \
		O=${WD}/bam/${sampleID}.${sampleDetail}.${downsample}.REGROUPED.bam \
		RGLB=${sampleID} RGSM=${sampleID} RGPU=Illumina RGPL=ILLUMINA RGID=${sampleID} \
		RGDS=${sampleID}.${sampleDetail}.${downsample}

	# create index
	samtools index ${WD}/bam/${sampleID}.${sampleDetail}.${downsample}.REGROUPED.bam

done
