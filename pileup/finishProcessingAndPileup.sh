#!/bin/bash

parallel=15 # number of threads to use for softclipping
date=18072022


# The database was copied from pg1: /Gustavo/Datasets/SGDP_americans/SGDP_0_preparation/
mkdir -p pileupFiles
if [ ! -f SGDP.database.pileup.bed ] || [ ! -f SGDP.database.pileup.bim ] || [ ! -f SGDP.database.pileup.fam ]; then
        scp madeleine@134.226.153.112://Gustavo/Datasets/SGDP_americans/SGDP_0_preparation/SGDP.database.pileup.bed SGDP.database.pileup.bed
        scp madeleine@134.226.153.112://Gustavo/Datasets/SGDP_americans/SGDP_0_preparation/SGDP.database.pileup.bim SGDP.database.pileup.bim
        scp madeleine@134.226.153.112://Gustavo/Datasets/SGDP_americans/SGDP_0_preparation/SGDP.database.pileup.fam SGDP.database.pileup.fam
fi
# NB: It only includes modern SGDP individuals



############################################################################################################
# PART 1:finish processing #


# INDEL #
for file in $(ls ../*/alignments/hs37d5/bam_files/*sorted.grouped.duprm.bam); do
	
	# Obtain directory and sample name
        pattern="bam_files"              # recognised regex
        path=${file%%$pattern*}          # path until hs37d5
        prefix=`echo ${#path} + 10 | bc` # numerical start of bam_files in path
        sample=${file:$prefix:-25}       # sample name with trimming name
	
	# Indel realignment (if required)
	if [[ ! -f ${path}final_files/${sample}.sorted.grouped.duprm.indel.bam ]]; then
		samtools index ${path}bam_files/${sample}.sorted.grouped.duprm.bam
        	java -jar /Software/GenomeAnalysisTK.jar \
                	-T RealignerTargetCreator \
                	-R /Reference_Genomes/hs37d5/hs37d5.fa \
			-I ${path}bam_files/${sample}.sorted.grouped.duprm.bam \
                	-o ${path}final_files/${sample}.sorted.grouped.duprm.indel.intervals
			# -allowPotentiallyMisencodedQuals
        	java -jar /Software/GenomeAnalysisTK.jar \
                	-T IndelRealigner \
                	-R /Reference_Genomes/hs37d5/hs37d5.fa \
                	-I ${path}bam_files/${sample}.sorted.grouped.duprm.bam \
                	-o ${path}final_files/${sample}.sorted.grouped.duprm.indel.bam \
                	-targetIntervals ${path}final_files/${sample}.sorted.grouped.duprm.indel.intervals
			# -allowPotentiallyMisencodedQuals
	fi
done




# SOFTCLIPPING #
rm -f listOfFiles # start a fresh file
for file in $(ls ../*/alignments/hs37d5/bam_files/*sorted.grouped.duprm.bam); do

        # Obtain directory and sample name
        pattern="bam_files"              # recognised regex
        path=${file%%$pattern*}          # path until hs37d5
        prefix=`echo ${#path} + 10 | bc` # numerical start of bam_files in path
        sample=${file:$prefix:-25}       # sample name with trimming name

        # Create a list of samples to map filter and softclip
	if [[ ! -f ${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam ]]; then
        	echo "${path}final_files/${sample}.sorted.grouped.duprm.indel.bam" >> listOfFiles
	fi
done

# Filter and softclip (if required)
if [[ -f listOfFiles ]]; then
	python /Programs/filter_bam.py -mq 25 \
        	                        -minbp 34 \
                	                -clip 2 \
                        	        -p ${parallel} \
                                	-i listOfFiles
	rm listOfFiles
fi
############################################################################################################


############################################################################################################
# Part 2: Pile-up with SGDP


for file in $(ls ../*/alignments/hs37d5/bam_files/*sorted.grouped.duprm.bam); do

        # Obtain directory and sample name
        pattern="bam_files"              # recognised regex
        path=${file%%$pattern*}          # path until hs37d5
        prefix=`echo ${#path} + 10 | bc` # numerical start of bam_files in path
        sample=${file:$prefix:-25}       # sample name with trimming name
	sampleInd=${sample%%.*}          # sample name without trimming name

        # QUALIMAP (if required)
        if [ ! -f ${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp_stats/qualimapReport.html ]; then
                /Software/qualimap_v2.1.1/qualimap bamqc --java-mem-size=50G \
                        -bam ${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam
        fi

	# Make sure they are indexed
	if [ ! -f ${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam.bai ]; then
		samtools index ${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam
	fi

	# Moves files to working directory otherwise filename in script will not correspond
        if [ -f ${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam ]; then
                mv ${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam \
			${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam
		mv ${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam.bai \
			${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam.bai
		echo "${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam" >> listOfBams
        fi

	# Modify name of previous piled up file to match the automatic output, so that files that already exist are not re-created (same time)
        [ ! -f pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.bed ] || \
		mv pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.bed ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.bed;
	[ ! -f pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.bim ] || \
		mv pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.bim ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.bim;
	[ ! -f pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.fam ] || \
		mv pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.fam ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.fam;
	[ ! -f pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.log ] || \
		mv pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.log ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.log;
	[ ! -f pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.map ] || \
		mv pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.map ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.map;
	[ ! -f pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.ped ] || \
		mv pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.ped ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.ped;
	[ ! -f pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.nosex ] || \
		mv pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.nosex ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.nosex;
done


# Ouput dataset with published ancient
# (mapping quality and read length is set to zero because already filtered with filter_bam.py
# and I want bam files and pile up to have the same filtering (there is small difference if you re-filter))
python /Programs/HapSNP_v9.py \
		-data /external_storage/madeleine/Downloaded_Data/shotgun/Pileup_SGDP/SGDP.database \
                -I listOfBams \
                -R /Reference_Genomes/hs37d5/hs37d5.fa \
                -p ${parallel} \
                -o PublishedShotgun \
                -mq 0 \
                -minbp 0 \
                -mem 250 \
                -gatk /Software/GenomeAnalysisTK.jar


rm -f listOfPileups
rm -f pileup_${date}.log
for file in $(ls ../*/alignments/hs37d5/bam_files/*sorted.grouped.duprm.bam); do

        # Obtain directory and sample name
        pattern="bam_files"              # recognised regex
        path=${file%%$pattern*}          # path until hs37d5
        prefix=`echo ${#path} + 10 | bc` # numerical start of bam_files in path
        sample=${file:$prefix:-25}       # sample name with trimming name
	sampleInd=${sample%%.*}          # sample name without trimming name

	# Modify pile up name for a better description
        [ ! -f ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.bed ] || \
		mv ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.bed pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.bed
        [ ! -f ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.bim ] || \
		mv ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.bim pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.bim
        [ ! -f ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.fam ] || \
		mv ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.fam pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.fam
        [ ! -f ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.log ] || \
		mv ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.log pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.log
        [ ! -f ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.map ] || \
		mv ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.map pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.map
        [ ! -f ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.ped ] || \
		mv ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.ped pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.ped
        [ ! -f ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.nosex ] || \
		mv ${sample}.SGDP.hs37d5.bq30.mq0.0bp.homozyg.nosex pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.nosex

        if [ -f ${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam ]; then
		if [ -f pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg.fam ]; then
        	        echo "pileupFiles/${sampleInd}.SGDP.hs37d5.bq30.mq25.clip2bp.34bp.homozyg" >> listOfPileups      # prepare a list to merge all pileup files
        	else
                	echo "${sampleInd} has not been included in final pileup up"               >> pileup_${date}.log # keep track of potential problems
        	fi
		# Move files back to processing directory
                mv ${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam \
			${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam
                mv ${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam.bai \
			${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam.bai
        fi
done

# Merge all files into one
if [ -f listOfPileups ]; then
        first=`head -1 listOfPileups`                                                                 # get the name of the first file
        sed 1,1d listOfPileups > listOfPileup                                                         # filter out first file (otherwie will be merged twice)
        plink1.9 --bfile ${first} --merge-list listOfPileup --make-bed --out publishedShotgun_${date} # finish with merging
        rm listOfPileup                                                                               # clean up directory
fi
############################################################################################################
