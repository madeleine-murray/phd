#!/bin/bash

parallel=15  # number of threads to use for softclipping
date=07122022



# The database was copied from pg1: /Gustavo/Datasets/SGDP_americans/SGDP_0_preparation/
if [ ! -f HOA.database.pileup.bed ] || [ ! -f HOA.database.pileup.bim ] || [ ! -f HOA.database.pileup.fam ]; then
	scp madeleine@134.226.153.112://home/niallc/Datasets/Lazaridis/Lazaridis.rmsnpids.90sRm.letters.Mathmissnprm.nosex.names_fixed/Lazaridis.rmsnpids.90sRm.letters.Mathmissnprm.nosex.names_fixed.bim HOA.database.pileup.bim
        scp madeleine@134.226.153.112://home/niallc/Datasets/Lazaridis/Lazaridis.rmsnpids.90sRm.letters.Mathmissnprm.nosex.names_fixed/Lazaridis.rmsnpids.90sRm.letters.Mathmissnprm.nosex.names_fixed.fam HOA.database.pileup.fam
	cat HOA.database.pileup.bim | awk '{print $1,$4-1,$4}' > HOA.database.pileup.bed # calculate starting SNP coordinates from end coordinates
fi


rm -f listOfBams
for file in $(ls ../*/alignments/hs37d5/final_files/*sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam); do

        # Obtain directory and sample name
        pattern="final_files"            # recognised regex
        path=${file%%$pattern*}          # path until hs37d5
        prefix=`echo ${#path} + 12 | bc` # numerical start of bam_files in path
        sample=${file:$prefix:-49}       # sample name with trimming name
        sampleInd=${sample%%.*}          # sample name without trimming name

        # Moves files to working directory otherwise filename in script will not correspond
        if [ -f ${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam ]; then
        	if [ ! -f ${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam.bai ]; then
                	samtools index ${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam
        	fi # Make sure files are index
                mv ${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam \
                        ${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam
                mv ${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam.bai \
                        ${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam.bai
		
		echo "${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam" >> listOfBams
	fi

	# Change pileup file names back to default to avoid re-piling up everything
	if [ -f pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg.fam ]; then
                mv pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg.bed ${sampleInd}.HOA.hs37d5.bq30.mq0.0bp.homozyg.bed
                mv pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg.bim ${sampleInd}.HOA.hs37d5.bq30.mq0.0bp.homozyg.bim
                mv pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg.fam ${sampleInd}.HOA.hs37d5.bq30.mq0.0bp.homozyg.fam
                mv pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg.log ${sampleInd}.HOA.hs37d5.bq30.mq0.0bp.homozyg.log
                mv pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg.map ${sampleInd}.HOA.hs37d5.bq30.mq0.0bp.homozyg.map
                mv pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg.nosex ${sampleInd}.HOA.hs37d5.bq30.mq0.0bp.homozyg.nosex
                mv pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg.ped ${sampleInd}.HOA.hs37d5.bq30.mq0.0bp.homozyg.ped
	fi
done


# Ouput dataset with published ancient
# (mapping quality and read length is set to zero because already filtered with filter_bam.py
# and I want bam files and pile up to have the same filtering (there is small difference if you re-filter))
python /Programs/HapSNP_v9.py \
                -data /external_storage/madeleine/Downloaded_Data/shotgun/Pileup_HOA/HOA.database \
                -I listOfBams \
                -R /Reference_Genomes/hs37d5/hs37d5.fa \
                -p ${parallel} \
                -o PublishedShotgun \
                -mq 0 \
                -minbp 0 \
                -mem 250 \
                -gatk /Software/GenomeAnalysisTK.jar



mkdir -p pileupFiles
rm -f errors_${date}
rm -f listOfPileups
# Move files back to original processing directory
for file in $(ls ../*/alignments/hs37d5/bam_files/*sorted.grouped.duprm.cram); do

        # Obtain directory and sample name
        pattern="bam_files"              # recognised regex
        path=${file%%$pattern*}          # path until hs37d5
        prefix=`echo ${#path} + 10 | bc` # numerical start of bam_files in path
        sample=${file:$prefix:-26}       # sample name with trimming name
        sampleInd=${sample%%.*}          # sample name without trimming name


	# Better name pileup files (filters were set to zero because already processed)

        # Moves files to working directory otherwise filename in script will not correspond
        if [ -f ${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam ]; then
                mv ${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam \
                        ${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam
                mv ${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam.bai \
                        ${path}final_files/${sample}.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam.bai

		# Check for potential problems (files missing)
		if [ ! -f ${sampleInd}.HOA.hs37d5.bq30.mq0.0bp.homozyg.fam ]; then
			echo "${sampleInd} with potential issues" >> errors_${date}                        # keep track of potential errors
		else
			
			# Filter samples who showed encoding errors at indel realignment
                	# filter out samples that show unexpected European ancestry on PCA and excluded from original paper probably due to contamination
                	if [ "${sample}" != "BC25.mixedTrimmed" ]  \
                	&& [ "${sample}" != "BC29.mixedTrimmed" ]  \
                	&& [ "${sample}" != "AM71.mixedTrimmed" ]  \
                	&& [ "${sample}" != "MA577.mixedTrimmed" ] \
                	&& [ "${sample}" != "895.mixedTrimmed" ]   \
                	&& [ "${sample}" != "PAPV_106_FE.alreadyTrimmed" ]  \
                	&& [ "${sample}" != "PAPV_109_FE.alreadyTrimmed" ]  \
                	&& [ "${sample}" != "PAPV_137_FE.alreadyTrimmed" ]  \
                	&& [ "${sample}" != "PAPV_146_FE.alreadyTrimmed" ]  \
                	&& [ "${sample}" != "PAPV_156_HU.alreadyTrimmed" ]  \
                	&& [ "${sample}" != "PAPV_26_M3SD.alreadyTrimmed" ] \
                	&& [ "${sample}" != "PAPV_27_LTP.alreadyTrimmed" ]  \
                	&& [ "${sample}" != "PAPV_52_FE.alreadyTrimmed" ]   \
                	&& [ "${sample}" != "PAPV_53_M3ID.alreadyTrimmed" ] \
                	&& [ "${sample}" != "PAPV_57_M1II.alreadyTrimmed" ] \
                	&& [ "${sample}" != "PAPV_61_FE.alreadyTrimmed" ]; then
				echo "pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg" >> listOfPileups # prepare a list to merge all pileup files
			fi

			# Moves files to working directory otherwise filename in script will not correspond
			mv ${sampleInd}.HOA.hs37d5.bq30.mq0.0bp.homozyg.bed pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg.bed
			mv ${sampleInd}.HOA.hs37d5.bq30.mq0.0bp.homozyg.bim pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg.bim
			mv ${sampleInd}.HOA.hs37d5.bq30.mq0.0bp.homozyg.fam pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg.fam
			mv ${sampleInd}.HOA.hs37d5.bq30.mq0.0bp.homozyg.log pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg.log
			mv ${sampleInd}.HOA.hs37d5.bq30.mq0.0bp.homozyg.map pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg.map
			mv ${sampleInd}.HOA.hs37d5.bq30.mq0.0bp.homozyg.nosex pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg.nosex
			mv ${sampleInd}.HOA.hs37d5.bq30.mq0.0bp.homozyg.ped pileupFiles/${sampleInd}.HOA.hs37d5.bq30.mq25.34bp.homozyg.ped
        	fi
	fi
done



# Merge all published shotgun files into one
if [ -f listOfPileups ]; then
        first=`head -1 listOfPileups`                                                                 # get the name of the first file
        sed 1,1d listOfPileups > listOfPileup                                                         # filter out first file (otherwie will be merged twice)
        plink1.9 --bfile ${first} --merge-list listOfPileup --make-bed --out publishedShotgun_${date} # finish with merging
        rm listOfPileup                                                                               # clean up directory
fi


# NB: run labelGroups.sh to assign groups
# it will also merge with modern individuals from HOA
