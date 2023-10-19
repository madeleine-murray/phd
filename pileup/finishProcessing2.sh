#!/bin/bash

# For files that were processed on the server Kanazawa because more threads
# but missing genome analysis toolkit to finish processing

parallel=4
sequencing="shotgun"

#directory="Posth_et_al_2018_Reconstructing_deep_history"
#directory="Nakatsuka_et_al_2020_Paleogenomic_reconstruction_deep_population_history_Andes"
#directory="Nakatsuka_et_al_2020_Ancient_genomes_in_South_Patagonia"
#directory="Fernandes_et_al_2020_Genetic_history_pre-contact_Caribbean"
#directory="Lindo_et_al_2022_Uruguay"
directory="Fuente_et_al_2018_Genomic_insights_of_late_maritime_hunter-gatherers_from_the_Chilean_Patagonia"

mkdir -p $sequencing/$directory/alignments/hs37d5/final_files
for file in $(ls $sequencing/$directory/alignments/hs37d5/bam_files/*alreadyTrimmed.sorted.grouped.duprm.bam); do

	# Obtain sample name
	pattern="bam_files"
	path=${file%%$pattern*}
	prefix=`echo ${#path} + 10 | bc`
	sample=${file:$prefix:-40}

	# Process post-merge
	samtools index ${sequencing}/${directory}/alignments/hs37d5/bam_files/${sample}.alreadyTrimmed.sorted.grouped.duprm.bam
	java -jar /Software/GenomeAnalysisTK.jar \
	        -T RealignerTargetCreator \
        	-R /Reference_Genomes/hs37d5/hs37d5.fa \
        	-I ${sequencing}/${directory}/alignments/hs37d5/bam_files/${sample}.alreadyTrimmed.sorted.grouped.duprm.bam \
        	-o ${sequencing}/${directory}/alignments/hs37d5/final_files/${sample}.alreadyTrimmed.sorted.grouped.duprm.indel.intervals
	java -jar /Software/GenomeAnalysisTK.jar \
        	-T IndelRealigner \
        	-R /Reference_Genomes/hs37d5/hs37d5.fa \
        	-I ${sequencing}/${directory}/alignments/hs37d5/bam_files/${sample}.alreadyTrimmed.sorted.grouped.duprm.bam \
        	-o ${sequencing}/${directory}/alignments/hs37d5/final_files/${sample}.alreadyTrimmed.sorted.grouped.duprm.indel.bam \
        	-targetIntervals ${sequencing}/${directory}/alignments/hs37d5/final_files/${sample}.alreadyTrimmed.sorted.grouped.duprm.indel.intervals

done


rm -f listOfFiles # start a fresh file
for file in $(ls ${sequencing}/${directory}/alignments/hs37d5/final_files/*alreadyTrimmed.sorted.grouped.duprm.indel.bam); do

        # Obtain sample name
        pattern="final_files"
        path=${file%%$pattern*}
        prefix=`echo ${#path} + 12 | bc`
        sample=${file:$prefix:-46}

        # Create list of samples
        echo "${sequencing}/${directory}/alignments/hs37d5/final_files/${sample}.alreadyTrimmed.sorted.grouped.duprm.indel.bam" >> listOfFiles
done


# Filter and softclip
python /Programs/filter_bam.py -mq 25 \
                                -minbp 34 \
                                -clip 2 \
                                -p ${parallel} \
                                -i listOfFiles
