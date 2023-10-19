#!/bin/bash

# For files that were processed on the server Kanazawa because more threads
# but missing genome analysis toolkit to finish processing

parallel=19
sequencing="shotgun"

#directory="Nagele_2020_Genomic_insight_early_peopling_Caribbean"
directory="Raghaven_et_al_2014_Genetic_prehistory_of_New_World_arctic"
#directory="Raghaven_et_al_2015_Recent_history_native_americans"

mkdir -p ${sequencing}/$directory/alignments/hs37d5/final_files
for file in $(ls ${sequencing}/$directory/alignments/hs37d5/bam_files/*mixedTrimmed.sorted.grouped.duprm.bam); do

	# Obtain sample name
	pattern="bam_files"
	path=${file%%$pattern*}
	prefix=`echo ${#path} + 10 | bc`
	sample=${file:$prefix:-38}

	# Process post-merge
	samtools index ${sequencing}/${directory}/alignments/hs37d5/bam_files/${sample}.mixedTrimmed.sorted.grouped.duprm.bam
	java -jar /Software/GenomeAnalysisTK.jar \
	        -T RealignerTargetCreator \
        	-R /Reference_Genomes/hs37d5/hs37d5.fa \
        	-I ${sequencing}/${directory}/alignments/hs37d5/bam_files/${sample}.mixedTrimmed.sorted.grouped.duprm.bam \
        	-o ${sequencing}/${directory}/alignments/hs37d5/final_files/${sample}.mixedTrimmed.sorted.grouped.duprm.indel.intervals
	java -jar /Software/GenomeAnalysisTK.jar \
        	-T IndelRealigner \
        	-R /Reference_Genomes/hs37d5/hs37d5.fa \
        	-I ${sequencing}/${directory}/alignments/hs37d5/bam_files/${sample}.mixedTrimmed.sorted.grouped.duprm.bam \
        	-o ${sequencing}/${directory}/alignments/hs37d5/final_files/${sample}.mixedTrimmed.sorted.grouped.duprm.indel.bam \
        	-targetIntervals ${sequencing}/${directory}/alignments/hs37d5/final_files/${sample}.mixedTrimmed.sorted.grouped.duprm.indel.intervals

done
exit

rm -f listOfFiles # start a fresh file
for file in $(ls ${sequencing}/${directory}/alignments/hs37d5/final_files/*mixedTrimmed.sorted.grouped.duprm.indel.bam); do

        # Obtain sample name
        pattern="final_files"
        path=${file%%$pattern*}
        prefix=`echo ${#path} + 12 | bc`
        sample=${file:$prefix:-44}

        # Create list of samples
        echo "${sequencing}/${directory}/alignments/hs37d5/final_files/${sample}.mixedTrimmed.sorted.grouped.duprm.indel.bam" >> listOfFiles
done

# Filter and softclip
python /Programs/filter_bam.py -mq 25 \
                                -minbp 34 \
                                -clip 2 \
                                -p ${parallel} \
                                -i listOfFiles
