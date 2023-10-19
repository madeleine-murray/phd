#!/bin/bash

: '
THIS SCRIPT MERGES SAMPLES INTO INDIVIDUALS
INDIVUDALS WERE INDIFIED FROM READ /st1/hdd/pg/human_data_processing/platforms/novaseq_macrogen/shigeki/READ/READ_results


# Merge samples into individuals (identical samples)
java -jar /Software/picard.jar MergeSamFiles \
	I=CpM02.trimmed25bp25q.bam \
	I=CpM03.trimmed25bp25q.bam \
	O=CpM0203.trimmed25bp25q.bam \
	2> CpM0203.trimmed25bp25q.log
java -jar /Software/picard.jar MergeSamFiles \
        I=CpM05.trimmed25bp25q.bam \
        I=CpM06.trimmed25bp25q.bam \
        O=CpM0506.trimmed25bp25q.bam \
        2> CpM0506.trimmed25bp25q.log
java -jar /Software/picard.jar MergeSamFiles \
        I=CpM42.trimmed25bp25q.bam \
        I=CpM43.trimmed25bp25q.bam \
        O=CpM4243.trimmed25bp25q.bam \
        2> CpM4243.trimmed25bp25q.log
'

individuals=("CpM0203 CpM0506 CpM4243")
for individual in ${individuals}; do
	
	# Sort and index
	samtools sort ${individual}.trimmed25bp25q.bam \
		-o ${individual}.trimmed25bp25q.sorted.bam
	samtools index ${individual}.trimmed25bp25q.sorted.bam
	samtools flagstat ${individual}.trimmed25bp25q.sorted.bam \
		> ${individual}.trimmed25bp25q.sorted.flagstat
	
	# Indel realignment
	java -jar /Software/GenomeAnalysisTK.jar \
        	-T RealignerTargetCreator \
        	-R /Reference_Genomes/hs37d5/hs37d5.fa \
        	-I ${individual}.trimmed25bp25q.sorted.bam \
        	-o ${individual}.trimmed25bp25q.sorted.indel.intervals
	java -jar /Software/GenomeAnalysisTK.jar \
        	-T IndelRealigner \
        	-R /Reference_Genomes/hs37d5/hs37d5.fa \
        	-I ${individual}.trimmed25bp25q.sorted.bam \
        	-o ${individual}.trimmed25bp25q.sorted.indel.bam \
        	-targetIntervals ${individual}.trimmed25bp25q.sorted.indel.intervals
	
	# Prepare sample list
	echo "${individual}.trimmed25bp25q.sorted.indel.bam" >> inputIndividuals
done

# Filter and softclip
python /Programs/filter_bam.py -mq 25 \
                                -minbp 34 \
                                -clip 2 \
                                -p 10 \
                                -i inputIndividuals

# Get coverage with qualimap
for individual in ${individuals}; do
	/Software/qualimap_v2.2.1/qualimap \
		bamqc --java-mem-size=50G -bam ${individual}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam &
done
