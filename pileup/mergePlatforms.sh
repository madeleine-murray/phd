#!/bin/bash

: '
THIS SCRIPT MERGES SAMPLES PROCESSED FROM NOVASEQ AND MISEQ
NOVASEQ WAS PROCESSED ON PG4 IN /st1/hdd/pg/human_data_processing/platforms/novaseq_macrogen/shigeki/HN00147677_HN00147680
MISEQ WAS PROCESSED ON PG1 IN /home/madeleine/maya/miseq/automatic


# Samples were processed as such:
python /Programs/seqPro_v5-1_extraAdapters.py \
	-R /Reference_Genomes/hs37d5/hs37d5.fa -n 10 -p 3 -o seqPro_maya_stats \
        -mq 0 -trimbp 25 \
        -mem 500 -trimq 25 -data WGS \
        -index Y -fish N \
        -end PE -md N \
        -dam N -lib Y -RGCN kanazawa \
        -RGPM NovaSeq -seqid none
'

directory="/st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/mergePlatform"
rm -f inputSamples

# Merge mi-seq and nova-seq
for file in $(ls *Mi_*trimmed25bp25q.sorted.grouped.duprm.bam); do
	sample=${file:0:5}
	
	# Merge and index
	java -jar /Software/picard.jar MergeSamFiles \
		I=${file} \
		I=${directory}/${sample}-A-EX1-Nova_S00.trimmed25bp25q.sorted.grouped.duprm.bam \
		O=${sample}.trimmed25bp25q.bam \
		2> ${sample}.trimmed25bp25q.log
	samtools sort ${sample}.trimmed25bp25q.bam \
		-o ${sample}.trimmed25bp25q.sorted.bam
	samtools index ${sample}.trimmed25bp25q.sorted.bam
	samtools flagstat ${sample}.trimmed25bp25q.sorted.bam \
		> ${sample}.trimmed25bp25q.sorted.flagstat
	
	# Indel realignment
	java -jar /Software/GenomeAnalysisTK.jar \
        	-T RealignerTargetCreator \
        	-R /Reference_Genomes/hs37d5/hs37d5.fa \
        	-I ${sample}.trimmed25bp25q.sorted.bam \
        	-o ${sample}.trimmed25bp25q.sorted.indel.intervals
	java -jar /Software/GenomeAnalysisTK.jar \
        	-T IndelRealigner \
        	-R /Reference_Genomes/hs37d5/hs37d5.fa \
        	-I ${sample}.trimmed25bp25q.sorted.bam \
        	-o ${sample}.trimmed25bp25q.sorted.indel.bam \
        	-targetIntervals ${sample}.trimmed25bp25q.sorted.indel.intervals
	
	# Prepare sample list
	echo "${sample}.trimmed25bp25q.sorted.indel.bam" >> inputSamples
done

# Filter and softclip
python /Programs/filter_bam.py -mq 25 \
                                -minbp 34 \
                                -clip 2 \
                                -p 10 \
                                -i inputSamples

# Get coverage with qualimap
for file in $(ls *trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam); do
	/Software/qualimap_v2.2.1/qualimap \
		bamqc --java-mem-size=50G -bam ${file} &
done
