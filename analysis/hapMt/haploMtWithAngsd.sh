#!/bin/bash


copan=("CpM0203 CpM0506 CpM11 CpM12 CpM13 CpM34 CpM4243")


# Copy Copan samples from pg4
#scp madeleine@134.226.153.138://st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/mergePlatform/CpM0203.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam /external_storage/madeleine/analysis/
#scp madeleine@134.226.153.138://st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/mergePlatform/CpM0506.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam /external_storage/madeleine/analysis/
#scp madeleine@134.226.153.138://st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/mergePlatform/CpM11.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam /external_storage/madeleine/analysis/
#scp madeleine@134.226.153.138://st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/mergePlatform/CpM12.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam /external_storage/madeleine/analysis/
#scp madeleine@134.226.153.138://st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/mergePlatform/CpM13.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam /external_storage/madeleine/analysis/
#scp madeleine@134.226.153.138://st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/mergePlatform/CpM34.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam /external_storage/madeleine/analysis/
#scp madeleine@134.226.153.138://st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/mergePlatform/CpM4243.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam /external_storage/madeleine/analysis/


# And index
#for sample in ${copan}; do
#	samtools index ../${sample}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
#done


# Extract only mtDNA from BAM (and then index)
for sample in ${copan}; do
	samtools view -b ../${sample}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam MT \
                > ${sample}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.mt.bam
        samtools index ${sample}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.mt.bam
done


# Create a mtDNA consensus with angsd for each sample
for sample in ${copan}; do
        /Software/angsd/angsd -nThreads 7 \
		-i ${sample}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.mt.bam \
		-minMapQ 25 -minQ 20 -setMinDepth 5 \
		-doFasta 2 -doCounts 1 \
		-out ${sample}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.mt.bam
done


# Extract the fasta from the gz file and modify the name inside the fasta file
for sample in ${copan}; do
        gzip -d ${sample}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.mt.bam.fa.gz

        sed s/MT/${sample}/g ${sample}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.mt.bam.fa > ${sample}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.mt.bam.fasta
	rm ${sample}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.mt.bam.fa
done


# Concatenate all samples and load to haplogrep online
cat *trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.mt.bam.fasta > Copan_mt.fasta
exit


# Run haplogrep
#for sample in ${copan}; do
#	/Software/haplogrep2.2.9/haplogrep classify \
#		--in ${sample}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.mt.bam.fasta --format fasta \
#		--out ${sample}-mt_haplogroup
#done


# Concatenate data
#echo -e "SampleID\tHaplogroug\tRank\tQuality\tRange" > Copan_haplogroup
#for sample in ${copan}; do tail -n 1 ${sample}-mt_haplogroup >> Copan_haplogroup; done

#awk '{print $1, $2, $3, $4}' Copan_haplogroup
