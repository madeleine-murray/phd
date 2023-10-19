#!/bin/bash

# Create list of individuals to determine sex
individuals=("CpM0203 CpM0506 CpM11 CpM12 CpM13 CpM34 CpM4243")
rm -f listIndividuals
for ind in $individuals; do
	echo "../mergePlatform/${ind}.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam" >> listIndividuals
done


# Determine sex using script: filters are set to zero because bam files were already filtered and it is not logical to re-filter
python /Programs/sampleSex_v2.py \
		-mq 0 \
		-i listIndividuals \
		-o sexing \
		-bp 0
