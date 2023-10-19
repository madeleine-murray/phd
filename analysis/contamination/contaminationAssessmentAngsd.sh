#!/bin/bash


# Control for contamination by comparing the proportion of X in Y chromosomes in males: CpM11, CpM12 and CpM13
# http://www.popgen.dk/angsd/index.php/Contamination


# Copy Copan males from pg4
scp madeleine@134.226.153.138://st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/mergePlatform/CpM11.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam /external_storage/madeleine/analysis/
scp madeleine@134.226.153.138://st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/mergePlatform/CpM12.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam /external_storage/madeleine/analysis/
scp madeleine@134.226.153.138://st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/mergePlatform/CpM13.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam /external_storage/madeleine/analysis/


# Index the bams
samtools index ../CpM11.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
samtools index ../CpM12.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
samtools index ../CpM13.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam


# Run angsd
/Software/angsd/angsd \
        -i ../CpM11.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam \
        -r X:5000000-154900000 \
        -doCounts 1 -iCounts 1 \
        -out CpM11_contamination
/Software/angsd/angsd \
        -i ../CpM12.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam \
        -r X:5000000-154900000 \
        -doCounts 1 -iCounts 1 \
        -out CpM12_contamination
/Software/angsd/angsd \
	-i ../CpM13.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam \
	-r X:5000000-154900000 \
	-doCounts 1 -iCounts 1 \
	-out CpM13_contamination


# Do jacknife in R
Rscript /Software/angsd/R/contamination.R \
        mapFile="/Software/angsd/RES/chrX.unique.gz" \
        hapFile="/Software/angsd/RES/HapMapChrX.gz" \
        countFile="/external_storage/madeleine/analysis/contamination/CpM13_contamination.icnts.gz" \
        mc.cores=24 > CpM11_contamination.results
Rscript /Software/angsd/R/contamination.R \
        mapFile="/Software/angsd/RES/chrX.unique.gz" \
        hapFile="/Software/angsd/RES/HapMapChrX.gz" \
        countFile="/external_storage/madeleine/analysis/contamination/CpM13_contamination.icnts.gz" \
        mc.cores=24 > CpM12_contamination.results
Rscript /Software/angsd/R/contamination.R \
	mapFile="/Software/angsd/RES/chrX.unique.gz" \
	hapFile="/Software/angsd/RES/HapMapChrX.gz" \
	countFile="/external_storage/madeleine/analysis/contamination/CpM13_contamination.icnts.gz" \
	mc.cores=24 > CpM13_contamination.results
