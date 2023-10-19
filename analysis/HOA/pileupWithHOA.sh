#!/bin/bash

date=07122022 # keep track of pileup version


############################################################
##### Pile up Copan with HOA and published Americans #######
############################################################

# Human Origin Array (HOA) was piled up with published shotgun from America on pg3:
if [ ! -f HOA_${date}.bed ] || [ ! -f HOA_${date}.pileup.bim ] || [ ! -f HOA_${date}.pileup.fam ]; then
        scp madeleine@134.226.152.253://external_storage/madeleine/Downloaded_Data/shotgun/Pileup_HOA/HOA-publishedShotgun_${date}.bim HOA_${date}.pileup.bim
        scp madeleine@134.226.152.253://external_storage/madeleine/Downloaded_Data/shotgun/Pileup_HOA/HOA-publishedShotgun_${date}.fam HOA_${date}.pileup.fam
	scp madeleine@134.226.152.253://external_storage/madeleine/Downloaded_Data/shotgun/Pileup_HOA/HOA-publishedShotgun_${date}.bed HOA_${date}.bed
	cat HOA_${date}.pileup.bim | awk '{print $1,$4-1,$4}' > HOA_${date}.pileup.bed # calculate starting SNP coordinates from end coordinates

	# Rename AG2 & MA1 from HOA to differ from AG2 & MA1 in published (same individual)
	sed -i 's/Siberian_Ice_Age AG2 0 0 1 1/Siberian_Ice_Age Ag2 0 0 1 1/g'                     HOA_${date}.pileup.fam
	sed -i 's/Siberian_Upper_Paleolithic MA1 0 0 1 1/Siberian_Upper_Paleolithic Ma1 0 0 1 1/g' HOA_${date}.pileup.fam
fi


# Prepare input Copan individuals
if [ ! -f ../mergePlatform/CpM0203.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam.bai ] || [ ! -f ../mergePlatform/CpM0506.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam.bai ] || [ ! -f ../mergePlatform/CpM11.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam.bai ] || [ ! -f ../mergePlatform/CpM12.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam.bai ] || [ ! -f ../mergePlatform/CpM13.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam.bai ] || [ ! -f ../mergePlatform/CpM34.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam.bai ] || [ ! -f ../mergePlatform/CpM4243.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam.bai ]; then
	samtools index ../mergePlatform/CpM0203.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
	samtools index ../mergePlatform/CpM0506.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
	samtools index ../mergePlatform/CpM11.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
	samtools index ../mergePlatform/CpM12.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
	samtools index ../mergePlatform/CpM13.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
	samtools index ../mergePlatform/CpM34.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
	samtools index ../mergePlatform/CpM4243.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam
fi
echo "../mergePlatform/CpM0203.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam"  > copanBamInput
echo "../mergePlatform/CpM0506.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam" >> copanBamInput
echo "../mergePlatform/CpM11.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam"   >> copanBamInput
echo "../mergePlatform/CpM12.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam"   >> copanBamInput
echo "../mergePlatform/CpM13.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam"   >> copanBamInput
echo "../mergePlatform/CpM34.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam"   >> copanBamInput
echo "../mergePlatform/CpM4243.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam" >> copanBamInput


# Pile up Copan individuals from HOA and published shotgun
# (mapping quality and read length is set to zero because already filtered in the previous directory)
python /Programs/HapSNP_v9_pg4.py \
                -data /st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/HOA/HOA_${date} \
                -I copanBamInput \
                -R /Reference_Genomes/hs37d5/hs37d5.fa \
                -p 7 \
                -o HOA-Copan_${date} \
                -mq 0 \
                -minbp 0 \
                -mem 250 \
                -gatk /Software/GenomeAnalysisTK.jar
		# HapSNP_v9.py to replace "plink1.9" to "plink"


# Assign information to Copan samples
echo "Copan CpM0203 0 0 2 1"  > HOA-Copan_${date}.HOA_${date}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM0506 0 0 2 1" >> HOA-Copan_${date}.HOA_${date}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM11 0 0 1 1"   >> HOA-Copan_${date}.HOA_${date}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM12 0 0 1 1"   >> HOA-Copan_${date}.HOA_${date}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM13 0 0 1 1"   >> HOA-Copan_${date}.HOA_${date}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM34 0 0 2 1"   >> HOA-Copan_${date}.HOA_${date}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM4243 0 0 2 1" >> HOA-Copan_${date}.HOA_${date}.hs37d5.bq30.mq0.0bp.homozyg.fam


# Merge Copan with HOA and published shotgun (rename bim and fam from HOA because different bed file than before)
echo "HOA_${date}" > file2
mv HOA_${date}.pileup.fam HOA_${date}.fam
mv HOA_${date}.pileup.bim HOA_${date}.bim
plink --bfile HOA-Copan_${date}.HOA_${date}.hs37d5.bq30.mq0.0bp.homozyg --merge-list file2 --make-bed --out HOA-PublishedShotgun-copan_${date}
rm file2


# Remove samples with no genotypes
# identified using the following command: plink --bfile HOA_${date} --recode --tab --missing --out HOA_${date}
echo "ASO CK01 0 0 0 -9"  > missingGenotypes.pop
echo "ASO CK02 0 0 0 -9" >> missingGenotypes.pop
plink --bfile HOA-PublishedShotgun-copan_${date} \
	--remove missingGenotypes.pop \
	--make-bed --out HOA-publishedShotgun-copan_${date}


# Remove files we won't be needing
rm HOA-PublishedShotgun-copan_${date}*
mv HOA_${date}.fam HOA_${date}.pileup.fam
mv HOA_${date}.bim HOA_${date}.pileup.bim
rm copanBamInput
rm -r CpM*
rm HOA-Copan_${date}.HOA_${date}.hs37d5.bq30.mq0.0bp*


# Remove individuals with too low coverage but keep Copan
#plink --bfile HOA-publishedShotgun-copan_${date} \
#	--mind 0.98 --make-bed \
#	--out HOA-publishedShotgun-copan_${date}_miss
awk '$6 > 0.98 && $1 != "Copan" {print $1,$2,"0","0","0","-9"}' HOA-publishedShotgun-copan_07122022_miss.imiss > lowCoverage.pop
plink --bfile HOA-publishedShotgun-copan_${date} \
        --remove lowCoverage.pop \
        --make-bed --out HOA-publishedShotgun-copan_${date}_miss

############################################################
############################################################
############################################################
