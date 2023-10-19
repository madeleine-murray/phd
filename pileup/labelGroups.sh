#!/bin/bash

date=07122022 # keep track of pileup version


########## Label published individuals ##########

# Correct manually some mislabeling
cat publishedShotgun_${date}.fam \
	| sed 's/B_003_1_B_003_1/B_003_1_SNPCAP_UDG_MIX/' \
	| sed 's/NiNg_1_A6_NiNg_1_A6/NiNg_1_A6_SNPCAP_UDG_MIX/' \
	| sed 's/XIV_C3_4_XIV_C3_4/XIV_C3_4_SNPCAP_UDG_MIX/' \
	| sed 's/Shuka_SNPCAP_UDG_MIX/Shuka-Kaa_SNPCAP_UDG_MIX/' \
	> publishedShotgun_${date}_corr.fam

# Crop sample name to only have sample name
awk '{print $1, substr($2,1,length($2)-15), $3, $4, $5, $6}' publishedShotgun_${date}_corr.fam > publishedShotgun_${date}_cor.fam
#rm publishedShotgun_${date}_corr.fam

# Link with sample names with group labels
rm -f publishedShotgun_${date}_lab.fam
while IFS="" read -r sample || [ -n "$sample" ]; do
	sampleName=`echo ${sample:3:-9}`
	awk '{if ($2=="'${sampleName}'") print($1,$2,"0","0","0","-9");}' groupsLabel.fam >> publishedShotgun_${date}_lab.fam
done < publishedShotgun_${date}_cor.fam

# Check if all pileup samples have been labeled and output a list of those who have not been labeled
rm -f publishedShotgun_${date}_lab.log
cat publishedShotgun_${date}_cor.fam groupsLabel.fam \
	| awk '{print $2}' | sort | uniq -c \
	| awk '{if ($1=="1") print($2, "has no label or was no processed !");}' >> publishedShotgun_${date}_lab.log

# Finish
rm publishedShotgun_${date}_cor.fam

: '
# Normal that following individuals are absent: too low number of reads (I think)
CK04
CK08
PS11
PS19
PS24
SN02
'

##############################################
########## Merge published with HOA ##########

# Prepare file name for HOA
if [ ! -f HOA.database.bed ]; then
	scp madeleine@134.226.153.112://home/niallc/Datasets/Lazaridis/Lazaridis.rmsnpids.90sRm.letters.Mathmissnprm.nosex.names_fixed/Lazaridis.rmsnpids.90sRm.letters.Mathmissnprm.nosex.names_fixed.bed HOA.database.bed
fi
mv HOA.database.pileup.bim HOA.database.bim
mv HOA.database.pileup.fam HOA.database.fam
echo "HOA.database" > file2

# Prepare file name for published shotgun to include label without overwriting original
mv publishedShotgun_${date}.bed publishedShotgun_${date}_lab.bed
mv publishedShotgun_${date}.bim publishedShotgun_${date}_lab.bim

# Merge HOA with published shotgun
plink1.9 --bfile publishedShotgun_${date}_lab --merge-list file2 --make-bed --out HOA-publishedShotgun_${date}

# Reset original names
rm file2
mv publishedShotgun_${date}_lab.bed publishedShotgun_${date}.bed
mv publishedShotgun_${date}_lab.bim publishedShotgun_${date}.bim
mv HOA.database.bim HOA.database.pileup.bim
mv HOA.database.fam HOA.database.pileup.fam
