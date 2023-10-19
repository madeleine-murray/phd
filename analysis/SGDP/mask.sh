#!/bin/bash

###########################################################################
#              Mask non-American segments from Mayan in SGDP              #
# Segments were identified with MOSAIC using HOA because more individuals #
#                refer to HOA folder for more information                 #
###########################################################################

date=010123

# Exclude the two Mayan individuals from SGDP
grep "Mayan" ../SGDP.database.Copan.${date}.fam > Mayan.fam
plink --bfile ../SGDP.database.Copan.${date} \
	--remove Mayan.fam \
	--make-bed --out SGDP.database.Copan.${date}.noMaya
grep "S_Mayan-1" ../SGDP.database.Copan.${date}.fam > Mayan1.fam
plink --bfile ../SGDP.database.Copan.${date} \
	--keep Mayan1.fam \
	--recode --out SGDP.database.Copan.${date}.justMaya1
grep "S_Mayan-2" ../SGDP.database.Copan.${date}.fam > Mayan2.fam
plink --bfile ../SGDP.database.Copan.${date} \
	--keep Mayan2.fam \
	--recode --out SGDP.database.Copan.${date}.justMaya2
rm Mayan.fam Mayan1.fam Mayan2.fam

# Identify non-American SNPs in Maya
python3 identifyNonAmericanSites.py
# change the input/output files to match individual of interest
# warning: script is not efficient and takes about 1 hour to run

# Mask non-American SNPs in Maya
awk '{print $2, "C1"}' mayan1.zero > maya1.zero
awk '{print $2, "C2"}' mayan2.zero > maya2.zero
awk '{print "Mayan", "S_Mayan-1", "C1"}' mayan1.zero | uniq > maya1.clst
awk '{print "Mayan", "S_Mayan-2", "C2"}' mayan2.zero | uniq > maya2.clst
plink --file SGDP.database.Copan.${date}.justMaya1 \
	--zero-cluster maya1.zero \
	--within maya1.clst \
	--make-bed --out SGDP.database.Copan.${date}.maskedMaya1
plink --file SGDP.database.Copan.${date}.justMaya2 \
	--zero-cluster maya2.zero \
	--within maya2.clst \
	--make-bed --out SGDP.database.Copan.${date}.maskedMaya2
#zzz.bwh.harvard.edu/plink/dataman.shtml#zero
#do each individual separetly for caution

# Merge masked Maya to the database
echo "SGDP.database.Copan.${date}.maskedMaya1" >  maya.list
echo "SGDP.database.Copan.${date}.maskedMaya2" >> maya.list
plink --bfile SGDP.database.Copan.${date}.noMaya \
	--merge-list maya.list \
	--allow-no-sex \
	--make-bed --out SGDP.database.Copan.${date}.masked

# Remove temporary files
rm SGDP.database.Copan.${date}.justMaya* SGDP.database.Copan.${date}.noMaya*
rm *zero *clst
rm maya.list SGDP.database.Copan.${date}.maskedMaya*

# Convert to human readable to check missingness
plink --bfile SGDP.database.Copan.${date}.masked --missing
grep "Mayan" plink.imiss
rm plink*
