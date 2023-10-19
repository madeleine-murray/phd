#!/bin/bash

###########################################################################
#              Mask non-American segments from Mayan in SGDP              #
# Segments were identified with MOSAIC using HOA because more individuals #
#                refer to HOA folder for more information                 #
###########################################################################

date=010123
# conda activate localAncestry # to load pandas in python script

# Exclude the two Mayan individuals from SGDP
grep "Mayan" ../SGDP.database.Copan.${date}.fam > Mayan.fam
plink --bfile ../SGDP.database.Copan.${date} \
	--remove Mayan.fam \
	--make-bed --out SGDP.database.Copan.${date}.noMaya
plink --bfile ../SGDP.database.Copan.${date} \
	--keep Mayan.fam \
	--recode --out SGDP.database.Copan.${date}.justMaya

# Mask non-American segements using a python script
python3 createMaskedPed.py

# Convert the masked ped to plink and merge for the rest
plink --file SGDP.database.Copan.${date}.justMaya \
	--make-bed --out SGDP.database.Copan.${date}.justMaya
echo "SGDP.database.Copan.${date}.justMaya" > maya.list
plink --bfile SGDP.database.Copan.${date}.noMaya \
	--merge-list maya.list \
	--allow-no-sex \
	--make-bed --out SGDP.database.Copan.${date}.masked
rm Mayan.fam maya.list
rm SGDP.database.Copan.${date}.justMaya*
rm SGDP.database.Copan.${date}.noMaya*

#awk '{for (i=1; i<=NF; i++) { if (i >= 111428 && i <= 122424) printf "%s ", $i } print ""}' SGDP.database.Copan.010123.maskedMaya.ped | tail -1
