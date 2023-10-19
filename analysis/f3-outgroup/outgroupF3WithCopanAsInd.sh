#!/bin/bash


datasetDate=090922

subDir="mainClusterWithCopanAsInd"
americans="CpM0203 CpM0506 CpM11 CpM12 CpM13 CpM34 CpM4243 PreLatePaleoindian_removed EarlyArchaic_removed LatePaleoindian LateArchaic_removed LateArchaic LatePaleoindian_removed Kennewick Paleoamericans Baja_Mexico Kaillachuro Pre-Columbian San_Francisco_Bay Pericu Taino Quechua Karitiana Chane SpiritCave Sumidouro_exclude Mainland_Chumash Island_Chumash LSN LSCI ESN TIW PuntaSantaAna Ayayema Sumidouro ancientPanama Surui Piapoco Mayan Mixe Zapotec Pima Lovelock Rio_Uncallane Aconcagua COR LUK ORU Fuego-Patagonian ASO_4000 ASO_Huron-Wendat ASO Anzick SMP BigBar LucyIslands Alaskan_Athabaskan Mixtec OldMissionPoint USR2"

# Make individuals from Copan-family as no-family
awk '{ if ($3=="Copan") print $1,$2,$1; else print $1,$2,$3; }' \
	../SGDP/SGDP.database.Copan.${datasetDate}.ind > SGDP.database.Copan.${datasetDate}_CopanAsInds.ind


# Make a list of all combinations for F3 with populations Y and X
mkdir -p ${subDir}
for y in ${americans}; do
	rm -f ${subDir}/${y}.X.Mbuti.F3.list
	for x in ${americans}; do
		echo ${y} ${x} Mbuti >> ${subDir}/${y}.X.Mbuti.F3.list
	done
done


# Create the parameter file
for y in ${americans}; do
	echo -e genotypename: ' \t ' ../SGDP/SGDP.database.Copan.${datasetDate}.eigenstratgeno > ${subDir}/${y}.X.Mbuti.F3.par
	echo -e snpname:      ' \t ' ../SGDP/SGDP.database.Copan.${datasetDate}.snp           >> ${subDir}/${y}.X.Mbuti.F3.par
	echo -e indivname:    ' \t ' SGDP.database.Copan.${datasetDate}_CopanAsInds.ind       >> ${subDir}/${y}.X.Mbuti.F3.par
	echo -e popfilename:  ' \t ' ${subDir}/${y}.X.Mbuti.F3.list                           >> ${subDir}/${y}.X.Mbuti.F3.par
done


# Run outgroup F3 statistics
for y in ${americans}; do
	/Software/AdmixTools-v6/bin/qp3Pop -p ${subDir}/${y}.X.Mbuti.F3.par > ${subDir}/${y}.X.Mbuti.F3.log
done


# Dipslay results
grep "result\|no data" ${subDir}/*.X.Mbuti.F3.log | awk '{print $3,$4,$6}' | column -t
grep "result\|no data" ${subDir}/*.X.Mbuti.F3.log | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10}' | column -t
