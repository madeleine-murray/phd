#!/bin/bash


datasetDate=090922

#subDir="all_SGDP"
#americans="Aleut Chane Chukchi Eskimo_Chaplin Eskimo_Naukan Eskimo_Sireniki Even Itelman Karitiana Mayan Mixe Mixtec Piapoco Pima Quechua Surui Tlingit Ulchi Zapotec Aconcagua AG2 Alaskan_Athabaskan ancientPanama Anzick ASO ASO_4000 ASO_Huron-Wendat Ayayema Baja_Mexico BigBar Birnirk colonialPanama Colonist Copan COR EarlyArchaic_removed ESN Fuego-Patagonian Island_Chumash Kaillachuro Kennewick LateArchaic LateArchaic_removed LatePaleoindian LatePaleoindian_removed Late_Dorset Lovelock LSCI LSN LucyIslands LUK MA1 Mainland_Chumash Middle_Dorset MummyMaderasEncoC2 Norton OldMissionPoint ORU Paleoamericans Pericu Pre-Columbian Pre-Dorset PreLatePaleoindian_removed PuntaSantaAna Rio_Uncallane San_Francisco_Bay Saqqaq SerradaCapivara SMP SpiritCave Sumidouro Sumidouro_exclude Taino Thule TIW TrailCreek USR1 USR2"
#cat SGDP/SGDP.database.Copan.090922.fam | grep " 0 0 0 " | awk '{print $1}' | uniq | tr '\n' ' ' # for ancient Americans

subDir="above15000snp_SGDP" #and with additional Asians: Han, Dai, Ami and Japanese
americans="Han Dai Ami Japanese Aleut Chane Chukchi Eskimo_Chaplin Eskimo_Naukan Eskimo_Sireniki Even Itelman Karitiana Mayan Mixe Mixtec Piapoco Pima Quechua Surui Tlingit Ulchi Zapotec Aconcagua AG2 Alaskan_Athabaskan ancientPanama Anzick ASO ASO_4000 ASO_Huron-Wendat Ayayema Baja_Mexico BigBar Colonist Copan COR ESN Fuego-Patagonian Island_Chumash Kaillachuro Late_Dorset Lovelock LSCI LSN LucyIslands LUK MA1 Mainland_Chumash Middle_Dorset MummyMaderasEncoC2 ORU Paleoamericans Pericu PuntaSantaAna Rio_Uncallane San_Francisco_Bay SpiritCave Sumidouro Sumidouro_exclude Taino TIW TrailCreek USR1"


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
	echo -e indivname:    ' \t ' ../SGDP/SGDP.database.Copan.${datasetDate}.ind           >> ${subDir}/${y}.X.Mbuti.F3.par
	echo -e popfilename:  ' \t ' ${subDir}/${y}.X.Mbuti.F3.list                           >> ${subDir}/${y}.X.Mbuti.F3.par
done

# Run outgroup F3 statistics
for y in ${americans}; do
	/Software/AdmixTools-v6/bin/qp3Pop -p ${subDir}/${y}.X.Mbuti.F3.par > ${subDir}/${y}.X.Mbuti.F3.log
done


# Dipslay results
grep "result\|no data" ${subDir}/*.X.Mbuti.F3.log | awk '{print $3,$4,$6}' | column -t
grep "result\|no data" ${subDir}/*.X.Mbuti.F3.log | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10}' | column -t
