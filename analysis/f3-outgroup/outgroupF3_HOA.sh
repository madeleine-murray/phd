#!/bin/bash


datasetDate="07122022_miss"
#cat ../HOA/HOA-publishedShotgun-copan_07122022.fam | grep " 0 0 0 " | awk '{print $1}' | uniq | tr '\n' ' ' # for ancient Americans

subDir="all_HOA" #and with additional Asians: Han, Dai, Ami and Japanese
americans="Han Dai Ami_Coriell Japanese Aleut Chukchi Chukchi_Reindeer Chukchi_Sir Eskimo_Chaplin Eskimo_Naukan Eskimo_Sireniki Even Itelmen Karitiana Mayan Mixe Mixtec Piapoco Pima Quechua_Coriell Surui Tlingit Ulchi Zapotec Copan Aconcagua AG2 ancientPanama ASO Athabaskan Ayayema Baja_Mexico BigBar Caribbean_Taino Clovis Colonist Coropuna Fuego-Patagonian_Yaghan Huron-Wendat Island_Chumash Kaillachuro Kennewick Lagoa_Santa Late_Dorset Lovelock LSCI Lucy_Islands Lukurmata MA1 Mainland_Chumash Middle_Dorset Old_Mission_Point Oruro Pericu Pericues PuntaSantaAna Rio_Uncallane San_Francisco_Bay San_Nicolas_Island_Early San_Nicolas_Island_Late Soro_Mikaya_Patjxa SpiritCave Tiwanaku TrailCreek USR1 USR2"

# Make a list of all combinations for F3 with populations Y and X
mkdir -p ${subDir}
for y in ${americans}; do
	rm -f ${subDir}/${y}.X.Mbuti.F3.list
	for x in ${americans}; do
		echo ${y} ${x} MbutiPygmy >> ${subDir}/${y}.X.Mbuti.F3.list
	done
done

# Create the parameter file
for y in ${americans}; do
	echo -e genotypename: ' \t ' ../HOA/HOA-publishedShotgun-copan_${datasetDate}.eigenstratgeno > ${subDir}/${y}.X.Mbuti.F3.par
	echo -e snpname:      ' \t ' ../HOA/HOA-publishedShotgun-copan_${datasetDate}.snp           >> ${subDir}/${y}.X.Mbuti.F3.par
	echo -e indivname:    ' \t ' ../HOA/HOA-publishedShotgun-copan_${datasetDate}.ind           >> ${subDir}/${y}.X.Mbuti.F3.par
	echo -e popfilename:  ' \t ' ${subDir}/${y}.X.Mbuti.F3.list                                 >> ${subDir}/${y}.X.Mbuti.F3.par
done

# Run outgroup F3 statistics
for y in ${americans}; do
	/Software/AdmixTools-v6/bin/qp3Pop -p ${subDir}/${y}.X.Mbuti.F3.par > ${subDir}/${y}.X.Mbuti.F3.log
done

exit
# Dipslay results
grep "result\|no data" ${subDir}/*.X.Mbuti.F3.log | awk '{print $3,$4,$6}'                    | column -t
grep "result\|no data" ${subDir}/*.X.Mbuti.F3.log | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10}' | column -t
