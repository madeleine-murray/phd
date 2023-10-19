#!/bin/bash

# F3 between the Copan individual as if they were from a different family. This will allow us to do a heatmap.
# Here, we are using the 7 ancient Mayans from Copan and the two modern Mayans from SGDP


datasetDate=090922


# Make individuals from Copan-family as no-family
awk '{ if ($3=="Copan" || $3=="Mayan") print $1,$2,$1; else print $1,$2,$3; }' \
	../SGDP/SGDP.database.Copan.${datasetDate}.ind > SGDP.database.Copan.${datasetDate}_CopanAsInds.ind


# Make a list of individuals from the Copan family
copanIndividuals="CpM0203 CpM0506 CpM11 CpM12 CpM13 CpM34 CpM4243 S_Mayan-1 S_Mayan-2"


# Create the list of combinations for the comparison of each individuals
rm -f CopanY.CopanX.Mbuti.heatmap.list
for y in ${copanIndividuals}; do
	for x in ${copanIndividuals}; do
		echo ${y} ${x} Mbuti >> CopanY.CopanX.Mbuti.heatmap.list
	done
done


# Create the parameter file
echo -e genotypename: ' \t ' ../SGDP/SGDP.database.Copan.${datasetDate}.eigenstratgeno > CopanY.CopanX.Mbuti.heatmap.par
echo -e snpname:      ' \t ' ../SGDP/SGDP.database.Copan.${datasetDate}.snp           >> CopanY.CopanX.Mbuti.heatmap.par
echo -e indivname:    ' \t ' SGDP.database.Copan.${datasetDate}_CopanAsInds.ind       >> CopanY.CopanX.Mbuti.heatmap.par
echo -e popfilename:  ' \t ' CopanY.CopanX.Mbuti.heatmap.list                         >> CopanY.CopanX.Mbuti.heatmap.par


# Run outgroup F3 statistics
/Software/AdmixTools-v6/bin/qp3Pop -p CopanY.CopanX.Mbuti.heatmap.par > CopanY.CopanX.Mbuti.heatmap.log
