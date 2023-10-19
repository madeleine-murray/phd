#!/bin/bash

date=20092022
thread=10


# Convert plink to ped
if [ ! -f ../../HOA/HOA-publishedShotgun-copan_${date}.map ]; then
	plink --bfile ../../HOA/HOA-publishedShotgun-copan_${date} --recode --tab --mind 0 --out ../../HOA/HOA-publishedShotgun-copan_${date}
fi


# Create a list of populations to keep
grep "Sardinian\|Han \|Chukchi\|Eskimo\|Itelmen\|Koryak\|Nganasan\|Karitiana\|Mayan\|Mixe\|Mixtec\|Piapoco\|Pima\|Quechua\|Surui\|Zapotec\|Copan\|Tiwanaku\|Lukurmata\|Coropuna\|Oruro\|ancientPanama\|Rio_Uncallane" ../../HOA/HOA-publishedShotgun-copan_${date}.fam > modernAndAncient1_${date}.poplist
# Europeans: 27 Sardinian
# East Asians: 33 Han (so excluding Han_NChina)
# Siberians: 23 Chukchi, 22 Eskimo, 6 Itelmen, 9 Koryak, 11 Nganasan
# Americans: 12 Karitiana, 18 Mayan, 10 Mixe, 10 Mixtec, 4 Piapoco, 14 Pima, 5 Quechua, 8 Surui, 10 Zapotec
# Ancients: 7 Copan, and the only ancient ancestry we can capture with admixture (8 Tiwanaku, 4 Lukurmata, 4 Coropuna, 1 Oruro, 9 ancientPanama, 5 Rio_Uncallane)


# and filter for the populations to keep
plink --file ../../HOA/HOA-publishedShotgun-copan_${date} \
	--keep modernAndAncient1_${date}.poplist \
	--make-bed --out HOA_copan_modernAndAncient1_${date}


# Pruning: removing redundant loci in linkage disequilibrium
plink --bfile HOA_copan_modernAndAncient1_${date} \
	--noweb --indep-pairwise 50 5 0.2 \
	--out HOA_copan_modernAndAncient1_${date}_PRUNE

plink --bfile HOA_copan_modernAndAncient1_${date} \
	--noweb --extract HOA_copan_modernAndAncient1_${date}_PRUNE.prune.in \
	--make-bed --out HOA_copan_modernAndAncient1_${date}_PRUNE


# Convert to eigenstrat
plink --bfile HOA_copan_modernAndAncient1_${date}_PRUNE \
	--recode --out HOA_copan_modernAndAncient1_${date}_PRUNE # create ped
awk '{print $2,$3,$4,$5,$6}' HOA_copan_modernAndAncient1_${date}_PRUNE.fam \
	| awk '{printf("%1d %s\n", NR, $0)}' > HOA_copan_modernAndAncient1_${date}_PRUNE.pedind
scp HOA_copan_modernAndAncient1_${date}_PRUNE.bim HOA_copan_modernAndAncient1_${date}_PRUNE.pedsnp

echo -e genotypename:    ' \t ' HOA_copan_modernAndAncient1_${date}_PRUNE.ped             > converter.par
echo -e snpname:         ' \t ' HOA_copan_modernAndAncient1_${date}_PRUNE.pedsnp         >> converter.par
echo -e indivname:       ' \t ' HOA_copan_modernAndAncient1_${date}_PRUNE.pedind         >> converter.par
echo -e outputformat:    ' \t ' EIGENSTRAT                                               >> converter.par
echo -e genotypeoutname: ' \t ' HOA_copan_modernAndAncient1_${date}_PRUNE.eigenstratgeno >> converter.par
echo -e snpoutname:      ' \t ' HOA_copan_modernAndAncient1_${date}_PRUNE.snp            >> converter.par
echo -e indivoutname:    ' \t ' HOA_copan_modernAndAncient1_${date}_PRUNE.ind            >> converter.par
echo -e familynames:     ' \t ' NO                                                       >> converter.par
convertf -p converter.par


# Set up varying K with 10 iterations each
for k in {2..13}; do
	mkdir -p K${k}
	cd K${k}
	for i in {1..30}; do 
		admixture -j${thread} --cv -s time ../HOA_copan_modernAndAncient1_${date}_PRUNE.bed ${k} \
			| tee HOA_copan_modernAndAncient1_${date}_PRUNE_K${k}_$i.log
		mv HOA_copan_modernAndAncient1_${date}_PRUNE.${k}.P HOA_copan_modernAndAncient1_${date}_PRUNE.K${k}_$i.P
		mv HOA_copan_modernAndAncient1_${date}_PRUNE.${k}.Q HOA_copan_modernAndAncient1_${date}_PRUNE.K${k}_$i.Q
	done
	cd ..
done


# Extract and summarise error information accross the iterations
for k in {2..13}; do
        rm -f HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors_K${k}.txt
        for i in {1..30}; do
                errorLine=`grep "CV error" K${k}/HOA_copan_modernAndAncient1_${date}_PRUNE_K${k}_$i.log`    # find error for each iteration
                error=${errorLine:16}
                echo -e ${k} '\t' ${error} >> HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors_K${k}.txt # output error for each K
        done
done
echo -e 'K2 \t K3 \t K4 \t K5 \t K6 \t K7 \t K8 \t K9 \t K10 \t K11 \t K12 \t K13' > HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors.txt # create header of output
paste HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors_K2.txt    \
	HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors_K3.txt  \
	HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors_K4.txt  \
	HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors_K5.txt  \
	HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors_K6.txt  \
	HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors_K7.txt  \
	HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors_K8.txt  \
	HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors_K9.txt  \
	HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors_K10.txt \
	HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors_K11.txt \
	HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors_K12.txt \
	HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors_K13.txt \
	| awk '{print $2, $4, $6, $8, $10, $12, $14, $16, $18, $20, $22, $24}' \
	>> HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors.txt
rm HOA_copan_modernAndAncient1_${date}_PRUNE_CV-errors_K*                                                   # remove temporary output files


# Extract and summarise error information accross the iterations
for k in {2..13}; do
        rm -f HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors_K${k}.txt
        for i in {1..30}; do
                errorLine=`grep "^Loglikelihood" K${k}/HOA_copan_modernAndAncient1_${date}_PRUNE_K${k}_$i.log`      # find error for each iteration
                error=${errorLine:15}
                echo -e ${k} '\t' ${error} >> HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors_K${k}.txt # output error for each K
        done
done
echo -e 'K2 \t K3 \t K4 \t K5 \t K6 \t K7 \t K8 \t K9 \t K10 \t K11 \t K12 \t K13' > HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors.txt # create header of output
paste HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors_K2.txt    \
        HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors_K3.txt  \
        HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors_K4.txt  \
        HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors_K5.txt  \
        HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors_K6.txt  \
        HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors_K7.txt  \
        HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors_K8.txt  \
        HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors_K9.txt  \
        HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors_K10.txt \
        HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors_K11.txt \
        HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors_K12.txt \
        HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors_K13.txt \
        | awk '{print $2, $4, $6, $8, $10, $12, $14, $16, $18, $20, $22, $24}' \
        >> HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors.txt
rm HOA_copan_modernAndAncient1_${date}_PRUNE_likelihood-errors_K*                                                   # remove temporary output files


exit
perl ../AncestryPainter.pl \
	-q K9/HOA_copan_modernAndAncient1_${date}_PRUNE.K9_2.Q \
	-i HOA_copan_modernAndAncient1_${date}_PRUNE.fam \
	-t Copan \
	-c painter/color.file \
	-p painter/pop.order \
	-o painter/admixture_modernAndAncient1_K9_i2
