#!/bin/bash




#######
# HOA #
#######

datasetDate="07122022_miss"

# Change some formatting to avoid errors
awk '{print $1,$2,0,$4,$5,$6}' ../HOA/HOA-publishedShotgun-copan_${datasetDate}.bim \
	> HOA-publishedShotgun-copan_${datasetDate}-zero.bim # set all genetic distance to zero to avoid snp list order error
awk '{print $1,$2,$3,$4,$5,1}' ../HOA/HOA-publishedShotgun-copan_${datasetDate}.fam \
	> HOA-PublishedShotgun-copan_${datasetDate}.fam # otherwise ancient will be ignored


echo -e "genotypename: ../HOA/HOA-publishedShotgun-copan_${datasetDate}.bed"                > parfile_convertf
echo -e "snpname: HOA-publishedShotgun-copan_${datasetDate}-zero.bim"                      >> parfile_convertf
echo -e "indivname: HOA-PublishedShotgun-copan_${datasetDate}.fam"                         >> parfile_convertf
echo -e "outputformat: EIGENSTRAT"                                                         >> parfile_convertf
echo -e "genotypeoutname: ../HOA/HOA-publishedShotgun-copan_${datasetDate}.eigenstratgeno" >> parfile_convertf
echo -e "snpoutname: ../HOA/HOA-publishedShotgun-copan_${datasetDate}.snp"                 >> parfile_convertf
echo -e "indivoutname: ../HOA/HOA-PublishedShotgun-copan_${datasetDate}.ind"               >> parfile_convertf
convertf -p parfile_convertf

# Reformat ind file
sed 's/:/\t/' ../HOA/HOA-PublishedShotgun-copan_${datasetDate}.ind | awk '{print $2,$3,$1}' > ../HOA/HOA-publishedShotgun-copan_${datasetDate}.ind
rm ../HOA/HOA-PublishedShotgun-copan_${datasetDate}.ind

printf "Karitiana\nMayan\nMixe\nMixtec\nPiapoco\nPima\nQuechua_Coriell\nSurui\nZapotec\n" > modernAmericans_HOA.pop
printf "Mayan\nMixe\nMixtec\nPiapoco\nQuechua_Coriell\nZapotec\n" > clusterAmericans_HOA.pop
printf "Karitiana\nMayan\nMixe\nMixtec\nPiapoco\nPima\nQuechua_Coriell\nSurui\nZapotec\nItelmen\nEskimo_Naukan\nEskimo_Chaplin\nEskimo_Sireniki\nChukchi\nChukchi_Reindeer\nChukchi_Sir\nDai\nAmi_Coriell\nHan\nJapanese\nUlchi\nTlingit\nAleut\nEven" > modernAmericans_HOA.pop

# Run PCA
echo -e "genotypename: ../HOA/HOA-publishedShotgun-copan_${datasetDate}.eigenstratgeno" > parfile_HOA # input file
echo -e "snpname: ../HOA/HOA-publishedShotgun-copan_${datasetDate}.snp"                >> parfile_HOA # input file
echo -e "indivname: ../HOA/HOA-publishedShotgun-copan_${datasetDate}.ind"              >> parfile_HOA # input file
echo -e "evecoutname: pca_HOA.evec"                                                    >> parfile_HOA # output file
echo -e "evaloutname: pca_HOA.eval"                                                    >> parfile_HOA # output file
echo -e "poplistname: modernAmericans_HOA.pop"                                         >> parfile_HOA # infer eigenvectors using only those populations
echo -e "numoutevec: 10"                                                               >> parfile_HOA # number of eigenvectors to output (default: 10)
echo -e "numoutlieriter: 0"                                                            >> parfile_HOA # max number of outliers to remove (default: 5) --> don't remove ancient
echo -e "maxpops: 500"                                                                 >> parfile_HOA # maximum number of population to be analysed (default: 100)
echo -e "lsqproject: YES"                                                              >> parfile_HOA # option to deal with missing data --> important for ancient data
echo -e "killr2: YES"                                                                  >> parfile_HOA # eliminate SNPs in LD with nearby SNPs (default: NO)
echo -e "autoshrink: YES"                                                              >> parfile_HOA # avoid streching the axes of cumputed compared to non-computed samples
echo -e "r2thresh: 0.2"                                                                >> parfile_HOA # SNPs pairs with r-squared>value will have 1 member removed (default -1.0)
smartpca2 -p parfile_HOA > parfile_HOA.log &

# Run only on the main cluster (don't use Karitiana, Surui and Pima on computation)
#echo -e "genotypename: ../HOA/HOA-publishedShotgun-copan_${datasetDate}.eigenstratgeno" > parfile_HOA_cluster # input file
#echo -e "snpname: ../HOA/HOA-publishedShotgun-copan_${datasetDate}.snp"                >> parfile_HOA_cluster # input file
#echo -e "indivname: ../HOA/HOA-publishedShotgun-copan_${datasetDate}.ind"              >> parfile_HOA_cluster # input file
#echo -e "evecoutname: pca_HOA_cluster.evec"                                            >> parfile_HOA_cluster # output file
#echo -e "evaloutname: pca_HOA_cluster.eval"                                            >> parfile_HOA_cluster # output file
#echo -e "poplistname: clusterAmericans_HOA.pop"                                        >> parfile_HOA_cluster # infer eigenvectors using only those populations
#echo -e "numoutevec: 10"                                                               >> parfile_HOA_cluster # number of eigenvectors to output (default: 10)
#echo -e "numoutlieriter: 0"                                                            >> parfile_HOA_cluster # max number of outliers to remove --> don't remove ancient
#echo -e "maxpops: 500"                                                                 >> parfile_HOA_cluster # maximum number of population to be analysed (default: 100)
#echo -e "lsqproject: YES"                                                              >> parfile_HOA_cluster # option to deal with missing data --> important for ancient data
#echo -e "killr2: YES"                                                                  >> parfile_HOA_cluster # eliminate SNPs in LD with nearby SNPs (default: NO)
#echo -e "autoshrink: YES"                                                              >> parfile_HOA_cluster # avoid streching the axes of cumputed compared to non-computed
#echo -e "r2thresh: 0.2"                                                                >> parfile_HOA_cluster # SNPs pairs with r-squared>value will have 1 member removed
#smartpca2 -p parfile_HOA_cluster > parfile_HOA_cluster.log
exit



########
# SGDP #
########

datasetDate=090922

printf "Chane\nKaritiana\nMayan\nMixe\nMixtec\nPiapoco\nPima\nQuechua\nSurui\nZapotec\n" > modernAmericans_SGDP.pop

# Run PCA
echo -e "genotypename: ../SGDP/SGDP.database.Copan.${datasetDate}.eigenstratgeno" > parfile_SGDP # input file
echo -e "snpname: ../SGDP/SGDP.database.Copan.${datasetDate}.snp"                >> parfile_SGDP # input file
echo -e "indivname: ../SGDP/SGDP.database.Copan.${datasetDate}.ind"              >> parfile_SGDP # input file
echo -e "evecoutname: pca_SGDP.evec"                                             >> parfile_SGDP # output file
echo -e "evaloutname: pca_SGDP.eval"                                             >> parfile_SGDP # output file
echo -e "poplistname: modernAmericans_SGDP.pop"                                  >> parfile_SGDP # infer eigenvectors using only those populations
echo -e "numoutevec: 10"                                                         >> parfile_SGDP # number of eigenvectors to output (default: 10)
echo -e "numoutlieriter: 0"                                                      >> parfile_SGDP # maximum number of outliers to remove (default: 5) --> don't remove ancient
echo -e "maxpops: 500"                                                           >> parfile_SGDP # maximum number of population to be analysed (default: 100)
echo -e "lsqproject: YES"                                                        >> parfile_SGDP # option to deal with missing data --> important for ancient data
echo -e "killr2: YES"                                                            >> parfile_SGDP # eliminate SNPs in LD with nearby SNPs (default: NO)
echo -e "autoshrink: YES"                                                        >> parfile_SGDP # avoid streching the axes of cumputed samples compared to non-computed samples
echo -e "r2thresh: 0.2"                                                          >> parfile_SGDP # pairs of SNPs with r-squared>value will have 1 member removed (default: -1.0)
smartpca2 -p parfile_SGDP > parfile_SGDP.log

# Remove an outlier identified in PCA
sed 's/S_Mixtec-2 F Mixtec/S_Mixtec-2 F Mixtec_outlier/' ../SGDP/SGDP.database.Copan.${datasetDate}.ind > SGDP.database.Copan.${datasetDate}_outlier.ind
echo -e "genotypename: ../SGDP/SGDP.database.Copan.${datasetDate}.eigenstratgeno" > parfile_SGDP_outlier # input file
echo -e "snpname: ../SGDP/SGDP.database.Copan.${datasetDate}.snp"                >> parfile_SGDP_outlier # input file
echo -e "indivname: SGDP.database.Copan.${datasetDate}_outlier.ind"              >> parfile_SGDP_outlier # input file
echo -e "evecoutname: pca_SGDP_outlier.evec"                                     >> parfile_SGDP_outlier # output file
echo -e "evaloutname: pca_SGDP_outlier.eval"                                     >> parfile_SGDP_outlier # output file
echo -e "poplistname: modernAmericans_SGDP.pop"                                  >> parfile_SGDP_outlier # infer eigenvectors using only those populations
echo -e "numoutevec: 10"                                                         >> parfile_SGDP_outlier # number of eigenvectors to output (default: 10)
echo -e "numoutlieriter: 0"                                                      >> parfile_SGDP_outlier # max number of outliers to remove (default: 5) --> don't remove ancient
echo -e "maxpops: 500"                                                           >> parfile_SGDP_outlier # maximum number of population to be analysed (default: 100)
echo -e "lsqproject: YES"                                                        >> parfile_SGDP_outlier # option to deal with missing data --> important for ancient data
echo -e "killr2: YES"                                                            >> parfile_SGDP_outlier # eliminate SNPs in LD with nearby SNPs (default: NO)
echo -e "autoshrink: YES"                                                        >> parfile_SGDP_outlier # avoid streching the axes of cumputed compared to non-computed samples
echo -e "r2thresh: 0.2"                                                          >> parfile_SGDP_outlier # pairs of SNPs with r-squared>value will have 1 member removed
smartpca2 -p parfile_SGDP_outlier > parfile_SGDP_outlier.log




#########
# Allen #
#########

datasetDate=15112022

echo -e "genotypename: ../Allen/Allen-1240-Copan-published_${datasetDate}.bed"                > parfile_convertf
echo -e "snpname: ../Allen/Allen-1240-Copan-published_${datasetDate}.bim"                    >> parfile_convertf
echo -e "indivname: ../Allen/Allen-1240-Copan-published_${datasetDate}.fam"                  >> parfile_convertf
echo -e "outputformat: EIGENSTRAT"                                                           >> parfile_convertf
echo -e "genotypeoutname: ../Allen/Allen-1240-Copan-published_${datasetDate}.eigenstratgeno" >> parfile_convertf
echo -e "snpoutname: ../Allen/Allen-1240-Copan-published_${datasetDate}.snp"                 >> parfile_convertf
echo -e "indivoutname: ../Allen/Allen-1240-copan-published_${datasetDate}.ind"               >> parfile_convertf
convertf -p parfile_convertf

# List modern americans (grep for all Americans keyword but filtered for moderns only with awk)
tail -n+2 ../Allen/v52.2_1240K_public.anno \
	| awk -F'\t' '{ if ($8=="Modern") print $13}' \
	| grep "ACB\|Argentina\|ASW\|Aymara\|Bahamas\|Mexico\|Belize\|Bolivia\|Brazil\|Canada\|Chane\|Chile\|Chipewyan\|CLM\|Cree\|Cuba\|Curacao\|Dominican\|GIH\|Greenland\|Guadeloupe\|Haiti\|Hawaiian\|Huichol\|Karitiana\|Piapoco\|Mayan\|Mixe\|Mixtec\|MXL\|Nahua\|Panama\|PEL\|Peru\|Piapoco\|Pima\|PuertoRico\|PUR\|Quechua\|StLucia\|Surui\|Tsimshian\|USA\|Venezuela\|Yukpa" \
	| sort | uniq \
	> modernAmericans_Allen.pop
sed -i 's/Ignore_Karitiana(discovery).DG/Ignore_Karitiana.DG/g' modernAmericans_Allen.pop
sed -i 's/Ignore_Karitiana(discovery).SDG/Ignore_Karitiana.SDG/g' modernAmericans_Allen.pop

# Reformat ind file
sed 's/:/\t/' ../Allen/Allen-1240-copan-published_${datasetDate}.ind | awk '{print $2,$3,$1}' > ../Allen/Allen-1240-Copan-published_${datasetDate}.ind
rm ../Allen/Allen-1240-copan-published_${datasetDate}.ind

# Run PCA
echo -e "genotypename: ../Allen/Allen-1240-Copan-published_${datasetDate}.eigenstratgeno" > parfile_Allen # input file
echo -e "snpname: ../Allen/Allen-1240-Copan-published_${datasetDate}.snp"                >> parfile_Allen # input file
echo -e "indivname: ../Allen/Allen-1240-Copan-published_${datasetDate}.ind"              >> parfile_Allen # input file
echo -e "evecoutname: pca_Allen.evec"                                                    >> parfile_Allen # output file
echo -e "evaloutname: pca_Allen.eval"                                                    >> parfile_Allen # output file
echo -e "poplistname: modernAmericans_Allen.pop"                                         >> parfile_Allen # infer eigenvectors using only those populations
echo -e "numoutevec: 10"                                                                 >> parfile_Allen # number of eigenvectors to output (default: 10)
echo -e "numoutlieriter: 0"                                                              >> parfile_Allen # maximum number of outliers to remove (default: 5)
echo -e "maxpops: 500"                                                                   >> parfile_Allen # maximum number of population to be analysed (default: 100)
echo -e "lsqproject: YES"                                                                >> parfile_Allen # option to deal with missing data --> important for ancient data
echo -e "killr2: YES"                                                                    >> parfile_Allen # eliminate SNPs in LD with nearby SNPs (default: NO)
echo -e "autoshrink: YES"                                                                >> parfile_Allen # avoid streching the axes of cumputed compared to non-computed samples
echo -e "r2thresh: 0.2"                                                                  >> parfile_Allen # pairs of SNPs with r-squared>value will have 1 removed (default=-1.0)
smartpca2 -p parfile_Allen > parfile_Allen.log

