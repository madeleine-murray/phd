#!/bin/bash

# Allen dataset downloaded from Reich's lab:
# https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data

datasetDate=15112022

if [ ! -f v52.2_1240K_public.anno ]; then
	wget https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V52/V52.2/SHARE/public.dir/v52.2_1240K_public.tar
	tar -xf v52.2_1240K_public.tar
	rm v52.2_1240K_public.tar
fi

# Create plink pileup files: format by excluding chromosomes 23 and 24 and re-organise columns
awk '{ if($2<23) printf("%d %d %d\n",                $2,$4-1,$4)       ;}' v52.2_1240K_public.snp > Allen-1240_${datasetDate}.pileup.bed
awk '{ if($2<23) printf("%d\t%d:%d\t0\t%d\t%s\t%s\n",$2,$2,$4,$4,$5,$6);}' v52.2_1240K_public.snp > Allen-1240_${datasetDate}.pileup.bim
tail -n+2 v52.2_1240K_public.anno | awk -F'\t' '{print $13,$3,"0","0","9","1"}'                   > Allen-1240_${datasetDate}.pileup.fam


################################################
# Pile up Copan individuals with Allen dataset #
################################################

# Input file with Copan
echo "../mergePlatform/CpM0203.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam"  > copanBamInput
echo "../mergePlatform/CpM0506.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam" >> copanBamInput
echo "../mergePlatform/CpM11.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam"   >> copanBamInput
echo "../mergePlatform/CpM12.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam"   >> copanBamInput
echo "../mergePlatform/CpM13.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam"   >> copanBamInput
echo "../mergePlatform/CpM34.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam"   >> copanBamInput
echo "../mergePlatform/CpM4243.trimmed25bp25q.sorted.indel.mq25.clip2bp.34bp.bam" >> copanBamInput
# (mapping quality and read length is set to zero because already filtered in the previous directory)
python /Programs/HapSNP_v9_pg4.py \
                -data /st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/Allen/Allen-1240_${datasetDate} \
                -I copanBamInput \
                -R /Reference_Genomes/hs37d5/hs37d5.fa \
                -p 7 \
                -o Allen-1240-Copan_${datasetDate} \
                -mq 0 \
                -minbp 0 \
                -mem 250 \
                -gatk /Software/GenomeAnalysisTK.jar
                # HapSNP_v9.py to replace "plink1.9" to "plink"

# Create binary bed file
echo -e "genotypename: v52.2_1240K_public.geno"                  > parfile_convertf
echo -e "snpname: v52.2_1240K_public.snp"                       >> parfile_convertf
echo -e "indivname: v52.2_1240K_public.ind"                     >> parfile_convertf
echo -e "outputformat: PACKEDPED"                               >> parfile_convertf
echo -e "genotypeoutname: Allen-1240_${datasetDate}.pileup.bed" >> parfile_convertf
echo -e "snpoutname: Allen-1240_${datasetDate}.bim"             >> parfile_convertf
echo -e "indivoutname: Allen-1240_${datasetDate}.fam"           >> parfile_convertf
convertf -p parfile_convertf
rm Allen-1240_${datasetDate}.bim Allen-1240_${datasetDate}.fam parfile_convertf

# Assign information to Copan samples
echo "Copan CpM0203 0 0 2 1"  > Allen-1240-Copan_${datasetDate}.Allen-1240_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM0506 0 0 2 1" >> Allen-1240-Copan_${datasetDate}.Allen-1240_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM11 0 0 1 1"   >> Allen-1240-Copan_${datasetDate}.Allen-1240_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM12 0 0 1 1"   >> Allen-1240-Copan_${datasetDate}.Allen-1240_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM13 0 0 1 1"   >> Allen-1240-Copan_${datasetDate}.Allen-1240_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM34 0 0 2 1"   >> Allen-1240-Copan_${datasetDate}.Allen-1240_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM4243 0 0 2 1" >> Allen-1240-Copan_${datasetDate}.Allen-1240_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam

# Filter for modern American (no African ancestry) from Allen and merge with Copan and other processed ancient
tail -n+2 v52.2_1240K_public.anno \
	#| awk -F'\t' '{ if ($15=="Argentina" || $15=="Bahamas" || $15=="Belize" || $15=="Bolivia" || $15=="Brazil" || $15=="Canada" || $15=="Chile" || $15=="Colombia" || $15=="Cuba" || $15=="Curacao" || $15=="Greenland" || $15=="Guadeloupe" || $15=="Haiti" || $15=="Mexico" || $15=="Panama" || $15=="Peru" || $15=="Puerto Rico" || $15=="St. Lucia" || $15=="USA" || $15=="Venezuela") print $0}' \
	| awk -F'\t' '{ if ($8=="Modern") print $13,$3}' \
	| grep "Argentina\|Aymara\|Bahamas\|Mexico\|Belize\|Bolivia\|Brazil\|Canada\|Chane\|Chile\|Chipewyan\|CLM\|Cree\|Cuba\|Curacao\|Dominican\|Greenland\|Guadeloupe\|Haiti\|Hawaiian\|Huichol\|Karitiana\|Piapoco\|Mayan\|Mixe\|Mixtec\|MXL\|Nahua\|Panama\|PEL\|Peru\|Piapoco\|Pima\|PuertoRico\|PUR\|Quechua\|StLucia\|Surui\|Tsimshian\|USA\|Venezuela\|Yukpa" \
	> modernAmericans_Allen.pop
awk '{print $1,$2}' Allen-1240-Copan_${datasetDate}.Allen-1240_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam >> modernAmericans_Allen.pop
#grep " -9" ../HOA/HOA-publishedShotgun-copan_07122022.fam | awk '{print $1}' | sort | uniq >> modernAmericans_Allen.pop # add processed ancients
echo "Allen-1240_${datasetDate}.pileup"            > file2
#echo "../HOA/HOA-publishedShotgun-copan_07122022" >> file2 # merge with HOA & processed
plink --bfile Allen-1240-Copan_${datasetDate}.Allen-1240_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg \
        --merge-list file2 \
	--make-bed --allow-no-sex \
        --keep modernAmericans_Allen.pop \
        --out Allen-1240-Copan-published_${datasetDate}
sed -i 's/Ignore_Karitiana(discovery).DG/Ignore_Karitiana.DG/g' Allen-1240-Copan-published_${datasetDate}.fam
sed -i 's/Ignore_Karitiana(discovery).SDG/Ignore_Karitiana.SDG/g' Allen-1240-Copan-published_${datasetDate}.fam
# Clean up directory
rm -r CpM* file2
rm Allen-1240-Copan_${datasetDate}.Allen-1240_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg*
exit


# Merge Copan with HOA and published shotgun (different bed file than before)
# and fiter for American populations because so many names to changes later on
grep "ACB\|Argentina\|ASW\|Aymara\|Bahamas\|Mexico\|Belize\|Bolivia\|Brazil\|Canada\|Chane\|Chile\|Chipewyan\|CLM\|Cree\|Cuba\|Curacao\|Dominican\|GIH\|Greenland\|Guadeloupe\|Haiti\|Hawaiian\|Huichol\|Karitiana\|Piapoco\|Mayan\|Mixe\|Mixtec\|MXL\|Nahua\|Panama\|PEL\|Peru\|Piapoco\|Pima\|PuertoRico\|PUR\|Quechua\|StLucia\|Surui\|Tsimshian\|USA\|Venezuela\|Yukpa" Allen-1240-Copan_${datasetDate}.fam | awk '{print $1,$2}' > americans.pop
awk '{print $1,$2}' Allen-1240-Copan_${datasetDate}.Allen-1240_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam >> americans.pop
echo "Allen-1240_${datasetDate}.pileup" > file2
plink --bfile Allen-1240-Copan_${datasetDate}.Allen-1240_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg \
	--merge-list file2 --make-bed --allow-no-sex \
	--keep americans.pop \
	--out Allen-1240-Copan_${datasetDate}

# Replace some family names otherwise too long and further analysis will crash
sed -i                                                 's/Argentina_Aconcagua_Inca_500BP.SG/Argentina_Inca_500BP.SG/g'              Allen-1240-Copan_${datasetDate}.fam
sed -i                                   's/Argentina_ArroyoSeco2_7200BP_1d.or.2d.rel.I7090/Argentina_ArroyoSeco2_7200BP.7090/g'    Allen-1240-Copan_${datasetDate}.fam #relative
sed -i                                               's/Argentina_ArroyoSeco2_7400BP_contam/Argentina_ArroyoSeco2_7400BP_ctam/g'    Allen-1240-Copan_${datasetDate}.fam # CONTAM
sed -i                                              's/Argentina_BeagleChannel_Yamana_100BP/Argentina_Yamana_100BP/g'               Allen-1240-Copan_${datasetDate}.fam
sed -i                                             's/Argentina_BeagleChannel_Yamana_1500BP/Argentina_Yamana_1500BP/g'              Allen-1240-Copan_${datasetDate}.fam
sed -i                                              's/Argentina_MitrePeninsula_Haush_400BP/Argentina_Haush_400BP/g'                Allen-1240-Copan_${datasetDate}.fam
sed -i                                    's/Argentina_NorthTierradelFiego_Selknam_100BP.SG/Argentina_Selknam_100BP.SG/g'           Allen-1240-Copan_${datasetDate}.fam
sed -i                                 's/Argentina_NorthTierradelFiego_Selknam_100BP_lc.SG/Argentina_Selknam_100BP_SG/g'           Allen-1240-Copan_${datasetDate}.fam
sed -i                                 's/Argentina_NorthTierradelFuego_LaArcillosa2_5800BP/Argentina_LaArcillosa2_5800BP/g'        Allen-1240-Copan_${datasetDate}.fam
sed -i                                       's/Argentina_NorthTierradelFuego_Selknam_500BP/Argentina_Selknam_500BP/g'              Allen-1240-Copan_${datasetDate}.fam
sed -i                                            's/Argentina_Selknam_500BP_brother.I12367/Argentina_Selknam_500BP.I12367/g'       Allen-1240-Copan_${datasetDate}.fam #relative
sed -i                                           's/Bahamas_EleutheraIsl_Ceramic_dup.I14878/Bahamas_Ceramic_dup.I14878/g'           Allen-1240-Copan_${datasetDate}.fam
sed -i                                              's/Belize_8800BP_possible.1d.rel.I19170/Belize_8800BP_rel.I19170/g'             Allen-1240-Copan_${datasetDate}.fam
sed -i                                              's/Bolivia_MonolitoDescabezado_Tiwanaku/Bolivia_Tiwanaku/g'                     Allen-1240-Copan_${datasetDate}.fam
sed -i                                           's/Brazil_LapaDoSanto_9900BP_1d.rel.Lapa25/Brazil_LapaDoSanto_9900BP.Lapa25/g'     Allen-1240-Copan_${datasetDate}.fam
sed -i                                                              's/LateDorset-XIV-H_126/XIV-H_126/g'                            Allen-1240-Copan_${datasetDate}.fam
sed -i                                           's/Chile_NorthTierradelFuego_Selknam_100BP/Chile_NorthTierradelFuego_100BP/g'      Allen-1240-Copan_${datasetDate}.fam
sed -i                                            's/Chile_SouthernContinent_Aonikenk_400BP/Chile_Aonikenk_400BP/g'                 Allen-1240-Copan_${datasetDate}.fam
sed -i                                          's/Chile_StraitOfMagellan_Kaweskar_100BP.SG/Chile_Kaweskar_100BP.SG/g'              Allen-1240-Copan_${datasetDate}.fam
sed -i                                       's/Chile_StraitOfMagellan_Kaweskar_100BP_lc.SG/Chile_Kaweskar_100BP_lc.SG/g'           Allen-1240-Copan_${datasetDate}.fam
sed -i                                        's/Chile_StraitOfMagellan_Kaweskar_100BP_o.SG/Chile_Kaweskar_100BP_o.SG/g'            Allen-1240-Copan_${datasetDate}.fam
sed -i                                        's/Chile_WesternArchipelago_Ayayema_5100BP.SG/Chile_Ayayema_5100BP.SG/g'              Allen-1240-Copan_${datasetDate}.fam
sed -i                                        's/Chile_WesternArchipelago_Kaweskar_800BP.SG/Chile_Kaweskar_800BP.SG/g'              Allen-1240-Copan_${datasetDate}.fam
sed -i                                       's/Chile_WesternArchipelago_Kaweskar_1200BP.SG/Chile_Kaweskar_1200BP.SG/g'             Allen-1240-Copan_${datasetDate}.fam
sed -i                                               's/Chile_Yamana_BeagleChannel_800BP.SG/Chile_Yamana_800BP.SG/g'                Allen-1240-Copan_${datasetDate}.fam
sed -i                                    's/Cuba_CanimarAbajo_Archaic_dup.or.1d.rel.CAO015/Cuba_Archaic_dup-rel.CAO015/g'          Allen-1240-Copan_${datasetDate}.fam
sed -i                                  's/Cuba_CanimarAbajo_Archaic_o2_possible.dup.CAO032/Cuba_Archaic_dup.CAO032/g'              Allen-1240-Copan_${datasetDate}.fam
sed -i                                          's/Dominican_Atajadizo_Ceramic_1d.rel.17903/Dominican_Ceramic_rel.17903/g'          Allen-1240-Copan_${datasetDate}.fam
sed -i                                             's/Dominican_ElSoco_Ceramic_sister.13189/Dominican_Ceramic_rel.13189/g'          Allen-1240-Copan_${datasetDate}.fam
sed -i                                  's/Dominican_JuanDolio_Ceramic_father.or.son.I23524/Dominican_Ceramic_rel.I23524/g'         Allen-1240-Copan_${datasetDate}.fam
sed -i                                          's/Dominican_LaCaleta_Ceramic_1d.rel.I27938/Dominican_Ceramic_rel.I27938/g'         Allen-1240-Copan_${datasetDate}.fam
sed -i                                         's/Dominican_LaCaleta_Ceramic_brother.I27322/Dominican_Ceramic_rel.I27322/g'         Allen-1240-Copan_${datasetDate}.fam
sed -i            's/Dominican_LaCaleta_Ceramic_daughter.I25656.sister.I27329.sister.I15590/Dominican.I25656.I27329.I15590/g'       Allen-1240-Copan_${datasetDate}.fam #relative
sed -i 's/Dominican_LaCaleta_Ceramic_father.I25694.son.I25656.brother.I16687.brother.I27329/Dominican25694.25656.16687.27329/g'     Allen-1240-Copan_${datasetDate}.fam #relative
sed -i                            's/Dominican_LaCaleta_Ceramic_mother.I27316.sister.I23544/Dominican_Cer_rel.I27316.I23544/g'      Allen-1240-Copan_${datasetDate}.fam
sed -i                                          's/Dominican_LaCaleta_Ceramic_sister.I27301/Dominican_Ceramic_rel.I27301/g'         Allen-1240-Copan_${datasetDate}.fam
sed -i                                          's/Dominican_LaCaleta_Ceramic_sister.I27309/Dominican_Ceramic_rel.I27309/g'         Allen-1240-Copan_${datasetDate}.fam
sed -i                                         's/Dominican_Macao_Ceramic_1.or.2d.rel.I7973/Dominican_Ceramic_rel.I7973/g'          Allen-1240-Copan_${datasetDate}.fam
sed -i                                         's/Greenland_EarlyNorse_1d.rel.VK1.SG.VK9.SG/Greenland_EarlyNorse_rel.VK1.VK9/g'     Allen-1240-Copan_${datasetDate}.fam
sed -i                                        's/Greenland_EarlyNorse_1d.rel.VK11.SG.VK1.SG/Greenland_EarlyNorse_rel.VK11.VK1/g'    Allen-1240-Copan_${datasetDate}.fam
sed -i                                                    's/Ignore_Karitiana(discovery).DG/Ignore_Karitiana.DG/g'                  Allen-1240-Copan_${datasetDate}.fam
sed -i                                                   's/Ignore_Karitiana(discovery).SDG/Ignore_Karitiana.SDG/g'                 Allen-1240-Copan_${datasetDate}.fam
sed -i                                                        's/Ignore_Karitiana_discovery/Ignore_Karitiana/g'                     Allen-1240-Copan_${datasetDate}.fam
sed -i                                     's/Panama_IsthmoColombian_Colonial_oAfrica_lc.SG/Panama_Colonial.SG/g'                   Allen-1240-Copan_${datasetDate}.fam
sed -i                                  's/Panama_IsthmoColombian_Colonial_oEastAsian_lc.SG/Panama_Colonial.SG/g'                   Allen-1240-Copan_${datasetDate}.fam
sed -i                                   's/Panama_IsthmoColombian_Colonial_oEuropean_lc.SG/Panama_Colonial.SG/g'                   Allen-1240-Copan_${datasetDate}.fam
sed -i                             's/Panama_IsthmoColombian_Colonial_oNativeAmerican_lc.SG/Panama_Colonial.SG/g'                   Allen-1240-Copan_${datasetDate}.fam
sed -i                                            's/Panama_IsthmoColombian_oEuropean_lc.SG/Panama_IsthmoColombian.SG/g'            Allen-1240-Copan_${datasetDate}.fam
sed -i                                             's/Panama_IsthmoColombian_PreColonial.SG/Panama_PreColonial.SG/g'                Allen-1240-Copan_${datasetDate}.fam
sed -i                                                      's/Peru_Cuncaicha_3300BP_contam/Peru_3300BP_contam/g'                   Allen-1240-Copan_${datasetDate}.fam
sed -i                                        's/PuertoRico_CanasColloresMonserrate_Ceramic/PuertoRico_Ceramic/g'                   Allen-1240-Copan_${datasetDate}.fam
sed -i                                                   's/PuertoRico_PasodelIndio_Ceramic/PuertoRico_Ceramic/g'                   Allen-1240-Copan_${datasetDate}.fam
sed -i                                                 's/PuertoRico_PuntaCandelero_Ceramic/PuertoRico_Ceramic/g'                   Allen-1240-Copan_${datasetDate}.fam
sed -i                                               's/USA_AK_Ancient_Athabaskan_1100BP.DG/USA_Athabaskan_1100BP.DG/g'             Allen-1240-Copan_${datasetDate}.fam
sed -i                                               's/USA_AK_Ancient_Athabaskan_1100BP.SG/USA_Athabaskan_1100BP.SG/g'             Allen-1240-Copan_${datasetDate}.fam
sed -i                              's/USA_AK_Ancient_Athabaskan_1100BP_father.or.son.I5319/USA_Athabaskan_1100BP_rel.I5319/g'      Allen-1240-Copan_${datasetDate}.fam
sed -i                                                   's/USA_Alaska_TrailCreek_9000BP.SG/USA_TrailCreek_9000BP.SG/g'             Allen-1240-Copan_${datasetDate}.fam
sed -i                                         's/USA_CA_Chumash_mother.or.daughter.PS18.SG/USA_Chumash_rel.PS18.SG/g'              Allen-1240-Copan_${datasetDate}.fam
sed -i                                             's/USA_CA_Late_SanNicolas.SG_1d.rel.SN50/USA_Late_SanNicolas.SG_rel.SN50/g'      Allen-1240-Copan_${datasetDate}.fam
sed -i                                             's/USA_CA_Late_SanNicolas_1d.rel.SN09.SG/USA_Late_SanNicolas_rel.SN09.SG/g'      Allen-1240-Copan_${datasetDate}.fam
sed -i                                              's/USA_CA_Mainland_Chumash_contam_lc.SG/USA_Mainland_Chumash_contam.SG/g'       Allen-1240-Copan_${datasetDate}.fam
sed -i                                's/USA_MarianaIslands_Latte_brother.I24499.son.I24297/USA_Latte_rel.I24499.I24297/g'          Allen-1240-Copan_${datasetDate}.fam
sed -i                                           's/USA_MarianaIslands_Latte_brother.I24502/USA_Latte_rel.I24502/g'                 Allen-1240-Copan_${datasetDate}.fam
sed -i                                           's/USA_MarianaIslands_Latte_brother.I24606/USA_Latte_rel.I24606/g'                 Allen-1240-Copan_${datasetDate}.fam
sed -i                              's/USA_MarianaIslands_Latte_mother.I24499.mother.I24503/USA_Latte_rel.I24499.I24503/g'          Allen-1240-Copan_${datasetDate}.fam
sed -i                                            's/USA_MarianaIslands_Latte_mother.I24593/USA_Latte_rel.I24593/g'                 Allen-1240-Copan_${datasetDate}.fam
sed -i                                                  's/USA_Nevada_LovelockCave_600BP.SG/USA_LovelockCave_600BP.SG/g'            Allen-1240-Copan_${datasetDate}.fam
sed -i                                                 's/USA_Nevada_LovelockCave_1850BP.SG/USA_LovelockCave_1850BP.SG/g'           Allen-1240-Copan_${datasetDate}.fam
sed -i                                                  's/USA_Nevada_SpiritCave_11000BP.SG/USA_SpiritCave_11000BP.SG/g'            Allen-1240-Copan_${datasetDate}.fam
sed -i        's/USA_NM_Chaco_mother.I2717_mother.I6215.with.different.father.sibling.I5597/USA_Chaco_rel.I2717.I6215.I5597/g'      Allen-1240-Copan_${datasetDate}.fam


# Clean up directory
rm -r CpM* file2
rm Allen-1240-Copan_${datasetDate}.Allen-1240_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg*

exit
##############################
# HO dataset has less sites than 1240 and no additional individuals for the Americas
if [ ! -f v52.2_HO_public.anno ]; then
        wget https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V52/V52.2/SHARE/public.dir/v52.2_HO_public.tar
        tar -xf v52.2_HO_public.tar
        rm v52.2_HO_public.tar
fi

# Create plink pileup files: format by excluding chromosomes 23 and 24 and re-organise columns
awk '{ if($2<23) printf("%d %d %d\n",                $2,$4-1,$4)       ;}' v52.2_HO_public.snp > Allen-HO_${datasetDate}.pileup.bed
awk '{ if($2<23) printf("%d\t%d:%d\t0\t%d\t%s\t%s\n",$2,$2,$4,$4,$5,$6);}' v52.2_HO_public.snp > Allen-HO_${datasetDate}.pileup.bim
tail -n+2 v52.2_HO_public.anno    | awk -F'\t' '{print $13,$3,"0","0","9","1"}'                > Allen-HO_${datasetDate}.pileup.fam

# Pileup Copan individual with Allen dataset
python /Programs/HapSNP_v9_pg4.py \
                -data /st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/Allen/Allen-HO_${datasetDate} \
                -I copanBamInput \
                -R /Reference_Genomes/hs37d5/hs37d5.fa \
                -p 7 \
                -o Allen-HO-Copan_${datasetDate} \
                -mq 0 \
                -minbp 0 \
                -mem 250 \
                -gatk /Software/GenomeAnalysisTK.jar

# Create binary bed file
echo -e "genotypename: v52.2_HO_public.geno"                   > parfile_convertf
echo -e "snpname: v52.2_HO_public.snp"                        >> parfile_convertf
echo -e "indivname: v52.2_HO_public.ind"                      >> parfile_convertf
echo -e "outputformat: PACKEDPED"                             >> parfile_convertf
echo -e "genotypeoutname: Allen-HO_${datasetDate}.pileup.bed" >> parfile_convertf
echo -e "snpoutname: Allen-HO_${datasetDate}.bim"             >> parfile_convertf
echo -e "indivoutname: Allen-HO_${datasetDate}.fam"           >> parfile_convertf
convertf -p parfile_convertf
rm Allen-HO_${datasetDate}.bim Allen-HO_${datasetDate}.fam parfile_convertf

# Assign information to Copan samples
echo "Copan CpM0203 0 0 2 1"  > Allen-HO-Copan_${datasetDate}.Allen-HO_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM0506 0 0 2 1" >> Allen-HO-Copan_${datasetDate}.Allen-HO_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM11 0 0 1 1"   >> Allen-HO-Copan_${datasetDate}.Allen-HO_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM12 0 0 1 1"   >> Allen-HO-Copan_${datasetDate}.Allen-HO_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM13 0 0 1 1"   >> Allen-HO-Copan_${datasetDate}.Allen-HO_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM34 0 0 2 1"   >> Allen-HO-Copan_${datasetDate}.Allen-HO_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam
echo "Copan CpM4243 0 0 2 1" >> Allen-HO-Copan_${datasetDate}.Allen-HO_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg.fam

# Merge Copan with HOA and published shotgun (different bed file than before)
# and fiter for American populations because so many names to changes later on
grep "ACB\|AfricanAmericans\|Argentina\|ASW\|Aymara\|Bahamas\|Mexico\|Belize\|Bolivia\|Brazil\|Canada\|Chane\|Chile\|Chipewyan\|CLM\|Cree\|Cuba\|Curacao\|Dominican\|GIH\|Greenland\|Guadeloupe\|Haiti\|Hawaiian\|Huichol\|Karitiana\|Piapoco\|Mayan\|Mixe\|Mixtec\|MXL\|Nahua\|Panama\|PEL\|Peru\|Piapoco\|Pima\|PuertoRico\|PUR\|Quechua\|StLucia\|Surui\|Tsimshian\|USA\|Venezuela\|Yukpa" Allen-HO-Copan_${datasetDate}.fam | awk '{print $1,$2}' > americans.pop
echo "Allen-HO_${datasetDate}.pileup"   > file2
plink --bfile Allen-HO-Copan_${datasetDate}.Allen-HO_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg \
        --merge-list file2 --make-bed --allow-no-sex \
        --keep americans.pop \
        --out Allen-HO-Copan_${datasetDate}

# Clean up directory
rm -r CpM* file2
rm Allen-1240-Copan_${datasetDate}.Allen-1240_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg*
rm Allen-HO-Copan_${datasetDate}.Allen-HO_${datasetDate}.hs37d5.bq30.mq0.0bp.homozyg*
