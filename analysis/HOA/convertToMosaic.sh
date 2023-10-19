#!/bin/bash

date=07122022 # keep track of pileup version

##############
# PHASE DATA #
##############

mkdir -p mosaic_input/
mv HOA_${date}.pileup.bim HOA_${date}.bim
mv HOA_${date}.pileup.fam HOA_${date}.fam
awk '{ if ($6!="-9") print $0}' HOA_${date}.fam > mosaic_input/include.txt
plink --bfile HOA_${date} --keep mosaic_input/include.txt --recode vcf --out mosaic_input/HOA_${date} # convert plink to vcf format
mv HOA_${date}.bim HOA_${date}.pileup.bim
mv HOA_${date}.fam HOA_${date}.pileup.fam
rm mosaic_input/include.txt

bcftools norm -m - -Ou --threads 4 -o mosaic_input/HOA_${date}.preprocessed.vcf mosaic_input/HOA_${date}.vcf # pre-processing for shapeit
bgzip mosaic_input/HOA_${date}.preprocessed.vcf # required for separating per chr
tabix -p vcf mosaic_input/HOA_${date}.preprocessed.vcf.gz # required for separating per chr

for chr in {1..22}; do

	bcftools view -r $chr -Ov -o mosaic_input/HOA_${date}_chr${chr}.vcf mosaic_input/HOA_${date}.preprocessed.vcf.gz # separate per chr for phasing (otherwise sorting issue)

	shapeit -V mosaic_input/HOA_${date}_chr${chr}.vcf -O mosaic_input/HOA_${date}_chr${chr}_phased --output-log mosaic_input/shapeit.log & # phase data
done

rm mosaic_input/HOA_${date}_chr*.vcf
rm mosaic_input/HOA_${date}.log mosaic_input/HOA_${date}.nosex HOA_${date}.vcf
rm mosaic_input/HOA_${date}.preprocessed.vcf.gz mosaic_input/HOA_${date}.preprocessed.vcf.gz.tbi

#####################
# CONVERT TO MOSAIC #
#####################


for chr in {1..22}; do
	
	echo "Creating MOSAIC input for chromosome ${chr}..." 

	# Create rate file
	if [ ! -f mosaic_input/rates.${chr} ]; then
		if [ ! -f 1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt ]; then
			wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz
			tar zxvf 1000GP_Phase3.tgz
		fi
		nbSites=`tail -n +2 1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt | wc -l`
		echo ":sites:"${nbSites} > mosaic_input/rates.${chr}
		tail -n +2 1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt | awk '{print $1}' | paste -sd" " >> mosaic_input/rates.${chr}
		tail -n +2 1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt | awk '{print $3}' | paste -sd" " >> mosaic_input/rates.${chr}
	fi
	
	# Create snp file
	if [ ! -f mosaic_input/snpfile.${chr} ]; then
		awk '{print $2,$1,0,$3,$4,$5}' mosaic_input/HOA_${date}_chr${chr}_phased.haps > mosaic_input/snpfile.${chr}
	fi

	# Create geno file
	group=`awk '{ if ($6!="-9") print $1}' HOA_${date}.pileup.fam | uniq`
	beginning=6
	finish=5
	for pop in ${group}; do
		nbPop=`awk '{ if ($1=="'${pop}'") print $0}' HOA_${date}.pileup.fam | wc -l`
		finish=$((2*$nbPop + $finish))
		awk 'BEGIN {first='$beginning'; last='$finish'} {for (i=first; i<last; i++) {printf("%s ", $i)} print $last}' mosaic_input/HOA_${date}_chr${chr}_phased.haps \
			| sed 's/ *//g' > mosaic_input/${pop}genofile.${chr}
		beginning=`expr $finish + 1`
	done

	# Merge geno files per continent
	paste -d '' mosaic_input/AA_Denvergenofile.${chr} \
		mosaic_input/Algeriangenofile.${chr} \
		mosaic_input/Bantu_SA_Hererogenofile.${chr} \
		mosaic_input/Bantu_SA_Ovambogenofile.${chr} \
		mosaic_input/Bantu_SA_Pedigenofile.${chr} \
		mosaic_input/Bantu_SA_S_Sothogenofile.${chr} \
		mosaic_input/Bantu_SA_Tswanagenofile.${chr} \
		mosaic_input/Bantu_SA_Zulugenofile.${chr} \
		mosaic_input/BantuKenyagenofile.${chr} \
		mosaic_input/BiakaPygmygenofile.${chr} \
		mosaic_input/Datoggenofile.${chr} \
		mosaic_input/Egyptian_Comasgenofile.${chr} \
		mosaic_input/Egyptian_Metspalugenofile.${chr} \
		mosaic_input/Esan_Nigeria_ESNgenofile.${chr} \
		mosaic_input/Ethiopian_Jewgenofile.${chr} \
		mosaic_input/Gambian_GWDgenofile.${chr} \
		mosaic_input/Hadza_Henngenofile.${chr} \
		mosaic_input/Ju_hoan_Northgenofile.${chr} \
		mosaic_input/Khomanigenofile.${chr} \
		mosaic_input/Kikuyugenofile.${chr} \
		mosaic_input/Libyan_Jewgenofile.${chr} \
		mosaic_input/Luhya_Kenya_LWKgenofile.${chr} \
		mosaic_input/Luogenofile.${chr} \
		mosaic_input/Mandenkagenofile.${chr} \
		mosaic_input/Masai_Ayodogenofile.${chr} \
		mosaic_input/Masai_Kinyawa_MKKgenofile.${chr} \
		mosaic_input/MbutiPygmygenofile.${chr} \
		mosaic_input/Mende_Sierra_Leone_MSLgenofile.${chr} \
		mosaic_input/Moroccan_Jewgenofile.${chr} \
		mosaic_input/Mozabitegenofile.${chr} \
		mosaic_input/Saharawigenofile.${chr} \
		mosaic_input/Somaligenofile.${chr} \
		mosaic_input/Tunisiangenofile.${chr} \
		mosaic_input/Tunisian_Jewgenofile.${chr} \
		mosaic_input/Yorubagenofile.${chr} \
		> mosaic_input/Africagenofile.${chr}
	paste -d '' mosaic_input/Bolivian_Cochabambagenofile.${chr} \
		mosaic_input/Bolivian_LaPazgenofile.${chr} \
		mosaic_input/Bolivian_Pandogenofile.${chr} \
		mosaic_input/Karitianagenofile.${chr} \
		mosaic_input/Mixegenofile.${chr} \
		mosaic_input/Mixtecgenofile.${chr} \
		mosaic_input/Piapocogenofile.${chr} \
		mosaic_input/Pimagenofile.${chr} \
		mosaic_input/Quechua_Coriellgenofile.${chr} \
		mosaic_input/Suruigenofile.${chr} \
		mosaic_input/Zapotecgenofile.${chr} \
		> mosaic_input/Americagenofile.${chr}
	paste -d '' mosaic_input/Aleutgenofile.${chr} \
		mosaic_input/Altaiangenofile.${chr} \
		mosaic_input/Chukchigenofile.${chr} \
		mosaic_input/Chukchi_Reindeergenofile.${chr} \
		mosaic_input/Chukchi_Sirgenofile.${chr} \
		mosaic_input/Dolgangenofile.${chr} \
		mosaic_input/Eskimo_Chaplingenofile.${chr} \
		mosaic_input/Eskimo_Naukangenofile.${chr} \
		mosaic_input/Eskimo_Sirenikigenofile.${chr} \
		mosaic_input/Evengenofile.${chr} \
		mosaic_input/Itelmengenofile.${chr} \
		mosaic_input/Kalmykgenofile.${chr} \
		mosaic_input/Koryakgenofile.${chr} \
		mosaic_input/Kyrgyzgenofile.${chr} \
		mosaic_input/Mansigenofile.${chr} \
		mosaic_input/Mongolagenofile.${chr} \
		mosaic_input/Nganasangenofile.${chr} \
		mosaic_input/Selkupgenofile.${chr} \
		mosaic_input/Tlingitgenofile.${chr} \
		mosaic_input/Tubalargenofile.${chr} \
		mosaic_input/Tuviniangenofile.${chr} \
		mosaic_input/Ulchigenofile.${chr} \
		mosaic_input/Yakutgenofile.${chr} \
		mosaic_input/Yukagir_Forestgenofile.${chr} \
		mosaic_input/Yukagir_Tundragenofile.${chr} \
		> mosaic_input/Siberiagenofile.${chr}
	paste -d '' mosaic_input/Ami_Coriellgenofile.${chr} \
		mosaic_input/Atayal_Coriellgenofile.${chr} \
		mosaic_input/Cambodiangenofile.${chr} \
		mosaic_input/Daigenofile.${chr} \
		mosaic_input/Daurgenofile.${chr} \
		mosaic_input/Hangenofile.${chr} \
		mosaic_input/Han_NChinagenofile.${chr} \
		mosaic_input/Hezhengenofile.${chr} \
		mosaic_input/Japanesegenofile.${chr} \
		mosaic_input/Kinh_Vietnam_KHVgenofile.${chr} \
		mosaic_input/Koreangenofile.${chr} \
		mosaic_input/Lahugenofile.${chr} \
		mosaic_input/Miaogenofile.${chr} \
		mosaic_input/Naxigenofile.${chr} \
		mosaic_input/Oroqengenofile.${chr} \
		mosaic_input/Shegenofile.${chr} \
		mosaic_input/Thaigenofile.${chr} \
		mosaic_input/Tugenofile.${chr} \
		mosaic_input/Tujiagenofile.${chr} \
		mosaic_input/Uygurgenofile.${chr} \
		mosaic_input/Xibogenofile.${chr} \
		mosaic_input/Yigenofile.${chr} \
		> mosaic_input/EastAsiagenofile.${chr}
	paste -d '' mosaic_input/Australian_ECCACgenofile.${chr} \
		mosaic_input/Bougainvillegenofile.${chr} \
		mosaic_input/Papuangenofile.${chr} \
		> mosaic_input/Oceaniagenofile.${chr}
	paste -d '' mosaic_input/Balochigenofile.${chr} \
		mosaic_input/Bengali_Bangladesh_BEBgenofile.${chr} \
		mosaic_input/Brahuigenofile.${chr} \
		mosaic_input/Burushogenofile.${chr} \
		mosaic_input/Cochin_Jewgenofile.${chr} \
		mosaic_input/GujaratiA_GIHgenofile.${chr} \
		mosaic_input/GujaratiB_GIHgenofile.${chr} \
		mosaic_input/GujaratiC_GIHgenofile.${chr} \
		mosaic_input/GujaratiD_GIHgenofile.${chr} \
		mosaic_input/Hazaragenofile.${chr} \
		mosaic_input/Kalashgenofile.${chr} \
		mosaic_input/Kusundagenofile.${chr} \
		mosaic_input/Makranigenofile.${chr} \
		mosaic_input/Pathangenofile.${chr} \
		mosaic_input/Punjabi_Lahore_PJLgenofile.${chr} \
		mosaic_input/Sindhigenofile.${chr} \
		> mosaic_input/SouthAsiagenofile.${chr}
	paste -d '' mosaic_input/Albaniangenofile.${chr} \
		mosaic_input/Ashkenazi_Jewgenofile.${chr} \
		mosaic_input/Basque_Frenchgenofile.${chr} \
		mosaic_input/Basque_Spanishgenofile.${chr} \
		mosaic_input/Belarusiangenofile.${chr} \
		mosaic_input/Bulgariangenofile.${chr} \
		mosaic_input/Croatiangenofile.${chr} \
		mosaic_input/Cypriotgenofile.${chr} \
		mosaic_input/Czechgenofile.${chr} \
		mosaic_input/English_Cornwall_GBRgenofile.${chr} \
		mosaic_input/English_Kent_GBRgenofile.${chr} \
		mosaic_input/Estoniangenofile.${chr} \
		mosaic_input/Finnish_FINgenofile.${chr} \
		mosaic_input/Frenchgenofile.${chr} \
		mosaic_input/French_Southgenofile.${chr} \
		mosaic_input/Greek_Comasgenofile.${chr} \
		mosaic_input/Greek_Coriellgenofile.${chr} \
		mosaic_input/Hungarian_Coriellgenofile.${chr} \
		mosaic_input/Hungarian_Metspalugenofile.${chr} \
		mosaic_input/Icelandicgenofile.${chr} \
		mosaic_input/Italian_Bergamogenofile.${chr} \
		mosaic_input/Italian_EastSiciliangenofile.${chr} \
		mosaic_input/Italian_Southgenofile.${chr} \
		mosaic_input/Italian_Tuscangenofile.${chr} \
		mosaic_input/Italian_WestSiciliangenofile.${chr} \
		mosaic_input/Lithuaniangenofile.${chr} \
		mosaic_input/Maltesegenofile.${chr} \
		mosaic_input/Mordoviangenofile.${chr} \
		mosaic_input/Norwegiangenofile.${chr} \
		mosaic_input/Orcadiangenofile.${chr} \
		mosaic_input/Russiangenofile.${chr} \
		mosaic_input/Saami_WGAgenofile.${chr} \
		mosaic_input/Sardiniangenofile.${chr} \
		mosaic_input/Scottish_Argyll_Bute_GBRgenofile.${chr} \
		mosaic_input/Spanish_Andalucia_IBSgenofile.${chr} \
		mosaic_input/Spanish_Aragon_IBSgenofile.${chr} \
		mosaic_input/Spanish_Baleares_IBSgenofile.${chr} \
		mosaic_input/Spanish_Canarias_IBSgenofile.${chr} \
		mosaic_input/Spanish_Cantabria_IBSgenofile.${chr} \
		mosaic_input/Spanish_Castilla_la_Mancha_IBSgenofile.${chr} \
		mosaic_input/Spanish_Castilla_y_Leon_IBSgenofile.${chr} \
		mosaic_input/Spanish_Cataluna_IBSgenofile.${chr} \
		mosaic_input/Spanish_Extremadura_IBSgenofile.${chr} \
		mosaic_input/Spanish_Galicia_IBSgenofile.${chr} \
		mosaic_input/Spanish_Murcia_IBSgenofile.${chr} \
		mosaic_input/Spanish_Pais_Vasco_IBSgenofile.${chr} \
		mosaic_input/Spanish_Valencia_IBSgenofile.${chr} \
		mosaic_input/Swedish_Motalagenofile.${chr} \
		mosaic_input/Ukrainian_Eastgenofile.${chr} \
		mosaic_input/Ukrainian_Westgenofile.${chr} \
		> mosaic_input/WestEurasiagenofile.${chr}
	paste -d '' mosaic_input/Abkhasiangenofile.${chr} \
		mosaic_input/Adygeigenofile.${chr} \
		mosaic_input/Armeniangenofile.${chr} \
		mosaic_input/Balkargenofile.${chr} \
		mosaic_input/BedouinAgenofile.${chr} \
		mosaic_input/BedouinBgenofile.${chr} \
		mosaic_input/Chechengenofile.${chr} \
		mosaic_input/Chuvashgenofile.${chr} \
		mosaic_input/Druzegenofile.${chr} \
		mosaic_input/Georgian_Jewgenofile.${chr} \
		mosaic_input/Georgian_Megrelsgenofile.${chr} \
		mosaic_input/Iraniangenofile.${chr} \
		mosaic_input/Iranian_Jewgenofile.${chr} \
		mosaic_input/Iraqi_Jewgenofile.${chr} \
		mosaic_input/Jordaniangenofile.${chr} \
		mosaic_input/Kumykgenofile.${chr} \
		mosaic_input/Lebanesegenofile.${chr} \
		mosaic_input/Lezgingenofile.${chr} \
		mosaic_input/Nogaigenofile.${chr} \
		mosaic_input/North_Ossetiangenofile.${chr} \
		mosaic_input/Palestiniangenofile.${chr} \
		mosaic_input/Saudigenofile.${chr} \
		mosaic_input/Syriangenofile.${chr} \
		mosaic_input/Tajik_Pomirigenofile.${chr} \
		mosaic_input/Turkishgenofile.${chr} \
		mosaic_input/Turkish_Adanagenofile.${chr} \
		mosaic_input/Turkish_Aydingenofile.${chr} \
		mosaic_input/Turkish_Balikesirgenofile.${chr} \
		mosaic_input/Turkish_Istanbulgenofile.${chr} \
		mosaic_input/Turkish_Jewgenofile.${chr} \
		mosaic_input/Turkish_Kayserigenofile.${chr} \
		mosaic_input/Turkish_Trabzongenofile.${chr} \
		mosaic_input/Turkmengenofile.${chr} \
		mosaic_input/Uzbekgenofile.${chr} \
		mosaic_input/Yemengenofile.${chr} \
		mosaic_input/Yemenite_Jewgenofile.${chr} \
		> mosaic_input/MiddleEastgenofile.${chr}

done
./assignContinent.sh # script to automise continent assignment for HOA moderns
awk '{ if ($2!="ancient") print $2}' HOA.pop > mosaic_input/sample.names
if [ -f 1000GP_Phase3 ]; then rm -r 1000GP_Phase3; fi
if [ -f 1000GP_Phase3.tgz ]; then rm 1000GP_Phase3.tgz; fi

group=`awk '{ if ($6!="-9") print $1}' HOA_${date}.pileup.fam | uniq` # delete temporary genofiles
for pop in ${group}; do rm mosaic_input/${pop}genofile*; done         # to save memory on disc

exit
########################################################################################
# NOTE:                                                                                #
# Difference between HOA-publishedShotgun-copan_${date}.fam and HOA_${date}.pileup.fam #
# are 7 Copan samples pileup up                                                        #
# minus 2 ancients which has 100% missing genotypes                                    #
########################################################################################

conda activate localAncestry
Rscript /home/madeleine/R/x86_64-pc-linux-gnu-library/4.1/MOSAIC/mosaic.R \
	Mayan mosaic_input/ \
	--ancestries 3 \
	--number 18 \
	--chromosomes 1:22 \
	--maxcores 6 \
	--panels "Africa America Siberia EastAsia SouthAsia WestEurasia"
conda deactivate
