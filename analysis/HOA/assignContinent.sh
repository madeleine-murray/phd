#!/bin/bash

###########################################
# WARNING: IN THE CONTEXT OF OUR ANALYSIS #
# MAYAN WERE ASSIGNED AS A SEPARATE GROUP #
###########################################

awk '{ if ($6=="1") print $1}' HOA_07122022.pileup.fam > HOA.pop # copy moderns only

awk -i inplace '$1 == "Mayan" { $2 = "Mayan" } 1' HOA.pop # group of interest

Africa="AA_Denver Algerian Bantu_SA_Herero Bantu_SA_Ovambo Bantu_SA_Pedi Bantu_SA_S_Sotho Bantu_SA_Tswana Bantu_SA_Zulu BantuKenya BiakaPygmy Datog Egyptian_Comas Egyptian_Metspalu Esan_Nigeria_ESN Ethiopian_Jew Gambian_GWD Hadza_Henn Ju_hoan_North Khomani Kikuyu Libyan_Jew Luhya_Kenya_LWK Luo Mandenka Masai_Ayodo Masai_Kinyawa_MKK MbutiPygmy Mende_Sierra_Leone_MSL Moroccan_Jew Mozabite Saharawi Somali Tunisian Tunisian_Jew Yoruba"
for group in ${Africa}; do
	awk -i inplace '$1 == "'${group}'" { $2 = "Africa" } 1' HOA.pop
done

America="Bolivian_Cochabamba Bolivian_LaPaz Bolivian_Pando Karitiana Mixe Mixtec Piapoco Pima Quechua_Coriell Surui Zapotec"
for group in ${America}; do
	awk -i inplace '$1 == "'${group}'" { $2 = "America" } 1' HOA.pop
done

Siberia="Aleut Altaian Chukchi Chukchi_Reindeer Chukchi_Sir Dolgan Eskimo_Chaplin Eskimo_Naukan Eskimo_Sireniki Even Itelmen Kalmyk Koryak Kyrgyz Mansi Mongola Nganasan Selkup Tlingit Tubalar Tuvinian Ulchi Yakut Yukagir_Forest Yukagir_Tundra"
for group in ${Siberia}; do
	awk -i inplace '$1 == "'${group}'" { $2 = "Siberia" } 1' HOA.pop
done

EastAsia="Ami_Coriell Atayal_Coriell Cambodian Dai Daur Han Han_NChina Hezhen Japanese Kinh_Vietnam_KHV Korean Lahu Miao Naxi Oroqen She Thai Tu Tujia Uygur Xibo Yi"
for group in ${EastAsia}; do
	awk -i inplace '$1 == "'${group}'" { $2 = "EastAsia" } 1' HOA.pop
done

Oceania="Australian_ECCAC Bougainville Papuan"
for group in ${Oceania}; do
	awk -i inplace '$1 == "'${group}'" { $2 = "Oceania" } 1' HOA.pop
done

SouthAsia="Balochi Bengali_Bangladesh_BEB Brahui Burusho Cochin_Jew GujaratiA_GIH GujaratiB_GIH GujaratiC_GIH GujaratiD_GIH Hazara Kalash Kusunda Makrani Pathan Punjabi_Lahore_PJL Sindhi"
for group in ${SouthAsia}; do
	awk -i inplace '$1 == "'${group}'" { $2 = "SouthAsia" } 1' HOA.pop
done

Colonist="English_Cornwall_GBR English_Kent_GBR French French_South Scottish_Argyll_Bute_GBR Spanish_Andalucia_IBS Spanish_Aragon_IBS Spanish_Baleares_IBS Spanish_Canarias_IBS Spanish_Cantabria_IBS Spanish_Castilla_la_Mancha_IBS Spanish_Castilla_y_Leon_IBS Spanish_Cataluna_IBS Spanish_Extremadura_IBS Spanish_Galicia_IBS Spanish_Murcia_IBS Spanish_Pais_Vasco_IBS Spanish_Valencia_IBS"
for group in ${Colonist}; do
	awk -i inplace '$1 == "'${group}'" { $2 = "Colonist" } 1' HOA.pop
done

WestEurasia="Albanian Ashkenazi_Jew Basque_French Basque_Spanish Belarusian Bulgarian Croatian Cypriot Czech Estonian Finnish_FIN Greek_Comas Greek_Coriell Hungarian_Coriell Hungarian_Metspalu Icelandic Italian_Bergamo Italian_EastSicilian Italian_South Italian_Tuscan Italian_WestSicilian Lithuanian Maltese Mordovian Norwegian Orcadian Russian Saami_WGA Sardinian Swedish_Motala Ukrainian_East Ukrainian_West"
for group in ${WestEurasia}; do
	awk -i inplace '$1 == "'${group}'" { $2 = "WestEurasia" } 1' HOA.pop
done

MiddleEast="Abkhasian Adygei Armenian Balkar BedouinA BedouinB Chechen Chuvash Druze Georgian_Jew Georgian_Megrels Iranian Iranian_Jew Iraqi_Jew Jordanian Kumyk Lebanese Lezgin Nogai North_Ossetian Palestinian Saudi Syrian Tajik_Pomiri Turkish Turkish_Adana Turkish_Aydin Turkish_Balikesir Turkish_Istanbul Turkish_Jew Turkish_Kayseri Turkish_Trabzon Turkmen Uzbek Yemen Yemenite_Jew"
for group in ${MiddleEast}; do
	awk -i inplace '$1 == "'${group}'" { $2 = "MiddleEast" } 1' HOA.pop
done

awk -i inplace '{ if ($2 == "") $2 = "ancient" } 1' HOA.pop
