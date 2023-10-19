#!/bin/bash

#DATE="Copan.010123" # date for reference pileup
#DATE="7CpM.010123"
DATE="CpM13.010123.masked"
#copan="CpM0203 CpM0506 CpM11 CpM12 CpM13 CpM34 CpM4243"
copan="Copan CpM13"
subdir="Copan_qpAdm2_CpM13_masked"
mkdir -p ${subdir}/

# Convert to eigenstrat
#sed 's/Copan CpM13 0 0 0 1/CpM13 CpM13 0 0 0 1/' ../SGDP/SGDP.database.Copan.010123.fam > ../SGDP/SGDP.database.${DATE}.fam
#echo -e "genotypename: ../SGDP/SGDP.database.Copan.010123.bed"           > ../SGDP/parfile_convertf
#echo -e "snpname: ../SGDP/SGDP.database.Copan.010123.bim"               >> ../SGDP/parfile_convertf
#echo -e "indivname: ../SGDP/SGDP.database.${DATE}.fam"                  >> ../SGDP/parfile_convertf
#echo -e "outputformat: EIGENSTRAT"                                      >> ../SGDP/parfile_convertf
#echo -e "genotypeoutname: ../SGDP/SGDP.database.${DATE}.eigenstratgeno" >> ../SGDP/parfile_convertf
#echo -e "snpoutname: ../SGDP/SGDP.database.${DATE}.snp"                 >> ../SGDP/parfile_convertf
#echo -e "indivoutname: ../SGDP/SGDP.Database.${DATE}.ind"               >> ../SGDP/parfile_convertf
#convertf -p ../SGDP/parfile_convertf
#sed 's/:/\t/' ../SGDP/SGDP.Database.${datasetDate}.ind | awk '{print $2,$3,$1}' > ../SGDP/SGDP.database.${datasetDate}.ind
#rm ../SGDP/parfile_convertf ../SGDP/SGDP.Database.${datasetDate}.ind

# Convert masked to eignstrat
#sed 's/Copan CpM13 0 0 0 1/CpM13 CpM13 0 0 0 1/' ../SGDP/mask/SGDP.database.Copan.010123.masked.fam > ../SGDP/SGDP.database.${DATE}.fam
#echo -e "genotypename: ../SGDP/mask/SGDP.database.${DATE}.bed"           > ../SGDP/parfile_convertf
#echo -e "snpname: ../SGDP/mask/SGDP.database.${DATE}.bim"               >> ../SGDP/parfile_convertf
#echo -e "indivname: ../SGDP/SGDP.database.${DATE}.fam"                  >> ../SGDP/parfile_convertf
#echo -e "outputformat: EIGENSTRAT"                                      >> ../SGDP/parfile_convertf
#echo -e "genotypeoutname: ../SGDP/SGDP.database.${DATE}.eigenstratgeno" >> ../SGDP/parfile_convertf
#echo -e "snpoutname: ../SGDP/SGDP.database.${DATE}.snp"                 >> ../SGDP/parfile_convertf
#echo -e "indivoutname: ../SGDP/SGDP.Database.${DATE}.ind"               >> ../SGDP/parfile_convertf
#convertf -p ../SGDP/parfile_convertf
#sed 's/:/\t/' ../SGDP/SGDP.Database.${DATE}.ind | awk '{print $2,$3,$1}' > ../SGDP/SGDP.database.${DATE}.ind
#rm ../SGDP/parfile_convertf ../SGDP/SGDP.Database.${DATE}.ind

# According the f4 statistics, Maya are not enough to explain ancestry in Copan, so test for admixture between "Copan/Mayan/source"=left vs right populations
source[0]="" # testing if additional pop source better explains the model
#source[1]="Lovelock"
#source[2]="Island_Chumash"
#source[3]="Taino"
#source[4]="PuntaSantaAna"
#source[5]="ESN"
#source[6]="Aconcagua"
source[7]="COR"
#source[8]="ancientPanama"
source[9]="LUK"
#source[10]="Pericu"
#source[11]="Fuego-Patagonian"
source[12]="Rio_Uncallane"
source[13]="Mainland_Chumash"
#source[14]="Baja_Mexico"
source[15]="TIW"
source[16]="LSN"
#source[17]="SMP"
source[18]="ORU"
source[19]="Sumidouro"
source[20]="Kaillachuro"
source[21]="LSCI"
# (Ayayema is not included in left because in right)

# Shouldn't include admixed population
right[0]="Mbuti Sardinian Han Sumidouro Ayayema PuntaSantaAna ASO_4000 MA1 SpiritCave Anzick Mixe"      # version 1 for right pop (Sumidouro)
right[1]="Mbuti Sardinian Han Ayayema PuntaSantaAna ASO_4000 MA1 SpiritCave Anzick Mixe"                # version 2 for right pop
right[2]="Mbuti Sardinian Han Sumidouro Ayayema PuntaSantaAna ASO_4000 MA1 USR1 SpiritCave Anzick Mixe" # version 3 for right pop (Sumidouro and USR1)
right[3]="Mbuti Sardinian Han Ayayema PuntaSantaAna ASO_4000 MA1 USR1 SpiritCave Anzick Mixe"           # version 4 for right pop (USR1)

mkdir -p ${subdir}
for target in ${copan}; do
	for ((k=0;k<${#right[@]};k++)); do
		for ((i=0;i<${#source[@]};i++)); do

			# Create left pops
			mkdir -p ${subdir}/${target}_Mayan_${source[i]}
                	echo "${target} Mayan ${source[i]}" | sed s/" "/"\n"/g > ${subdir}/${target}_Mayan_${source[i]}/left.pops

			# Create right pops
			echo ${right[k]} | sed s/" "/"\n"/g > ${subdir}/right_${k}.pops
	
			# Create parameter file
			echo -e "genotypename: ../SGDP/SGDP.database.${DATE}.eigenstratgeno"  > ${subdir}/${target}_Mayan_${source[i]}/par.qpAdm.right_${k}
			echo -e "snpname: ../SGDP/SGDP.database.${DATE}.snp"                 >> ${subdir}/${target}_Mayan_${source[i]}/par.qpAdm.right_${k}
			echo -e "indivname: ../SGDP/SGDP.database.${DATE}.ind"               >> ${subdir}/${target}_Mayan_${source[i]}/par.qpAdm.right_${k}
			echo -e "popleft: ${subdir}/${target}_Mayan_${source[i]}/left.pops"  >> ${subdir}/${target}_Mayan_${source[i]}/par.qpAdm.right_${k}
			echo -e "popright: ${subdir}/right_${k}.pops"                        >> ${subdir}/${target}_Mayan_${source[i]}/par.qpAdm.right_${k}
			echo -e "details: YES"                                               >> ${subdir}/${target}_Mayan_${source[i]}/par.qpAdm.right_${k}
			echo -e "allsnps: YES"                                               >> ${subdir}/${target}_Mayan_${source[i]}/par.qpAdm.right_${k}

			# Run program
			/Software/AdmixTools-v6/bin/qpAdm -p ${subdir}/${target}_Mayan_${source[i]}/par.qpAdm.right_${k} \
				> ${subdir}/output.qpAdm.right_${k}.${target}.Mayan.${source[i]}
		done
	done
done
exit
grep "00  0\|01  1\|10  1" ${subdir}/output.qpAdm.right* | column -t

grep "std. errors:" output.qpAdm.right*Copan* | awk '{print $1}'
grep "best coefficients:" ${subdir}/output.qpAdm.right* | awk '{print $4}'
grep "std. errors:" output.qpAdm.right*Copan* | awk '{print $4}'
