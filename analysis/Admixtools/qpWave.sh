#!/bin/bash

DATE=010123 # date for reference pileup

# According the f4 statistics, Maya are not enough to explain ancestry in Copan, so test for admixture between "Copan/Mayan/source"=left vs right poppulations
source[0]="" # testing if additional pop source better explains the model
source[1]="Lovelock"
source[2]="Island_Chumash"
source[3]="Taino"
source[4]="PuntaSantaAna"
source[5]="ESN"
source[6]="Aconcagua"
source[7]="COR"
source[8]="ancientPanama"
source[9]="LUK"
source[10]="Pericu"
source[11]="Fuego-Patagonian"
source[12]="Rio_Uncallane"
source[13]="Mainland_Chumash"
source[14]="Baja_Mexico"
source[15]="TIW"
source[16]="LSN"
source[17]="SMP"
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

mkdir -p Copan_qpWave
for ((k=0;k<${#right[@]};k++)); do
	for ((i=0;i<${#source[@]};i++)); do

		# Create left pops
		mkdir -p Copan_qpWave/Copan_Mayan_${source[i]}
                echo "Copan Mayan ${source[i]}" | sed s/" "/"\n"/g > Copan_qpWave/Copan_Mayan_${source[i]}/left.pops

		# Create right pops
		echo ${right[k]} | sed s/" "/"\n"/g > Copan_qpWave/right_${k}.pops
	
		# Create parameter file
		echo -e "genotypename: ../SGDP/SGDP.database.Copan.${DATE}.eigenstratgeno" \
			> Copan_qpWave/Copan_Mayan_${source[i]}/par.qpWave.right_${k}
		echo -e "snpname: ../SGDP/SGDP.database.Copan.${DATE}.snp"                 \
			>> Copan_qpWave/Copan_Mayan_${source[i]}/par.qpWave.right_${k}
		echo -e "indivname: ../SGDP/SGDP.database.Copan.${DATE}.ind"               \
			>> Copan_qpWave/Copan_Mayan_${source[i]}/par.qpWave.right_${k}
		echo -e "popleft: Copan_qpWave/Copan_Mayan_${source[i]}/left.pops"         \
			>> Copan_qpWave/Copan_Mayan_${source[i]}/par.qpWave.right_${k}
		echo -e "popright: Copan_qpWave/right_${k}.pops"                           \
			>> Copan_qpWave/Copan_Mayan_${source[i]}/par.qpWave.right_${k}
		echo -e "details: YES"                                                     \
			>> Copan_qpWave/Copan_Mayan_${source[i]}/par.qpWave.right_${k}
		echo -e "allsnps: YES"                                                     \
			>> Copan_qpWave/Copan_Mayan_${source[i]}/par.qpWave.right_${k}

		# Run program
		/Software/AdmixTools-v6/bin/qpWave -p Copan_qpWave/Copan_Mayan_${source[i]}/par.qpWave.right_${k} \
			> Copan_qpWave/output.qpWave.right_${k}.Copan.Mayan.${source[i]}
	done
done
exit
grep "00  0\|01  1\|10  1" Copan_qpWave/output.qpWave.right* | column -t
