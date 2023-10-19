#!/bin/bash


# Input
DATE=010123                                                                                        # date for reference pileup
source[0]="Mayan"                                                                                  # source for combinations for left pop
source[1]="ancientPanama"                                                                          # source for combinations for left pop
source[2]="TIW"                                                                                    # source for combinations for left pop
source[3]="LSN"                                                                                    # source for combunations for left pop
source[4]="COR"                                                                                    # source for ocmbinations for left pop
right[0]="Mbuti Sardinian Han Sumidouro Ayayema PuntaSantaAna ASO_4000 MA1 SpiritCave Anzick Mixe" # version 1 for right pop
right[1]="Mbuti Sardinian Han Ayayema PuntaSantaAna ASO_4000 MA1 SpiritCave Anzick Mixe"           # version 2 for right pop


mkdir -p Copan_qpAdm
for ((k=0;k<${#right[@]};k++)); do
	for ((i=0;i<${#source[@]};i++)); do
        	for ((j=0;j<${#source[@]};j++)); do
                	if (($i != $j)); then

				# Create left pops
                        	mkdir -p Copan_qpAdm/Copan_${source[i]}_${source[j]}
                        	echo "Copan ${source[i]} ${source[j]}" | sed s/" "/"\n"/g > Copan_qpAdm/Copan_${source[i]}_${source[j]}/left.pops

				# Create right pops
				echo ${right[k]} | sed s/" "/"\n"/g > Copan_qpAdm/right_${k}.pops
	
				# Create parameter file
				echo -e "genotypename: ../SGDP/SGDP.database.Copan.${DATE}.eigenstratgeno" \
					> Copan_qpAdm/Copan_${source[i]}_${source[j]}/par.qpAdm.right_${k}
				echo -e "snpname: ../SGDP/SGDP.database.Copan.${DATE}.snp"                 \
					>> Copan_qpAdm/Copan_${source[i]}_${source[j]}/par.qpAdm.right_${k}
				echo -e "indivname: ../SGDP/SGDP.database.Copan.${DATE}.ind"               \
					>> Copan_qpAdm/Copan_${source[i]}_${source[j]}/par.qpAdm.right_${k}
				echo -e "popleft: Copan_qpAdm/Copan_${source[i]}_${source[j]}/left.pops"   \
					>> Copan_qpAdm/Copan_${source[i]}_${source[j]}/par.qpAdm.right_${k}
				echo -e "popright: Copan_qpAdm/right_${k}.pops"                            \
					>> Copan_qpAdm/Copan_${source[i]}_${source[j]}/par.qpAdm.right_${k}
				echo -e "details: YES"                                                     \
					>> Copan_qpAdm/Copan_${source[i]}_${source[j]}/par.qpAdm.right_${k}
				echo -e "allsnps: YES"                                                     \
					>> Copan_qpAdm/Copan_${source[i]}_${source[j]}/par.qpAdm.right_${k}

				# Run program
				/Software/AdmixTools-v6/bin/qpAdm -p Copan_qpAdm/Copan_${source[i]}_${source[j]}/par.qpAdm.right_${k} \
					> Copan_qpAdm/output.qpAdm.right_${k}.Copan.${source[i]}.${source[j]}
			fi
		done
	done
done
exit
grep "00  0\|01  1\|10  1" Copan_qpAdm/output.qpAdm.right* | column -t
