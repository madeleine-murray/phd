#!/bin/bash

#############################################
#               DAMAGE PATTERN              #
# Access damage patterns on non-UDG treated #
#############################################

# Run map damage
for file in $(ls ../mergePlatform/CpM*Mi*bam); do
        sample=${file:17:5}
	mapDamage \
		-i ${file} \
		-d results_${sample} \
		-r /Reference_Genomes/hs37d5/hs37d5.fa \
		> ${sample}.mapDamage.log \
		2> ${sample}.mapDamage_error.log
	mv results_${sample}/3pGtoA_freq.txt               results_${sample}/${sample}_3pGtoA_freq.txt
	mv results_${sample}/5pCtoT_freq.txt               results_${sample}/${sample}_5pCtoT_freq.txt
	mv results_${sample}/dnacomp.txt                   results_${sample}/${sample}_dnacomp.txt
	mv results_${sample}/Fragmisincorporation_plot.pdf results_${sample}/${sample}_fragmisincorporation_plot.pdf
	mv results_${sample}/Length_plot.pdf               results_${sample}/${sample}_length_plot.pdf
	mv results_${sample}/lgdistribution.txt            results_${sample}/${sample}_lgdistribution.txt
	mv results_${sample}/misincorporation.txt          results_${sample}/${sample}_misincorporation.txt
	mv results_${sample}/Runtime_log.txt               results_${sample}/${sample}_runtime_log.txt
done
