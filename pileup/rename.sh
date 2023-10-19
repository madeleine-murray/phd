#!/bin/bash


directory="shotgun/Raghaven_et_al_2015_Recent_history_native_americans/alignments/hs37d5/final_files/"


####################                   SET UP FILE NAMES                   ####################
list=`ls -lS ${directory}*sorted.grouped.duprm.indel.bam | awk '{print $9}'` # order list by file size so that next step starts only when biggest file has finished
sample=()
for file in $list; do
        sample+=$(echo $file | rev | cut -f7-20 -d. | rev); # get sample name from file
        sample+=' '                                         # add space to seperate elements in list
done
array=($sample)                                             # convert list to array (for iteration later-on)
###############################################################################################

for ((i=0; i<${#array[@]}; i++)); do
	mv ${array[i]}.trimmedbp25q.sorted.grouped.duprm.indel.mq25.clip2bp.34bp_stats ${array[i]}.mixedTrimmed.sorted.grouped.duprm.indel.mq25.clip2bp.34bp_stats
	mv ${array[i]}.trimmedbp25q.sorted.grouped.duprm.indel.bai ${array[i]}.mixedTrimmed.sorted.grouped.duprm.indel.bai
	mv ${array[i]}.trimmedbp25q.sorted.grouped.duprm.indel.bam ${array[i]}.mixedTrimmed.sorted.grouped.duprm.indel.bam
	mv ${array[i]}.trimmedbp25q.sorted.grouped.duprm.indel.intervals ${array[i]}.mixedTrimmed.sorted.grouped.duprm.indel.intervals
	mv ${array[i]}.trimmedbp25q.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam ${array[i]}.mixedTrimmed.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam
	mv ${array[i]}.trimmedbp25q.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam.bai ${array[i]}.mixedTrimmed.sorted.grouped.duprm.indel.mq25.clip2bp.34bp.bam.bai
done
