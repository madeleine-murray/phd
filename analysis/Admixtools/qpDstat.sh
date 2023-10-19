#!/bin/bash

DATE=010123
popdir="Copan_qpDstats/"
set1='Copan Mixe'
set2='Copan SpiritCave'
set3='Copan ASO_4000'
set4='Copan Sumidouro'
set5='Copan Anzick'
set6='Copan Ayayema'
set7='Copan PuntaSantaAna'
set8='Copan USR1'
set9='Copan Mayan'
set10='Copan Taino'
F4="YES"
out="Mbuti"

population=`awk '{print $1}' ../SGDP/SGDP.database.Copan.${DATE}.fam | uniq | tr "\n" " "`


mkdir -p ${popdir}
for k in {1..10}; do

	setvar=\$set${k}
	set=`eval echo ${setvar}`
	target=`echo ${set} | awk '{print $1}'`
	ref=`echo ${set}    | awk '{print $2}'`

	mkdir -p /${popdir}${out}.${ref}.${target}
	cd /${popdir}${out}.${ref}.${target}

	if test "${F4}" != "YES"; then
		rm -f par.${out}.${ref}.${target}.Dstats.result
	else
		rm -f par.${out}.${ref}.${target}.Fstats.result
	fi
	
	for pop in ${population}; do
		echo "${out} ${pop} ${ref} ${target}" > ${out}.${pop}.${ref}.${target}.list
	
		echo -e "genotypename: ../../../SGDP/SGDP.database.Copan.${DATE}.eigenstratgeno" > par.${out}.${pop}.${ref}.${target}.Dstats
		echo -e "snpname: ../../../SGDP/SGDP.database.Copan.${DATE}.snp"                >> par.${out}.${pop}.${ref}.${target}.Dstats
		echo -e "indivname: ../../../SGDP/SGDP.database.Copan.${DATE}.ind"              >> par.${out}.${pop}.${ref}.${target}.Dstats
		echo -e "popfilename: ${out}.${pop}.${ref}.${target}.list"                      >> par.${out}.${pop}.${ref}.${target}.Dstats
		echo -e "printsd: YES"                                                          >> par.${out}.${pop}.${ref}.${target}.Dstats
		echo -e "f4mode: ${F4}"                                                         >> par.${out}.${pop}.${ref}.${target}.Dstats

	done
done
exit
		if test "${F4}" != "YES"; then
			/Software/AdmixTools-v6/bin/qpDstat -p par.${out}.${pop}.${ref}.${target}.Dstats > par.${out}.${pop}.${ref}.${target}.Dstats.log
			grep "result:" par.${out}.${pop}.${ref}.${target}.Dstats.log \
				| awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11;}' - >> \
				par.${out}.${ref}.${target}.Dstats.result
		else
			/Software/AdmixTools-v6/bin/qpDstat -p par.${out}.${pop}.${ref}.${target}.Dstats > par.${out}.${pop}.${ref}.${target}.Fstats.log
			grep "result:" par.${out}.${pop}.${ref}.${target}.Fstats.log \
				| awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11;}' - >> \
				par.${out}.${ref}.${target}.Fstats.result
		fi
	done
done
