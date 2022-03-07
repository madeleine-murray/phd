#!/bin/bash

##################################################################################################################################################################
                                                        # DESCRIPTION #                                                    #                       #             #
                                      # After the testing imputation accuracy of JpKa6904                                  #                     ###             #
                                          # from script glimpseScript_testAccuracy.sh                                      #                 #####               #
                                           # (that only calculated for transversions)                                      #               ##   # ##             #
                                            # here we calculate for transitions only                                       #              ##   #   ##            #
                                        # (the previous script must have been run before)                                  #              ##  #    ##            #
                                                                                                                           #               ###    ##             #
                                                                                                                           #               # #####               #
                                                       # PARAMETERS #                                                      #              #                      #
CHROM=("1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22")  # list of chromosomes to impute upon, autosomes only   #######################################
RUN=/Software/bcftools-1.10.2/                                      # software path                                                                              #
dir=/st1/ssd/shigeki/glimpse                                        # working directory                                                                          #
DOWNSAMPLE=("0X01 0X05 0X1 0X5 1X0 2X0 3X0 4X0 5X0")                # downsample information                                                                     #
sample=Yamnaya                                                      # name of sample ID                                                                          #
name=sorted.grouped.mq25.rmdup.indel.softclip.rl34                  # useless but hopefully make the script easier to read                                       #
DIP=${sample}_1000G_GENOTYPE_30_DEPT_10                             # prefix of original diploid file (bed, bim and fam)                                         #
# corrected diploid: awk '{ $2 = ($1":"$4); print $1,$2,$3,$4,$5,$6}' F23_GENOTYPE_30_DEPT_10.bim > F23_1000G_GENOTYPE_30_DEPT_10.bim                            #
REFALT=${sample}_1000G_GENOTYPE_30_DEPT_10                          # information for alternative and reference alleles                                          #
mafRefFile=/st1/hdd/pg/ancient_human_resources/reference_datasets/modern_wgs_studies/1000G_Phase_3/plink/whole_genome/MAF_filters/1000G_autosomes.worldmaf01.bim #
threshold=99                                                        # imputed snp certainty threshold in pourcentage                                             #
output=${dir}/finalAccuracy/${sample}                               # name of output file, tab seperated                                                         #
FILTER=("trv trs trv.maf01 trs.maf01")                              # filtering type: transversions/transitions, maf                                             #
##################################################################################################################################################################


# CHECKING FOR POTENTIAL ERRORS

# Check that original diploid file is present (from step 9)
mkdir -p ${dir}/glimpse_plink
if [ ! -f ${dir}/glimpse_plink/${DIP}.bed ]; then
        echo "Reference diploid file is missing: bed, bim and fam must be present."
        echo "${DIP}"
        echo "Please add files to ${dir}/glimpse_plink/"
        exit
fi

# Check for the maf reference file (from step 10)
for filter in ${FILTER}; do
        if [ ${filter} = "trv.maf01" ] || [ ${filter} = "trs.maf01" ]; then
                if [ ! -f ${mafRefFile} ]; then
                        echo "The file reference for MAF is not found."
                        echo "Please enter a valid path for 'mafRefFile' in the parameters."
                        exit
                fi
        fi
done

# Extract information for alternative and reference alleles (step 11)
mkdir -p ${dir}/glimpse_vcf
mkdir -p ${dir}/glimpse_comparison
for CHR in ${CHROM}; do
        if [ ! -f ${dir}/glimpse_comparison/${REFALT}.chr${CHR}.refalt ]; then
                if [ ! -f ${dir}/glimpse_vcf/${REFALT}.vcf.gz ]; then
                        echo "Information for alterative and reference alleles is missing."
                        echo "Please create the following file:"
                        echo "${dir}/glimpse_vcf/${REFALT}.vcf.gz"
                        exit
                else
                        zcat ${dir}/glimpse_vcf/${REFALT}.vcf.gz           \
                                | awk '{if($1=="chr'${CHR}'") printf("'${CHR}':%d %c %c\n",$2,$4,$5);}' \
                                > ${dir}/glimpse_comparison/${REFALT}.chr${CHR}.refalt;
                fi
        fi
done


# STEP 9: EXTRACT COMMON SNP
for downsample in ${DOWNSAMPLE}; do
        for CHR in ${CHROM}; do

                # Extract the second column (i.e. CHR:POS) from the diploid bim
                awk '{if($1=="'${CHR}'")print $2;}' ${dir}/glimpse_plink/${DIP}.bim > ${dir}/glimpse_comparison/aaa

                # Extract SNP sites in "aaa" from the imputed data
                plink --bfile ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup \
                        --noweb --extract ${dir}/glimpse_comparison/aaa \
                        --out ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip --make-bed

                # Extract the second column (i.e. CHR:POS) from imputed JpKa6904
                cut -f 2 ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.bim \
                        > ${dir}/glimpse_comparison/bbb
                # Extract SNP sites in "bbb" from the diploid data
                plink --bfile ${dir}/glimpse_plink/${DIP} \
                        --noweb --extract ${dir}/glimpse_comparison/bbb \
                        --out ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp --make-bed

        done
done
rm ${dir}/glimpse_comparison/aaa
rm ${dir}/glimpse_comparison/bbb

# STEP 10: CONVERT TO TPED WITH TRANSVERIONS OR TRANSITIONS, WITH OR WITHOUT MAF FILTERING
# This analysis only focuses on the transversion SNP sites to minimize the impact of DNA damages in ancient DNA data
for downsample in ${DOWNSAMPLE}; do
        for CHR in ${CHROM}; do
                for filter in ${FILTER}; do

                        if [ ${filter} = "trv.maf01" ] || [ ${filter} = "trs.maf01" ]; then
                                cut -f 2 ${mafRefFile} > ${dir}/glimpse_comparison/snp.1000G.maf01.list

                                plink --bfile ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp \
                                        --noweb --extract ${dir}/glimpse_comparison/snp.1000G.maf01.list \
                                        --make-bed --out ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.maf
                                plink --bfile ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip \
                                        --noweb --extract ${dir}/glimpse_comparison/snp.1000G.maf01.list \
                                        --make-bed --out ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.maf
                                awk '{if(($5=="G" && $6=="A") || ($5=="A" && $6=="G") || ($5=="C" && $6=="T") || ($5=="T" && $6=="C"))print;}' \
                                        ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.maf.bim \
                                        | cut -f 2 > ${dir}/glimpse_comparison/trans.list
                        else
                                awk '{if(($5=="G" && $6=="A") || ($5=="A" && $6=="G") || ($5=="C" && $6=="T") || ($5=="T" && $6=="C"))print;}' \
                                        ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.bim \
                                        | cut -f 2 > ${dir}/glimpse_comparison/trans.list
                        fi

                        if [ ${filter} = "trv" ] || [ ${filter} = "trv.maf01" ]; then
                                if [ ${filter} = "trv.maf01" ]; then
                                plink --bfile ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.maf \
                                        --exclude ${dir}/glimpse_comparison/trans.list \
                                        --recode --transpose \
                                        --out ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.${filter}
                                plink --bfile ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.maf \
                                        --exclude ${dir}/glimpse_comparison/trans.list \
                                        --recode --transpose \
                                        --out ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.${filter}
                                else
                                plink --bfile ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp \
                                        --exclude ${dir}/glimpse_comparison/trans.list \
                                        --recode --transpose \
                                        --out ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.${filter}
                                plink --bfile ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip \
                                        --exclude ${dir}/glimpse_comparison/trans.list \
                                        --recode --transpose \
                                        --out ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.${filter}
                                fi
                                elif [ ${filter} = "trs" ] || [ ${filter} = "trs.maf01" ]; then
                                if [ ${filter} = "trs.maf01" ]; then
                                plink --bfile ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.maf \
                                        --extract ${dir}/glimpse_comparison/trans.list \
                                        --recode --transpose \
                                        --out ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.${filter}
                                plink --bfile ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.maf \
                                        --extract ${dir}/glimpse_comparison/trans.list \
                                        --recode --transpose \
                                        --out ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.${filter}
                                else
                                plink --bfile ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp \
                                        --extract ${dir}/glimpse_comparison/trans.list \
                                        --recode --transpose \
                                        --out ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.${filter}
                                plink --bfile ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip \
                                        --extract ${dir}/glimpse_comparison/trans.list \
                                        --recode --transpose \
                                        --out ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.${filter}
                                fi
                        else
                                echo "Error: Please do not change the filtering parameters"
                                echo "You can remove the undesired outcomes, but do no add your 'inventions' "
                                rm ${dir}/glimpse_comparison/trans.list
                                exit
                        fi
                done
        done
done
rm ${dir}/glimpse_comparison/trans.list


# STEP 11: COMPARE AND OUTPUT RESULTS
mkdir -p ${dir}/testAccuracy
for filter in ${FILTER}; do
        if [ ! -f ${output}.${filter}.tsv ]; then
                echo "Downsample"       "Chr"   "hetero_rate"   "homo_ref_rate" "homo_alt_rate" "hetero_count"  "homo_ref_count"        "homo_alt_count"        "hetero_total"  "homo_ref_total"        "homo_alt_total" > ${output}.${filter}.tsv
        fi # create header
done

for downsample in ${DOWNSAMPLE}; do
        for CHR in ${CHROM}; do
                for filter in ${FILTER}; do

                        # STEP 11.1: COMPARE WITH NON-IMPUTED DATA

                        allSNP=`paste ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.${filter}.tped \
                                ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.${filter}.tped \
                                | awk -F' ' '{if(($5==$11 && $6==$12) || ($5==$12 && $6==$11))print;}' | wc -l`

                        # Heterozygote
                        heteroSite=`paste ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.${filter}.tped \
                                ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.${filter}.tped \
                                | awk -F' ' '{if($11!=$12 && (($5==$11 && $6==$12) || ($5==$12 && $6==$11))) print;}' | wc -l`
                        heteroTotal=`awk -F' ' '{if($5!=$6)print;}' ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.${filter}.tped \
                                | wc -l`

                        # Homozygote: retrieve reference and alternative information and get SNPs present in reference only
                        cut -d " " -f 2 ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.${filter}.tped > ${dir}/glimpse_comparison/new.list
                        fgrep -wf ${dir}/glimpse_comparison/new.list ${dir}/glimpse_comparison/${REFALT}.chr${CHR}.refalt \
                                > ${dir}/glimpse_comparison/${DIP}.imp.${filter}.refalt
                        cut -d " " -f 1 ${dir}/glimpse_comparison/${DIP}.imp.${filter}.refalt > ${dir}/glimpse_comparison/new_list
                        fgrep -wf ${dir}/glimpse_comparison/new_list \
                                ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.${filter}.tped \
                                > ${dir}/glimpse_comparison/diploid_common_with_imputed_and_reference.tped
                        fgrep -wf ${dir}/glimpse_comparison/new_list \
                                ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.${filter}.tped \
                                > ${dir}/glimpse_comparison/imputed_common_with_diploid_and_reference.tped
                        rm ${dir}/glimpse_comparison/new.list
                        rm ${dir}/glimpse_comparison/new_list

                        homoRefSite=`paste -d " " ${dir}/glimpse_comparison/${DIP}.imp.${filter}.refalt \
                                ${dir}/glimpse_comparison/diploid_common_with_imputed_and_reference.tped \
                                ${dir}/glimpse_comparison/imputed_common_with_diploid_and_reference.tped \
                                | awk -F' ' '{if($14==$15 && $2==$14 && (($8==$14 && $9==$15) || ($8==$15 && $9==$14)))print;}' | wc -l`
                        homoRefTotal=`paste -d " " ${dir}/glimpse_comparison/${DIP}.imp.${filter}.refalt \
                                ${dir}/glimpse_comparison/imputed_common_with_diploid_and_reference.tped \
                                | awk -F' ' '{if($8==$9 && $2==$8)print;}' | wc -l`
                        homoAltSite=`paste -d " " ${dir}/glimpse_comparison/${DIP}.imp.${filter}.refalt \
                                ${dir}/glimpse_comparison/diploid_common_with_imputed_and_reference.tped \
                                ${dir}/glimpse_comparison/imputed_common_with_diploid_and_reference.tped \
                                | awk -F' ' '{if($14==$15 && $3==$14 && (($8==$14 && $9==$15) || ($8==$15 && $9==$14)))print;}' | wc -l`
                        homoAltTotal=`paste -d " " ${dir}/glimpse_comparison/${DIP}.imp.${filter}.refalt \
                                ${dir}/glimpse_comparison/imputed_common_with_diploid_and_reference.tped \
                                | awk -F' ' '{if($8==$9 && $3==$8)print;}' | wc -l`

                        rm ${dir}/glimpse_comparison/diploid_common_with_imputed_and_reference.tped
                        rm ${dir}/glimpse_comparison/imputed_common_with_diploid_and_reference.tped
                        
                        # STEP 11.2: CALCULATE CONCORDANCE RATES AND OUTPUT

                        heteroRate=`echo  "scale=2; ${heteroSite}  * 100 / ${heteroTotal}"  | bc`
                        homoRefRate=`echo "scale=2; ${homoRefSite} * 100 / ${homoRefTotal}" | bc`
                        homoAltRate=`echo "scale=2; ${homoAltSite} * 100 / ${homoAltTotal}" | bc`

                        echo "${downsample}     ${CHR}  ${heteroRate}   ${homoRefRate}  ${homoAltRate}  ${heteroSite}   ${homoRefSite}  ${homoAltSite}  ${heteroTotal}  ${homoRefTotal} ${homoAltTotal}" >> ${output}.${filter}.tsv

                        #rm ${dir}/glimpse_comparison/${DIP}.imp.${filter}.refalt
                done
        done
done
exit

# STEP 12: REMOVE USELESS FILES
# I seperated each file so that it is easy to tweak
for downsample in ${DOWNSAMPLE}; do
        for CHR in ${CHROM}; do
                rm -f ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.bed
                rm -f ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.bim
                rm -f ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.fam
                rm -f ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.log
                rm -f ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.nosex

                rm -f ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.maf.bed
                rm -f ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.maf.bim
                rm -f ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.maf.fam
                rm -f ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.maf.log
                rm -f ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.maf.nosex

                rm -f ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.bed
                rm -f ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.bim
                rm -f ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.fam
                rm -f ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.log
                rm -f ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.nosex

                rm -f ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.maf.bed
                rm -f ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.maf.bim
                rm -f ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.maf.fam
                rm -f ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.maf.log
                rm -f ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.maf.nosex
                
                for filter in ${FILTER}; do
                        rm ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.${filter}.tped
                        rm ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.${filter}.tfam
                        rm ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.${filter}.log
                        rm ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.${filter}.nosex

                        rm ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.${filter}.tped
                        rm ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.${filter}.tfam
                        rm ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.${filter}.log
                        rm ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.${filter}.nosex
                done
        done
done
rm ${dir}/glimpse_comparison/${REFALT}.chr*.refalt
for filter in ${FILTER}; do
        echo "Your results have been outputed to ${output}.${filter}.tsv"
done
echo "DONE: all files have been processed."
