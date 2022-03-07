#!/bin/bash

############################################################################################################################
                                                       # DESCRIPTION #                                                     #
                                   # Imputing downsampled versions of genomes with glimpse                                 #
                                   # To assess imputation accuracy in ancient individuals                                  #
                                  # https://odelaneau.github.io/GLIMPSE/tutorial_b38.html                                  #
                         # In this new version, we use of a mapping quality of 25 instead of 20                            #
               # And we re-set genotype likelihoods to zero for transitions to ignore deamination and DNA damage           #
                               # Using the script: /Programs/filter_vcf_for_imputation_v5.py                               #
                                                                                                                           #
                                                                                                                           #
                                                       # PARAMETERS #                                                      #
thread=15                                                           # number of threads for parallelisation                #
CHROM=("1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22")  # list of chromosomes to impute upon, autosomes only   #
RUN=/Software/bcftools-1.10.2/                                      # software path                                        #
dir=/st1/ssd/shigeki/glimpse                                        # working directory                                    #
DOWNSAMPLE=("0X01 0X05 0X1 0X5 1X0 2X0 3X0 4X0 5X0")                # downsample information                               #
sample=Yamnaya                                                      # name of sample ID                                    #  
name=sorted.grouped.mq25.rmdup.indel.softclip.rl34                  # useless but hopefully make the script easier to read #
BAM=${dir}/glimpse_bam/${sample}.${name}                            # downsampled JpKa6904                                 #
DIP=${sample}_1000G_GENOTYPE_30_DEPT_10                             # prefix of original diploid file (bed, bim and fam)   #
REFALT=${sample}_1000G_GENOTYPE_30_DEPT_10                          # information for alternative and reference alleles    #
REF=${dir}/glimpse_ref_panel/1000GP                                 # reference panel                                      #
threshold=99                                                        # imputed snp certainty threshold in pourcentage       #
output=${dir}/testAccuracy/${sample}_results_test.tsv               # name of output file, tab seperated                   #
############################################################################################################################


# CHECKING FOR POTENTIAL ERRORS

# Check that original diploid file is present (from step 9)
mkdir -p ${dir}/glimpse_plink
if [ ! -f ${dir}/glimpse_plink/${DIP}.bed ]; then
        echo "Reference diploid file is missing: bed, bim and fam must be present."
        echo "${DIP}"
        echo "Please add files to ${dir}/glimpse_plink/"
        exit
fi

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


# STEP 2.2: CREATE REFERENCE PANEL
mkdir -p ${dir}/glimpse_ref_panel
for CHR in ${CHROM}; do
        if [ ! -f ${REF}.chr${CHR}.bcf ] || [ ! -f ${REF}.chr${CHR}.bcf.csi ]; then
                for chr in {1..23} X; do
                        echo ${chr} chr${chr}
                done >> ${dir}/glimpse_ref_panel/chr_names.txt

                ${RUN}bcftools annotate --rename-chrs ${dir}/glimpse_ref_panel/chr_names.txt \
                        /st1/hdd/pg/ancient_human_resources/modern_genotype_datasets/1000G_Phase_3/vcfs/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -Ou | \
                        ${RUN}bcftools view -m 2 -M 2 -v snps --threads ${thread} -Ob -o ${REF}.chr${CHR}.bcf
                ${RUN}bcftools index -f ${REF}.chr${CHR}.bcf
        fi
done
rm -f ${dir}/glimpse_ref_panel/chr_names.txt # temporary file only for re-naming


# STEP 3: COMPUTATION OF GENOTYPE LIKELIHOOD
# GLIMPSE requires input data to take the form of Genotype Likelihoods (GLs).
# GLs need to be computed at all target individuals and all variant sites
# present in the reference panel of haplotypes used for the imputation.

# STEP 3.1: EXTRACTING VARIBLES IN THE REFERENCE PANEL
for CHR in ${CHROM}; do
        if [ ! -f ${REF}.chr${CHR}.sites.tsv.gz ] || [ ! -f ${REF}.chr${CHR}.sites.tsv.gz.tbi ]; then
                ${RUN}bcftools view -G -m 2 -M 2 -v snps ${REF}.chr${CHR}.bcf -Oz \
                        -o ${REF}.chr${CHR}.sites.vcf.gz
                ${RUN}bcftools index -f ${REF}.chr${CHR}.sites.vcf.gz
                ${RUN}bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ${REF}.chr${CHR}.sites.vcf.gz \
                        | bgzip -c > ${REF}.chr${CHR}.sites.tsv.gz
                tabix -s1 -b2 -e2 ${REF}.chr${CHR}.sites.tsv.gz
        fi
done

# STEP 3.2: COMPUTING GENOTYPE LIKELIHOODS FOR A SINGLE INDIVIDUAL
mkdir -p ${dir}/glimpse_bam
for downsample in ${DOWNSAMPLE}; do
        if   [ ! -f ${BAM}.${downsample}.REGROUPED.bam ]; then
                echo "ERROR: the input file is missing."
                echo "Please create the following file:"
                echo "${BAM}.${downsample}.REGROUPED.bam"
                exit
        elif [ ! -f ${BAM}.${downsample}.REGROUPED.bam.bai ]; then
                samtools index ${BAM}.${downsample}.REGROUPED.bam
        fi

        for CHR in ${CHROM}; do
                if [ ! -f ${dir}/glimpse_vcf/${sample}.${name}.${downsample}.chr${CHR}.REGROUPED.vcf.gz ]; then
                        ${RUN}bcftools mpileup -f /Reference_Genomes/hg19_Cristina_with_rCRS/hg19_Cristina_with_rCRS.fasta -I -E -a 'FORMAT/DP' \
                                -T ${REF}.chr${CHR}.sites.vcf.gz \
                                -r chr${CHR} ${BAM}.${downsample}.REGROUPED.bam -Ou | \
                                ${RUN}bcftools call -Aim -C alleles \
                                -T ${REF}.chr${CHR}.sites.tsv.gz -Oz \
                                -o ${dir}/glimpse_vcf/${sample}.${name}.${downsample}.chr${CHR}.REGROUPED.vcf.gz
                fi

                if [ ! -f ${dir}/glimpse_vcf/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.vcf.gz.csi ]; then
                        gunzip -f ${dir}/glimpse_vcf/${sample}.${name}.${downsample}.chr${CHR}.REGROUPED.vcf.gz
                        python /Programs/filter_vcf_for_imputation_v6.py ${dir}/glimpse_vcf/${sample}.${name}.${downsample}.chr${CHR}.REGROUPED.vcf \
                                | bgzip -c > ${dir}/glimpse_vcf/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.vcf.gz
                        bcftools index -f ${dir}/glimpse_vcf/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.vcf.gz
                fi
        done
done

# STEP 3.3: MERGING GENOTYPE LIKELIHOOD OF MULTIPLE INDIVIDUALS
#bcftools merge -m none -r ${CHR} -Oz -o merged.chr${CHR}.${downsample}.vcf.gz -l list.txt


# STEP 4: SPLIT THE GENOME INTO CHUNKS
for CHR in ${CHROM}; do
        if [ ! -f ${dir}/glimpse_ref_panel/chunks.chr${CHR}.txt ]; then
                GLIMPSE_chunk --input ${REF}.chr${CHR}.sites.vcf.gz \
                        --region chr${CHR} --window-size 2000000 --buffer-size 200000 \
                        --output ${dir}/glimpse_ref_panel/chunks.chr${CHR}.txt
        fi
done


# STEP 5: IMPUTE AND PHASE
mkdir -p ${dir}/glimpse_imputed

for downsample in ${DOWNSAMPLE}; do
        for CHR in ${CHROM}; do

                # Extract variable from chunk file
                IDvar=`cut -f 1 ${dir}/glimpse_ref_panel/chunks.chr${CHR}.txt`
                IRGvar=`cut -f 3 ${dir}/glimpse_ref_panel/chunks.chr${CHR}.txt`
                ORGvar=`cut -f 4 ${dir}/glimpse_ref_panel/chunks.chr${CHR}.txt`

                # Convert string variable to list (for looping afterwards)
                IDlistA=($IDvar); IDlistB=($IDvar); IRGlist=($IRGvar); ORGlist=($ORGvar);
                for ((i=0;i<10;i++)); do IDlistB[$i]="0"${IDlistB[$i]}; done # convert the single digits to double digits

                # Parallelise the corresponding number of threads
                # (If ID is not dividable by the number of threads or equal to zero, then run in background, else run in foregound:
                # this allows to run jobs on a controlled number threads. And wait for the last runs to finish.)
                for i in ${IDlistA[@]}; do
                        if [ $(expr $i % $thread ) -ne 0 ] || [ $i -eq 0 ]; then
                                GLIMPSE_phase --input ${dir}/glimpse_vcf/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.vcf.gz \
                                        --reference ${REF}.chr${CHR}.bcf \
                                        --map ${dir}/glimpse_map/chr${CHR}.b37.gmap.gz \
                                        --input-region ${IRGlist[$i]} --output-region ${ORGlist[$i]} \
                                        --output ${dir}/glimpse_imputed/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.${IDlistB[$i]}.bcf &
                                # wait for the last chunk to finish before indexing
                                if [ $i -eq ${IDlistA[-1]} ]; then
                                        wait
                                fi
                        else
                                GLIMPSE_phase --input ${dir}/glimpse_vcf/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.vcf.gz \
                                        --reference ${REF}.chr${CHR}.bcf \
                                        --map ${dir}/glimpse_map/chr${CHR}.b37.gmap.gz \
                                        --input-region ${IRGlist[$i]} --output-region ${ORGlist[$i]} \
                                        --output ${dir}/glimpse_imputed/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.${IDlistB[$i]}.bcf
                        fi
                done

                # Index only after all imputation has finished
                for i in ${IDlistA[@]}; do
                        ${RUN}bcftools index -f ${dir}/glimpse_imputed/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.${IDlistB[$i]}.bcf
                done
        done
done

# STEP 6: LIGATE MULTIPLE CHUNKS
mkdir -p ${dir}/glimpse_ligated
for downsample in ${DOWNSAMPLE}; do
        for CHR in ${CHROM}; do
                ls ${dir}/glimpse_imputed/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.*.bcf > ${dir}/glimpse_ligated/list.chr${CHR}.txt
                GLIMPSE_ligate --input ${dir}/glimpse_ligated/list.chr${CHR}.txt \
                        --output ${dir}/glimpse_ligated/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.bcf
                ${RUN}bcftools index -f ${dir}/glimpse_ligated/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.bcf
        done
done


# STEP 7: CONVERT TO VCF
for downsample in ${DOWNSAMPLE}; do
        for CHR in ${CHROM}; do
                # STEP 7.1: CONVERT BCF TO VCF
                ${RUN}bcftools view -Ov \
                        -o ${dir}/glimpse_ligated/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.vcf \
                        ${dir}/glimpse_ligated/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.bcf
                ${RUN}bcftools index -f ${dir}/glimpse_ligated/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.vcf

                # STEP 7.2: FILTER IMPUTED SNP FOR 99% CERTAINTY
                python /Programs/filter_vcf_GP_v2.py \
                        -i ${dir}/glimpse_ligated/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.vcf \
                        -o ${dir}/glimpse_ligated/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.vcf \
                        -t 0.${threshold}
        done
done


# STEP 8: CONVERT TO PLINK

# STEP 8.1: CONVERT VCF TO PLINK
for downsample in ${DOWNSAMPLE}; do
        for CHR in ${CHROM}; do
                plink --vcf ${dir}/glimpse_ligated/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.vcf \
                        --make-bed \
                        --out ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold} \
                        --geno 0
        done
done

# STEP 8.2: CORRECT FORMAT OF SNP ID COLUMN (CHR:POS)
for downsample in ${DOWNSAMPLE}; do
        for CHR in ${CHROM}; do

                # Save original imputed file
                if [ ! -f ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.bim.org ]; then
                        mv ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.bim \
                                ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.bim.org
                        awk '{printf("%d\t%d:%d\t0\t%d\t%c\t%c\n",$1,$1,$4,$4,$5,$6);}' \
                                ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.bim.org \
                                > ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.bim
                                # correct format of SNP ID column (CHR:POS)
                fi

                # Remove duplicated SNP sites
                cut -f 2 ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.bim \
                        | sort | uniq -d > ${dir}/glimpse_plink/snp.remove
                plink --bfile ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold} \
                        --noweb --exclude ${dir}/glimpse_plink/snp.remove \
                        --out ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup --make-bed
                rm ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.bed
                rm ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.bim
                rm ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.fam
                rm ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.log
                rm ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nosex
        done
done
rm ${dir}/glimpse_plink/snp.remove

exit # USING ANOTHER SCRIPT FOR EASIER MODIFICATIONS

# STEP 9: EXTRACT COMMON SNP
for downsample in ${DOWNSAMPLE}; do
        for CHR in ${CHROM}; do

                # Extract the second column (i.e. CHR:POS) from the diploid bim
                awk '{if($1=="'${CHR}'")print $2;}'  ${dir}/glimpse_plink/${DIP}.bim > ${dir}/glimpse_comparison/aaa

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


# STEP 10: CONVERT TO TPED WITH TRANSVERSIONS ONLY
# This analysis only focuses on the transversion SNP sites to minimize the impact of DNA damages in ancient DNA data
for downsample in ${DOWNSAMPLE}; do
        for CHR in ${CHROM}; do
                awk '{if(($5=="G" && $6=="A") || ($5=="A" && $6=="G") || ($5=="C" && $6=="T") || ($5=="T" && $6=="C"))print;}' \
                        ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.bim \
                        | cut -f 2 - > ${dir}/glimpse_comparison/transitions.list
                plink --bfile ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp \
                        --exclude ${dir}/glimpse_comparison/transitions.list \
                        --recode --transpose --out ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.trv
                plink --bfile ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip \
                        --exclude ${dir}/glimpse_comparison/transitions.list \
                        --recode --transpose --out ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.trv
        done
done
rm ${dir}/glimpse_comparison/transitions.list

# STEP 11: COMPARE AND OUTPUT RESULTS
mkdir -p ${dir}/testAccuracy
if [ ! -f ${output} ]; then
        echo "Downsample"       "Chr"   "hetero_rate"   "homo_ref_rate" "homo_alt_rate" "hetero_count"  "homo_ref_count"        "homo_alt_count"        "hetero_total"  "homo_ref_total"        "homo_alt_total" > ${output}
fi # create header

for downsample in ${DOWNSAMPLE}; do
        for CHR in ${CHROM}; do

                # STEP 11.1: COMPARE WITH NON-IMPUTED DATA

                allSNP=`paste ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.trv.tped \
                        ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.trv.tped \
                        | awk -F' ' '{if(($5==$11 && $6==$12) || ($5==$12 && $6==$11))print;}' | wc -l`

                # Heterozygote
                heteroSite=`paste ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.trv.tped \
                        ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.trv.tped \
                        | awk -F' ' '{if($11!=$12 && (($5==$11 && $6==$12) || ($5==$12 && $6==$11))) print;}' | wc -l`
                heteroTotal=`awk -F' ' '{if($5!=$6)print;}' ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.trv.tped | wc -l`


                # Homozygote: retrieve reference and alternative information and get SNPs present in reference only
                cut -d " " -f 2 ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.trv.tped > ${dir}/glimpse_comparison/new.list
                fgrep -wf ${dir}/glimpse_comparison/new.list ${dir}/glimpse_comparison/${REFALT}.chr${CHR}.refalt \
                        > ${dir}/glimpse_comparison/${DIP}.imp.$trv.refalt
                cut -d " " -f 1 ${dir}/glimpse_comparison/${DIP}.imp.trv.refalt > ${dir}/glimpse_comparison/new_list
                fgrep -wf ${dir}/glimpse_comparison/new_list \
                        ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.trv.tped \
                        > ${dir}/glimpse_comparison/diploid_common_with_imputed_and_reference.tped
                fgrep -wf ${dir}/glimpse_comparison/new_list \
                        ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.trv.tped \
                        > ${dir}/glimpse_comparison/imputed_common_with_diploid_and_reference.tped
                rm ${dir}/glimpse_comparison/new.list
                rm ${dir}/glimpse_comparison/new_list
                
                homoRefSite=`paste -d " " ${dir}/glimpse_comparison/${DIP}.imp.trv.refalt \
                        ${dir}/glimpse_comparison/diploid_common_with_imputed_and_reference.tped \
                        ${dir}/glimpse_comparison/imputed_common_with_diploid_and_reference.tped \
                        | awk -F' ' '{if($14==$15 && $2==$14 && (($8==$14 && $9==$15) || ($8==$15 && $9==$14)))print;}' | wc -l`
                homoRefTotal=`paste -d " " ${dir}/glimpse_comparison/${DIP}.imp.trv.refalt \
                        ${dir}/glimpse_comparison/imputed_common_with_diploid_and_reference.tped \
                        | awk -F' ' '{if($8==$9 && $2==$8)print;}' | wc -l`
                homoAltSite=`paste -d " " ${dir}/glimpse_comparison/${DIP}.imp.trv.refalt \
                        ${dir}/glimpse_comparison/diploid_common_with_imputed_and_reference.tped \
                        ${dir}/glimpse_comparison/imputed_common_with_diploid_and_reference.tped \
                        | awk -F' ' '{if($14==$15 && $3==$14 && (($8==$14 && $9==$15) || ($8==$15 && $9==$14)))print;}' | wc -l`
                homoAltTotal=`paste -d " " ${dir}/glimpse_comparison/${DIP}.imp.trv.refalt \
                        ${dir}/glimpse_comparison/imputed_common_with_diploid_and_reference.tped \
                        | awk -F' ' '{if($8==$9 && $3==$8)print;}' | wc -l`

                rm ${dir}/glimpse_comparison/diploid_common_with_imputed_and_reference.tped
                rm ${dir}/glimpse_comparison/imputed_common_with_diploid_and_reference.tped


                # STEP 11.2: CALCULATE CONCORDANCE RATES AND OUTPUT

                heteroRate=`echo  "scale=2; ${heteroSite}  * 100 / ${heteroTotal}"  | bc`
                homoRefRate=`echo "scale=2; ${homoRefSite} * 100 / ${homoRefTotal}" | bc`
                homoAltRate=`echo "scale=2; ${homoAltSite} * 100 / ${homoAltTotal}" | bc`

                echo "${downsample}     ${CHR}  ${heteroRate}   ${homoRefRate}  ${homoAltRate}  ${heteroSite}   ${homoRefSite}  ${homoAltSite}  ${heteroTotal}  ${homoRefTotal} ${homoAltTotal}" >> ${output}

        done
done
rm ${dir}/glimpse_comparison/${DIP}.imp.trv.refalt

# STEP 12: REMOVE USELESS FILES
# I seperated each file so that it is easy to tweak
for downsample in ${DOWNSAMPLE}; do
        for CHR in ${CHROM}; do
                rm ${dir}/glimpse_vcf/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.vcf*
                rm ${dir}/glimpse_imputed/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED*

                rm ${dir}/glimpse_ligated/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.bcf*
#               rm ${dir}/glimpse_ligated/${sample}.${name}.${downsample}.chr${CHR}.likehood.REGROUPED.merged.vcf
                rm ${dir}/glimpse_ligated/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.vcf

#               rm ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.bed
#               rm ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.bim
#               rm ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.fam
                rm ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.log
                rm ${dir}/glimpse_plink/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.nosex

#               rm ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.bed
#               rm ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.bim
#               rm ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.fam
#               rm ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.log
#               rm ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.nosex
#               
#               rm ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.bed
#               rm ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.bim
#               rm ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.fam
#               rm ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.log
#               rm ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.nosex
#               
#               rm ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.trv.tped
#               rm ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.trv.tfam
#               rm ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.trv.log
#               rm ${dir}/glimpse_comparison/${DIP}_${downsample}_chr${CHR}.imp.trv.nosex
#               
#               rm ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.trv.tped
#               rm ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.trv.tfam
#               rm ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.trv.log
#               rm ${dir}/glimpse_comparison/${sample}.${name}.${downsample}.chr${CHR}.likelihood.REGROUPED.merged.GP${threshold}.nodup.dip.trv.nosex
        done
done
#rm ${dir}/glimpse_comparison/${REFALT}.chr*.refalt

echo "END: Your results have been outputed to ${output}."
