#!/bin/bash

###############################################################################################
# UPDATED PROCESSING PIPELINE FOR HS37D5                                                      #
# BRING MORE TO LINE WITH THE UPDATED SEQPRO FORMAT                                           #
                                                                                              #
# This script will run from its current folder, usually "raw_data"                            #
# and will create additional folders on its level                                             #
# such as "trimmed" and "alignements"                                                         #
                                                                                              #
# PARAMETERS #                                                                                #
threads=25   # number of threads to be used                                                   #
quaBase=25   # base quality for trimming (will only be performed if trimming is required)     #
minLen=25    # minimum length for trimming (will only be performed if trimming is required)   #
mapQua=20    # mapping quality (can be set to 0 if you want to filter later with HapSNP)      #
###############################################################################################


####################                   DOWNLOAD RAW DATA                   ####################
echo "DOWNLOADING RAW DATA..." # download fastq from ENA
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR286/007/ERR2868047/ERR2868047.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR286/008/ERR2868048/ERR2868048.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR286/000/ERR2868050/ERR2868050.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR286/001/ERR2868051/ERR2868051.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR286/009/ERR2868035/ERR2868035.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR286/006/ERR2868046/ERR2868046.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR286/005/ERR2868045/ERR2868045.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR286/008/ERR2868038/ERR2868038.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR286/009/ERR2868039/ERR2868039.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR286/006/ERR2868036/ERR2868036.fastq.gz
#mv ERR2868047.fastq.gz TrailCreek.fastq.gz
#mv ERR2868048.fastq.gz 19651.fastq.gz
#mv ERR2868050.fastq.gz AHUR_2064.fastq.gz
#mv ERR2868051.fastq.gz Lovelock1.fastq.gz
#mv ERR2868035.fastq.gz Lovelock2.fastq.gz
#mv ERR2868046.fastq.gz Lovelock3.fastq.gz
#mv ERR2868045.fastq.gz Lovelock4.fastq.gz
#mv ERR2868038.fastq.gz SA5832.fastq.gz
#mv ERR2868039.fastq.gz A460.fastq.gz
#mv ERR2868036.fastq.gz Aconcagua.fastq.gz
###############################################################################################


####################                   SET UP FILE NAMES                   ####################
list=`ls -lS *fastq.gz | awk '{print $9}'` # order list by file size so that next step starts only when biggest file has finished
sample=()
for file in $list; do
        sample+=$(echo $file | rev | cut -f3-20 -d. | rev); # get sample name from file
        sample+=' '                                         # add space to seperate elements in list
done
array=($sample)                                             # convert list to array (for iteration later-on)
###############################################################################################

####################                        TRIMMING                       ####################
echo "TRIMMING ADAPTERS..."
mkdir -p ../trimmed

# trim all files but one in the background
for ((i=1; i<${#array[@]}; i++)); do
        gunzip -c ${array[i]}.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths) {print l, lengths[l]}; print "total reads: " counter}' \
                > ${array[i]}.read_lengths_distribution.txt & # count length of each read
done
gunzip -c ${array[0]}.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths) {print l, lengths[l]}; print "total reads: " counter}' \
        > ${array[0]}.read_lengths_distribution.txt & # count length of each read

for ((i=1; i<${#array[@]}; i++)); do
        # if all reads are of the same length, then data was not trimmed
        if [ $(wc -l ${array[i]}.read_lengths_distribution.txt | cut -d' ' -f1) -eq 1 ]; then
                /home/niallc/anaconda3/bin/AdapterRemoval \
                        --qualitybase $quaBase  \
                        --file1 ${array[i]}.fastq.gz \
                        --basename ../trimmed/${array[i]}.trimmed${minLeng}bp${quaBase}q \
                        --trimns \
                        --trimqualities \
                        --minquality 2 \
                        --gzip \
                        --minlength $minLen \
                        > ../trimmed/${array[i]}.trimmed${minLeng}bp${quaBase}q.AdapterRemoval.errors &
        else
                scp ${array[i]}.fastq.gz ../trimmed/${array[i]}.trimmed${minLeng}bp${quaBase}q.truncated.gz
        fi
done
# trim last file in foreground
if [ $(wc -l ${array[0]}.read_lengths_distribution.txt | cut -d' ' -f1) -eq 1 ]; then
        /home/niallc/anaconda3/bin/AdapterRemoval \
                --qualitybase $quaBase  \
                --file1 ${array[0]}.fastq.gz \
                --basename ../trimmed/${array[i]}.trimmed${minLeng}bp${quaBase}q \
                --trimns \
                --trimqualities \
                --minquality 2 \
                --gzip \
                --minlength $minLen \
                > ../trimmed/${array[0]}.trimmed${minLeng}bp${quaBase}q.AdapterRemoval.errors &
else
        scp ${array[0]}.fastq.gz ../trimmed/${array[0]}.trimmed${minLeng}bp${quaBase}q.truncated.gz
fi

for ((i=0; i<${#array[@]}; i++)); do
        if [ $(wc -l ${array[i]}.read_lengths_distribution.txt | cut -d' ' -f1) -eq 1 ]; then
                echo "Trimming required for ${array[i]}"
        else
                echo "Trimming NOT required for ${array[i]}"
        fi
done

# remove temporary file
for ((i=0; i<${#array[@]}; i++)); do rm ${array[i]}.read_lengths_distribution.txt; done
###############################################################################################


####################                       ALIGNMENTS                      ####################
echo "ALIGNING: CREATING SAI FILES..."
mkdir -p ../alignments
mkdir -p ../alignments/hs37d5
mkdir -p ../alignments/hs37d5/sai_files
mkdir -p ../alignments/hs37d5/bam_files
mkdir -p ../alignments/hs37d5/final_files


####################                       SAI TO BAM                      ####################
for ((i=0; i<${#array[@]}; i++)); do
        bwa aln -l 16500 -n 0.01 -o 2 -t $threads /Reference_Genomes/hs37d5/hs37d5.fa \
                ../trimmed/${array[i]}.trimmed${minLeng}bp${quaBase}q.truncated.gz \
                > ../alignments/hs37d5/sai_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sai \
                2> ../alignments/hs37d5/sai_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.bwa.log
done
###############################################################################################


####################                PROCESS BAM (pre-merge)                ####################
# Make BAM from SAI
for ((i=1; i<${#array[@]}; i++)); do
        bwa samse /Reference_Genomes/hs37d5/hs37d5.fa ../alignments/hs37d5/sai_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sai \
                ${array[i]}.fastq.gz | samtools view -Sb  -F4 \
                > ../alignments/hs37d5/bam_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.bam &
done
bwa samse /Reference_Genomes/hs37d5/hs37d5.fa ../alignments/hs37d5/sai_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sai \
        ${array[0]}.fastq.gz | samtools view -Sb  -F4 \
        > ../alignments/hs37d5/bam_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.bam

# Sort BAM file
for ((i=1; i<${#array[@]}; i++)); do
        samtools sort ../alignments/hs37d5/bam_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.bam \
                -o ../alignments/hs37d5/bam_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.bam &
done
samtools sort ../alignments/hs37d5/bam_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.bam \
        -o ../alignments/hs37d5/bam_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.bam

# Add read groups: something to identify the original downloaded file if there are multiple for the one individual
for ((i=1; i<${#array[@]}; i++)); do
        group=$(echo ${array[i]}.fastq.gz | cut -f1 -d-); lib=$(echo ${array[i]}.fastq.gz | cut -f4-6 -d-); seq=$(echo ${array[i]}.fastq.gz | cut -f1 -d. | cut -f7-8 -d-);
        java -Xmx5g -jar /Software/AddOrReplaceReadGroups.jar \
                I=../alignments/hs37d5/bam_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.bam \
                O=../alignments/hs37d5/bam_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.bam \
                RGLB=$lib RGSM=$group RGPU=Illumina RGPL=ILLUMINA RGID=$seq RGDS=${array[i]}.fastq.gz &
done
group=$(echo ${array[0]}.fastq.gz | cut -f1 -d-); lib=$(echo ${array[0]}.fastq.gz | cut -f4-6 -d-); seq=$(echo ${array[0]}.fastq.gz | cut -f1 -d. | cut -f7-8 -d-);
java -Xmx5g -jar /Software/AddOrReplaceReadGroups.jar \
        I=../alignments/hs37d5/bam_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.bam \
        O=../alignments/hs37d5/bam_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.bam \
        RGLB=$lib RGSM=$group RGPU=Illumina RGPL=ILLUMINA RGID=$seq RGDS=${array[0]}.fastq.gz

# Duplicate removal
for ((i=1; i<${#array[@]}; i++)); do
        java -jar /Software/picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE \
                I=../alignments/hs37d5/bam_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.bam \
                O=../alignments/hs37d5/bam_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.bam \
                M=marked_dup_metrics.txt &
done
java -jar /Software/picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE \
        I=../alignments/hs37d5/bam_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.bam \
        O=../alignments/hs37d5/bam_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.bam \
        M=marked_dup_metrics.txt
rm ../alignments/hs37d5/bam_files/marked_dup_metrics.txt
# Note: there are further optimal options found in https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
###############################################################################################


####################                PROCESS BAM (post-merge)               ####################
for ((i=1; i<${#array[@]}; i++)); do
        samtools index ../alignments/hs37d5/bam_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.bam &
done
samtools index ../alignments/hs37d5/bam_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.bam

for ((i=1; i<${#array[@]}; i++)); do
        java -jar /Software/GenomeAnalysisTK.jar \
                -T RealignerTargetCreator \
                -R /Reference_Genomes/hs37d5/hs37d5.fa \
                -I ../alignments/hs37d5/bam_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.bam \
                -o ../alignments/hs37d5/final_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.intervals &
done
java -jar /Software/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R /Reference_Genomes/hs37d5/hs37d5.fa \
        -I ../alignments/hs37d5/bam_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.bam \
        -o ../alignments/hs37d5/final_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.intervals

for ((i=1; i<${#array[@]}; i++)); do
        java -jar /Software/GenomeAnalysisTK.jar \
                -T IndelRealigner \
                -R /Reference_Genomes/hs37d5/hs37d5.fa \
                -I ../alignments/hs37d5/bam_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.bam \
                -o ../alignments/hs37d5/final_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.bam \
                -targetIntervals ../alignments/hs37d5/final_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.intervals &
done
java -jar /Software/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R /Reference_Genomes/hs37d5/hs37d5.fa \
        -I ../alignments/hs37d5/bam_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.bam \
        -o ../alignments/hs37d5/final_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.bam \
        -targetIntervals ../alignments/hs37d5/final_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.intervals

# Mapping: not needed if we want to use mapping 0 and filter during HapSNP
for ((i=1; i<${#array[@]}; i++)); do
        samtools view -b -q${mapQua} ../alignments/hs37d5/final_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.bam \
                > ../alignments/hs37d5/final_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.mq${mapQua}.bam &
done
samtools view -b -q${mapQua} ../alignments/hs37d5/final_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.bam \
        > ../alignments/hs37d5/final_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.mq${mapQua}.bam
       
# Softclipping
for ((i=1; i<${#array[@]}; i++)); do
        samtools view -h ../alignments/hs37d5/final_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.mq${mapQua}.bam \
                | python /Programs/softclip.py \
                | samtools view -Sb \
                > ../alignments/hs37d5/final_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.mq${mapQua}.softclip.bam
done
samtools view -h ../alignments/hs37d5/final_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.mq${mapQua}.bam \
        | python /Programs/softclip.py \
        | samtools view -Sb \
        > ../alignments/hs37d5/final_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.mq${mapQua}.softclip.bam

for ((i=1; i<${#array[@]}; i++)); do
        samtools index ../alignments/hs37d5/final_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.mq${mapQua}.softclip.bam
done
samtools index ../alignments/hs37d5/final_files/${array[0]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.mq${mapQua}.softclip.bam

# Check the BAM: get coverage stats and run read length distribution
for ((i=0; i<${#array[@]}; i++)); do

        /Software/qualimap_v2.1.1/qualimap bamqc --java-mem-size=50G \
                -bam ../alignments/hs37d5/final_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.mq${mapQua}.softclip.bam &

        samtools view ../alignments/hs37d5/final_files/${array[i]}.trimmed${minLeng}bp${quaBase}q.sorted.grouped.duprm.indel.mq${mapQua}.softclip.bam \
                | awk '{print length($10);}' | sort | uniq -c > ../alignments/hs37d5/final_files/${array[i]}_read_length_distribution.txt &
done
