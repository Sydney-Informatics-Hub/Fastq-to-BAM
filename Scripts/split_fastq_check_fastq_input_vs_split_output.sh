#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: this script is run by split_fastq_check_fastq_input_vs_split_output_run_parallel.pbs
# Version: 2.0
#
# For more details see: https://github.com/Sydney-Informatics-Hub/Fastq-to-BAM
#
# If you use this script towards a publication, please acknowledge the
# Sydney Informatics Hub (or co-authorship, where appropriate).
#
# Suggested acknowledgement:
# The authors acknowledge the support provided by the Sydney Informatics Hub,
# a Core Research Facility of the University of Sydney. This research/project
# was undertaken with the assistance of resources and services from the National
# Computational Infrastructure (NCI), which is supported by the Australian
# Government, and the Australian BioCommons which is enabled by NCRIS via
# Bioplatforms Australia funding.
#
##########################################################################

#1 - Get original number of pairs from fastQC data.txt
#2 - Get number of pairs processed from fastp logs
#3 - Get number of pairs output from line-counting the split fastq
#4 - print warning if mismatch detected at any step; print "sampe OK" if no errors

fastq_name=$1
fastq=$(basename $fastq_name)
out=../Check_split_fastq/${fastq}.check
rm -rf $out

err=0
#1 - Get original number of pairs from fastQC data.txt
qc1=$(ls ../FastQC_batch2/${fastq}*R1*fastqc/fastqc_data.txt)
qc2=$(ls ../FastQC_batch2/${fastq}*R2*fastqc/fastqc_data.txt)

fastqc_1=$(grep "Total Sequences" $qc1 | awk '{print $3}')
fastqc_2=$(grep "Total Sequences" $qc2 | awk '{print $3}')

# Check number of reads in R1 and R2 match
if [[ $fastqc_1 -ne $fastqc_2 ]]
then
        printf "$fastq error: fastqc R1 and R2 read counts do not match - $fastqc_1 and $fastqc_2 respectively\n" > $out
        ((err++))
fi

#2 - Get number of pairs processed from fastp logs
fastp_log=./Logs/split_fastq/${fastq}.log

fastp_1=$(grep -A 1 "Read1 before filtering" $fastp_log | tail -1 | awk '{print $3}')
fastp_2=$(grep -A 1 "Read2 before filtering" $fastp_log | tail -1 | awk '{print $3}')
fastp_tot=$(grep -m 1 "reads passed filter" $fastp_log | awk '{print $4}')
fastp_sum=$(expr $fastp_1 + $fastp_2)

if [[ $fastp_sum -ne $fastp_tot ]]
then
        printf "$fastq error: fastp R1 and R2 input sum does not match fastp total reads read - $fastp_sum and $fastp_tot respectively\n" >> $out
        ((err++))
fi

if [[ $fastp_1 -ne $fastp_2 ]]
then
        printf "$fastq error: fastp R1 and R2 read counts do not match - $fastp_1 and $fastp_2 respectively\n" >> $out
        ((err++))
fi

#3 - Get number of pairs output from line-counting the split fastq
splits=../Split_fastq

split_1_lines=$(ls ${splits}/*${fastq}_R1*f*q.gz | parallel -j $NCPUS --will-cite "zcat {} | wc -l " | awk '{s+=$1} END {print s}') # check regex
split_2_lines=$(ls ${splits}/*${fastq}_R2*f*q.gz | parallel -j $NCPUS --will-cite "zcat {} | wc -l " | awk '{s+=$1} END {print s}')

if [[ $split_1_lines -ne $split_2_lines ]]
then
        printf "$fastq error: fastp split R1 and R2 line counts do not match - $split_1_lines and $split_1_lines respectively\n" >> $out
        ((err++))
fi

split_pairs=$(expr $split_1_lines \/ 4)

#Check 3 sources
if [[ $fastqc_1 -ne $fastp_1 ]]
then
        printf "$fastq error: fastQC and fastp pair counts do not match - $fastqc_1 and $fastp_1 respectively\n" >> $out
        ((err++))
fi

if [[ $fastp_1 -ne $split_pairs ]]
then
        printf "$fastq error: fastp input and fastp split output pair counts do not match - $fastp_1 and $split_pairs respectively\n" >> $out
        ((err++))
fi

if [[ $err -eq 0 ]]
then
        printf "$fastq has passed all checks\n" >> $out
fi

