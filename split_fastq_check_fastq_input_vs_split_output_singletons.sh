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
fastq=$(basename $fastq_name | sed 's/.fastq.gz//')
out=./Check_fastq_split/${fastq}.check
rm -rf $out

err=0
#1 - Get original number of reads from fastQC data.txt
qc1=$(ls ./FastQC/${fastq}*fastqc/fastqc_data.txt)
fastqc_1=$(grep "Total Sequences" $qc1 | awk '{print $3}')

#2 - Get number of reads processed from fastp logs
fastp_log=./Logs/Fastp/${fastq}.log

fastp_1=$(grep -A 1 "Read1 before filtering" $fastp_log | tail -1 | awk '{print $3}')
fastp_tot=$(grep -m 1 "reads passed filter" $fastp_log | awk '{print $4}')

if [[ $fastp_1 -ne $fastp_tot ]]
then
        printf "$fastq error: fastp reads before and after filtering (no filtering applied) do not match - $fastp_1 and $fastp_tot respectively\n" >> $out
        ((err++))
fi

if [[ $fastp_1 -ne $fastqc_1 ]]
then
        printf "$fastq error: fastp and fastQC read counts do not match - $fastp_1 and $fastqc_1 respectively\n" >> $out
        ((err++))
fi

#3 - Get number of reads output from line-counting the split fastq
splits=./Fastq_split

split_1_lines=$(ls ${splits}/*${fastq}*f*q.gz | parallel -j $NCPUS --will-cite "zcat {} | wc -l " | awk '{s+=$1} END {print s}') # check regex
split_pairs=$(expr $split_1_lines \/ 4)

#Check 3 sources

if [[ $fastp_1 -ne $split_pairs ]]
then
        printf "$fastq error: fastp input and fastp split output read counts do not match - $fastp_1 and $split_pairs respectively\n" >> $out
        ((err++))
fi

if [[ $err -eq 0 ]]
then
        printf "$fastq has passed all checks\n" >> $out
fi


