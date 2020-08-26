#!/bin/bash 

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: Split each fastq pair into 2,000,000 lines 
# (500,000 reads) with fastp, executed in parallel
# Usage: this script is executed by split_fastq_run_parallel.pbs
# Details:
# 	Check that the regex at 'ls' matches your data. Edit as required.
# Author: Cali Willet
# cali.willet@sydney.edu.au
# Date last modified: 26/08/2020
#
# If you use this script towards a publication, please acknowledge the
# Sydney Informatics Hub (or co-authorship, where appropriate).
#
# Suggested acknowledgement:
# The authors acknowledge the scientific and technical assistance 
# <or e.g. bioinformatics assistance of <PERSON>> of Sydney Informatics
# Hub and resources and services from the National Computational 
# Infrastructure (NCI), which is supported by the Australian Government
# with access facilitated by the University of Sydney.
# 
#########################################################

fastp=/scratch/<project>/apps/fastp/0.20.0/fastp #or replace with module load if appropriate

fqpair=`echo $1 | cut -d ',' -f 1`
file=$(basename $fqpair)

fq1=$(ls ${fqpair}*R1.f*q.gz) #Must check regex for each batch.
fq2=$(ls ${fqpair}*R2.f*q.gz) 

log=./Logs/Fastp/${file}.log

$fastp -i ${fq1} \
	-I ${fq2} \
	-AGQL \
	-w $NCPUS \
	-S 2000000 \
	-d 0 \
	--out1 Fastq_split/${file}_R1.fastq.gz \
	--out2 Fastq_split/${file}_R2.fastq.gz 2>${log}
