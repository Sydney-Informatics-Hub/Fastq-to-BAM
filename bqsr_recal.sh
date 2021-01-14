#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: run GATK Base Recalibrator over parallel tasks
# Usage: this script is executed by bqsr_recal_run_parallel.pbs
# Author: Cali Willet
# cali.willet@sydney.edu.au
# Date last modified: 24/07/2020
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

module load gatk/4.1.2.0 #loads also java/jdk1.8.0_60

labSampleID=`echo $1 | cut -d ',' -f 1` 
intNum=`echo $1 | cut -d ',' -f 2`
interval=`echo $1 | cut -d ',' -f 3`

bam=./Dedup_sort/${labSampleID}.coordSorted.dedup.bam 
log=./Logs/BQSR_recal/${labSampleID}.${intNum}.log
err=./Error_capture/BQSR_recal/${labSampleID}.${intNum}.err
out=./BQSR_tables/${labSampleID}.${intNum}.recal_data.table

ref=./Reference/hs38DH.fasta
thousandG_indels=./Reference/Known_vars/Homo_sapiens_assembly38.known_indels.vcf
gold_standard_indels=./Reference/Known_vars/Mills_and_1000G_gold_standard.indels.hg38.vcf
dbsnp=./Reference/Known_vars/Homo_sapiens_assembly38.dbsnp138.vcf

gatk BaseRecalibrator \
	--java-options "-Xmx6G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	-R $ref \
	-L $interval \
	-I $bam  \
	--known-sites $dbsnp \
	--known-sites $thousandG_indels \
	--known-sites $gold_standard_indels \
	-O $out >> $log 2>&1


if grep -q -i error $log
then 
        printf "Error in GATK log ${log}\n" >> $err
fi 

if grep -q Exception $log
then 
        printf "Exception in GATK log ${log}\n" >> $err
fi
	
	
