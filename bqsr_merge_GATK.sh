#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: merge recalibrated split BAM files into a final 
# BAM per sample with GATK GatherBamFiles, parallel by sample
# Usage: this script is executed by bqsr_merge_GATK_run_parallel.pbs
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

set -e 


labSampleID=$1

jvm=28 # Set this variable to 28 for hugemem queue or 10 for normal queue


module load gatk/4.1.2.0
module load samtools/1.10

log=./Logs/BQSR_merge/${labSampleID}.log
err=./Error_capture/BQSR_merge/${labSampleID}.err
bam_out=./Final_bams/${labSampleID}.final.bam
ref=./Reference/hs38DH.fasta


\rm -rf $bam_out #attempt to set stripe on existing file will cause fatal error
lfs setstripe -c 15 $bam_out


gatk GatherBamFiles \
	--java-options "-Xmx${jvm}G" \
	-I BQSR_merge_lists/${labSampleID}.list \
	-O $bam_out \
	-R $ref \
	--CREATE_INDEX=true >> $log 2>&1
		
	#\
	#--CREATE_MD5_FILE=true #Do not create md5 in this job, as it will cost 3x SU


if ! samtools quickcheck $bam_out
then 
        printf "Corrupted or missing BAM\n" > $err  
fi

if grep -q -i error $log
then 
        printf "Error in GATK log ${log}\n" >> $err
fi 

if grep -q -i exception $log
then 
        printf "Exception in GATK log ${log}\n" >> $err
fi
