#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: run GATK ApplyQSR over parallel tasks
# Usage: this script is executed by bqsr_apply_run_parallel.pbs
# Details:
# 	Compression needs to be applied at the ApplyBQSR step if merging with 
# 	GATK. SAMbamba merge makes compression=5 BAMs, but GATK merge can not
# 	compress (despite flags to that effect)
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
module load samtools/1.10

# special chracters in contig names requires this specific method to read in positional arguments:
IFS=','
params=($1)
labSampleID=${params[0]}
contig=${params[1]}
counter=${params[2]}


ref=./Reference/hs38DH.fasta
table=./BQSR_tables/${labSampleID}.recal_data.table
bam_in=./Dedup_sort/${labSampleID}.coordSorted.dedup.bam
bam_out=./BQSR_apply/${labSampleID}.${counter}.recal.bam
log=./Logs/BQSR_apply/${labSampleID}.${counter}.log
err=./Error_capture/BQSR_apply/${labSampleID}.${counter}.err

gatk ApplyBQSR \
	--java-options "-Xmx6G -Dsamjdk.compression_level=5" \
	-R $ref \
	-I $bam_in  \
	--bqsr-recal-file $table \
	-L $contig \
	-O $bam_out >> $log 2>&1


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
