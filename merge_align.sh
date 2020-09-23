#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: Merge the split BAMs for each sample into 
# sample-level, queryname sorted BAMs
# Usage: this script is executed by merge_align_run_parallel.pbs
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

module load samtools/1.10
module load sambamba/0.7.1

sample=`echo $1 | cut -d ',' -f 1` 
labSampleID=`echo $1 | cut -d ',' -f 2` 


cd ./Align_split


err=../Error_capture/Align_merge/${labSampleID}.err
rm -rf $err

bams=$(find . -name "*${sample}*.bam" | xargs echo) #Check regex per batch - for CSCC, this is OK, all sample names are unique with this simple regex

rm -rf ../Align_merged/${labSampleID}.merged.nameSorted.bam   #attempt to set stripe on existing file will cause fatal error
lfs setstripe -c 15 ../Align_merged/${labSampleID}.merged.nameSorted.bam


sambamba merge -t $NCPUS ../Align_merged/${labSampleID}.merged.nameSorted.bam $bams


if ! samtools quickcheck ../Align_merged/${labSampleID}.merged.nameSorted.bam
then 
	printf "Corrupted or missing BAM\n" > $err  
fi
