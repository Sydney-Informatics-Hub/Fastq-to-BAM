#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: create BAM lists for merging recalibrated split 
# BAM files with GATK GatherBamFiles
# Usage: bash bqsr_merge_make_bamLists.sh <cohort_name>
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

cohort=$1

mkdir -p ./Inputs/BQSR_merge_lists

samples=$(awk 'NR > 1 {print $2}' samples.config)
samples=($samples)

for (( s = 0; s < ${#samples[@]}; s++ ))
do
	labSampleID=${samples[$s]}
	list=./Inputs/BQSR_merge_lists/${labSampleID}.list
	\rm -rf $list
	for (( i = 0; i < 3366; i++ ))
	do 
		printf "BQSR_apply/${labSampleID}.${i}.recal.bam\n" >> $list		
	done
	printf "BQSR_apply/${labSampleID}.unmapped.recal.bam\n" >> $list
done
