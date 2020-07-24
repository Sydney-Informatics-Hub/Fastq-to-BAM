#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: check the output BAM size is as expected after 
# sorting and duplicate read marking
# Usage: bash dedup_sort_check_sorted_bam_sizes.sh <cohort_name>
# Details:
#	Expected ratio of merged unsorted BAM to sorted BAM is ~ 0.55 
#	(sorted BAM smaller than unsorted). Supply cohort name as 
#	first positional argument to bash if desired. Output will be
# 	in ./Logs/Dedup_sort
#
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

output=./Logs/Dedup_sort/check_unsorted_vs_sorted_bam_output_sizes.txt

if [ "$cohort" ]
then
	output=./Logs/Dedup_sort/${cohort}_check_unsorted_vs_sorted_bam_output_sizes.txt
fi

rm -rf $output

printf "#ID\tAlignMerged\tDedupSort\n" >> $output

awk 'NR>1' samples.config | while read LINE
do 
	labSampleID=`echo $LINE | cut -d ' ' -f 2`
	merged_bam_size=$(ls -lh Align_merged/${labSampleID}.merged.nameSorted.bam | awk '{print $5}' | sed 's/G$//')
        sorted_bam_size=$(ls -lh Dedup_sort/${labSampleID}.coordSorted.dedup.bam | awk '{print $5}' | sed 's/G$//')	
	ratio=$(bc -l <<< "scale=2; ${sorted_bam_size}/${merged_bam_size}")
        printf "${labSampleID}\t${merged_bam_size}\t${sorted_bam_size}\t${ratio}\n" >> $output
done
