#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: check that the merged alignments meet expected size
# Usage: bash merge_align_check_fastq_input_vs_bam_output_sizes.sh <cohort_name>
# Details:
#	 Expected ratio is ~ 1.4 - 1.5 (merged BAM larger than fastq)
# 	This can deviate from expected if a preprocessing step has output fastq
# 	in a different compression level
# 	Supply cohort name as positional argument to bash if desired
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

output=./Logs/Align_merge/merge_align_check_fastq_input_vs_bam_output_sizes.txt

if [ "$cohort" ]
then
	output=./Logs/Align_merge/${cohort}_merge_align_check_fastq_input_vs_bam_output_sizes.txt	
fi

rm -rf $output

printf "#ID\tFastq\tAlignMerged\tRatio\n" >> $output

awk 'NR>1' samples.config | while read LINE
do 
        sample=`echo $LINE | cut -d ' ' -f 1`
        labSampleID=`echo $LINE | cut -d ' ' -f 2`
	fq_size=$(find  ./Fastq -type f -name *${sample}* -print | xargs ls -lh | awk '{print $5}' | sed 's/G$//'| awk '{s+=$1} END {print s}')
	bam_size=$(ls -lh Align_merged/${labSampleID}.merged.nameSorted.bam | awk '{print $5}' | sed 's/G$//')
	ratio=$(bc -l <<< "scale=2; ${bam_size}/${fq_size}")
	printf "${labSampleID}\t${fq_size}\t${bam_size}\t${ratio}\n" >> $output
done


