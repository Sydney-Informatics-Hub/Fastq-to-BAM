#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: make inputs file for parallel exectuion of merging
# split BAM files per sample
# Usage: bash merge_align_make_input.sh <cohort_name>
# Details:
# 	Create input list for merging parallel alignment output
# 	Read sample info from <cohort>.config
# 	Assumes all non-cancer samples are designated 'N' (normal)
# 	or 'B' (blood), and are lower coverage. All other phenotype 
# 	IDs are assigned to 'tumour'. Splitting in this way is 
# 	recommended when there is a binomial difference in coverage, 
# 	eg 30X/60X cancer projects. For cohorts where the input data 
# 	varies by a large degree and is not really binomial, the 
# 	grouping can be done on file size (see 'merge_align_make_input_by_fastq_size.sh')
# 	If no binomial grouping is desired, change group=true to group=false
# 	Provide cohort name as argument
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

group=true

t_input=./Inputs/merge_align.inputs-tumour
n_input=./Inputs/merge_align.inputs-normal
input=./Inputs/merge_align.inputs

rm -f $t_input
rm -f $n_input
rm -f $input

awk 'NR>1' ${cohort}.config | while read LINE
do 
	sample=`echo $LINE | cut -d ' ' -f 1`
	labSampleID=`echo $LINE | cut -d ' ' -f 2`

	if [[ $group = true ]]
	then
        	if [[ $labSampleID = *-N|B ]]
		then
			printf "${sample},${labSampleID}\n" >> $n_input
		else
			printf "${sample},${labSampleID}\n" >> $t_input
		fi
	else
		printf "${sample},${labSampleID}\n" >> $input		
	fi					
done	

if [ -f $t_input ]
then
	tasks=`wc -l < $t_input`
	printf "Number of tumour sample merge tasks to run: ${tasks}\n"
fi

if [ -f $n_input ]
then
	tasks=`wc -l < $n_input`
	printf "Number of normal sample merge tasks to run: ${tasks}\n"
fi

if [ -f $input ]
then
	tasks=`wc -l < $input`
	printf "Number of merge tasks to run: ${tasks}\n"
fi
