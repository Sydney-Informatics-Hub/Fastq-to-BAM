#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: Create input list for bam metrics
# Usage: bash bam_metrics_make_input.sh <cohort_name>
# Details:
# 	Read sample info from <cohort>.config. Provide cohort name as argument.
# 	Assumes all non-cancer samples are designated 'N' (normal)
# 	or 'B' (blood), and are lower coverage. All other phenotype 
# 	IDs are assigned to 'tumour'. Splitting in this way is 
# 	recommended when there is a binomial difference in coverage, 
# 	eg 30X/60X cancer projects. If no binomial grouping is desired, 
#	change group=true to group=false
#
# Author: Cali Willet
# cali.willet@sydney.edu.au
# Date last modified: 12/08/2020
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

t_input=./Inputs/bam_metrics.inputs-tumour
n_input=./Inputs/bam_metrics.inputs-normal
input=./Inputs/bam_metrics.inputs

rm -f $t_input
rm -f $n_input
rm -f $input

awk 'NR>1' ${cohort}.config | while read LINE
do         
	labSampleID=`echo $LINE | cut -d ' ' -f 2`	
	
	if [[ $group = true ]]
	then
		if [[ $labSampleID = *-N ]]
       		then
                	printf "${labSampleID}\n" >> $n_input
        	else
                	printf "${labSampleID}\n" >> $t_input
        	fi
	else
		printf "${labSampleID}\n" >> $input
	fi					
done	


if [ -f $t_input ]
then
	tasks=`wc -l < $t_input`
	printf "Number of tumour sample bam metrics tasks to run: ${tasks}\n"
fi

if [ -f $n_input ]
then
	tasks=`wc -l < $n_input`
	printf "Number of normal sample bam metrics tasks to run: ${tasks}\n"
fi


if [ -f $input ]
then
	tasks=`wc -l < $input`
	printf "Number of bam metrics tasks to run: ${tasks}\n"
fi
