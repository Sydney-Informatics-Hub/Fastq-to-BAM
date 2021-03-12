#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: check unique ID matching to fastq
# Usage: bash check_unique_fastq_ids.sh <config>
# Details:
#       Run this before submitting the alignment job. 
# 	FastQC and split fastq can be run before this.
# 	Check that each sample in column 1 of config has 
#	a unique name, ie that the sample ID can be used 
#	to find the fastq for the sample without incorrectly 
#	matching fastq belonging to other samples
#
# Author: Cali Willet
# cali.willet@sydney.edu.au
# Date last modified: 12/3/21
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


if [ -z "$1" ]
then
        echo "Please run this script with the name of your config e.g. bash check_unique_fastq_ids.sh <config>"
        exit
fi

config=$1

samples=$(awk 'NR>1' $config | cut -f 1)
samples=($samples)

f=0

for ((i = 0; i < ${#samples[@]}; i++ ))
do 
	sample=${samples[$i]}
	m=$(awk -v sample="$sample" '$1 ~ sample' $config | wc -l)
	if [[ $m -gt 1 ]]
	then
		printf "\nFAIL: $sample has matched $m IDs in column 1\nIt is imperative to rectify this before running align!\n"
		printf "Matches are:\n"
		awk -v sample="$sample" '$1 ~ sample' $config
		((f++))
	fi
done

if [[ $f == 0 ]]
then
	printf "\n################################################################\n"
	printf "\nPASS: All sample IDs in column 1 of $config are unique.\nPlease run check_unique_fastq_ids.sh to verify\n"
	printf "\n################################################################\n\n"
else 	
	printf "\n################################################################\n"
	printf "\nFAIL: there are $f samples with non-unique IDs in column 1 of ${config}\n"
	printf "\n################################################################\n"
	printf "\nIt is IMPERATIVE that this is rectified before running alignment!\n"
	printf "\nTo fix this, edit the nonunique IDs in column 1 to a longer match to\nthe fastq. If this is not possible, rename the fastq.\n"
	printf "\n################################################################\n"
	printf "\nAfter correcting the nonunique IDs, run this script until no fails\nare found, then run check_unique_fastq_ids.sh to verify\n\n"
	printf "\n################################################################\n\n"
fi	  
