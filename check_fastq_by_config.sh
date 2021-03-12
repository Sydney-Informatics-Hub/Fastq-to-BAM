#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: check unique ID matching to fastq part 2
# Usage: bash check_fastq_by_config.sh <config>
# Details:
#       Run this before submitting the alignment job. 
# 	FastQC and split fastq can be run before this.
# 	First run check_unique_fastq_ids.sh. After correcting
#	any nonunique IDs in column 1 of config, run this
#	to verify that no nonunique fastq matches are occurring
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
        echo "Please run this script with the name of your config e.g. bash check_fastq_by_config.sh <config>"
        exit
fi

config=$1

samples=$(awk 'NR>1' $config | cut -f 1)
samples=($samples)

c=0
s=0

fq=$(ls -1 ./Fastq/*f*q.gz | wc -l) # check this suffix matches your dataset, edit as required

echo $fq

for ((i = 0; i < ${#samples[@]}; i++ ))
do 
	sample=${samples[$i]}
	((s++))
	found=$(ls -1 ./Fastq/*${sample}* | wc -l)
	((c+=$found))
	
	if [[ $found == 0 ]]
	then 
		printf "FAIL: No fastq found for ${sample}\n"
	else
		printf "\t$found fastq files found for sample $sample\n"
	fi
done

printf "\n################################################################\n"

if [[ $fq -ne $c ]]
then
	printf "\nFAIL: the number of fastq in ./Fastq (${fq}) does not match the\nnumber of fastq matched to sample IDs in ${config} (${c})\n"
	printf "\n################################################################\n"
	printf "\nIt is IMPERATIVE that this is rectified before running alignment!\n"
	printf "\n################################################################\n\n"
else
	printf "\nPASS: the number of fastq in ./Fastq (${fq}) matches the number of\nfastq matched to sample IDs in ${config}\n"
	printf "\n################################################################\n\n"
fi
