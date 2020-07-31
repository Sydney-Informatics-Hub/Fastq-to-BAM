#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: make inputs for parallel concatenation of HLA gene fastq files
# split parallel alignment into one HLA gene per sample
# Usage: bash cat_HLA_make_input.sh <cohort_name> 
# Details:
#       For tumour/normal, run time is ~ half for normal so you may split
# 	this into 2 jobs if preferred. If so, this script assumes all 
#	non-cancer samples are designated 'N' ('normal'), and are lower 
#	coverage. All other phenotype IDs are assigned to 'tumour'
#       If no binomial grouping is desired, change group=true to group=false
#       Provide cohort name as argument. Sample info is read from <cohort>.config
#
# Author: Cali Willet
# cali.willet@sydney.edu.au
# Date last modified: 31/07/2020
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

group=false

t_input=./Inputs/cat_HLA.inputs-tumour
n_input=./Inputs/cat_HLA.inputs-normal
input=./Inputs/cat_HLA.inputs

rm -f $t_input
rm -f $n_input
rm -f $input

awk 'NR>1' ${cohort}.config | while read CONFIG
do 
        labSampleID=`echo $CONFIG | cut -d ' ' -f 2`
	sample=`echo $CONFIG | cut -d ' ' -f 1`
	
        if [[ $group = true ]]
        then
                if [[ $labSampleID = *-N ]]
                then
                        printf "${sample},${labSampleID}\n" >> $n_input
                else 
                        printf "${sample},${labSampleID}\n" >> $t_input
                fi
        else
                printf "${sample},${labSampleID}\n" >> $input
        fi
done 


if [ -f $input ]
then
        tasks=`wc -l < $input`
        printf "Number of cat HLA tasks to run: ${tasks}\n"
fi

if [ -f $n_input ]
then
        tasks=`wc -l < $n_input`
        printf "Number of cat HLA normal sample tasks to run: ${tasks}\n"
fi

if [ -f $t_input ]
then
        tasks=`wc -l < $t_input`
        printf "Number of cat HLA tumour sample tasks to run: ${tasks}\n"
fi


