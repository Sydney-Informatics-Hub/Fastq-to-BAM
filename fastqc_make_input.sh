#! /bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: make input for fastQC
# Usage: bash fastqc_make_inputs.sh
# Details:
# 	Create input file for fastqc_run_parallel.pbs
# 	Sort by input file size largest to smallest 
# 	VIP !!! Check that the find pattern matches, eg fastq/fq, gz vs not gz
# 	VIP !!! Manually verify that the inputs list is the same line length
# 	as the number of files in the fastq directory

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


inputs=./Inputs/fastqc.inputs

find ./Fastq -name "*fastq*" -print | xargs ls -l | sort -rnk 5 | awk '{print $9}'> $inputs

tasks=`wc -l < $inputs`

printf "Number of fastQC tasks to run: ${tasks}\n"

