#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: make inputs file for parallel exectuion of fastp fastq splitting
# Usage: bash split_fastq_make_input.sh
# Details:
# 	To generate split_fastq.input, run 'bash split_fastq_make_input.sh'
# 	AFTER checking that the regex to find the fastq files matches the regex in
#	this script, and in split_fastq.sh (eg if your data ends in something 
#	other than fq.gz or fastq.gz, or if it does not use _R1_/_R2_, _R1./_R2., 
#	or .R1./.R2. to designate paired files, you will need to edit the regexes.
#	Assumes all input data is paired and resides in a single directory 
#	named './Fastq'. After running, check the number of inputs printed matches
#	your number of fastq pairs. For samples with multiple pairs, please do
#	not concatenate them into one massive pair per sample. That is general
#	advice that extends beyond the scope of this analysis :-)
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

inputs=./Inputs/split_fastq.inputs

mkdir -p ./Inputs

rm -f $inputs

ls ./Fastq/*.f*q.gz | sed 's/_R1.*\|_R2.*\|_R1_*\|_R2_*\|.R1.*\|.R2.*//' | uniq > $inputs

tasks=`wc -l < $inputs`
printf "Number of fastq pairs to split: ${tasks}\n"
