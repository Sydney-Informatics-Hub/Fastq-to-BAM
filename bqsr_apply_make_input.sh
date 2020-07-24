#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: make inputs file for parallel exectuion of GATK ApplyBQSR
# Usage: bash bqsr_apply_mke_input.sh <cohort_name>
# Details:
# 	Job can be run as separate tumour/normal jobs, or as one job. 
#	The contigs do take longer to print for tumour compared to normal re 
#	more data to print, but the impact of input sample size on effiency is 
#	lower than for other jobs, as there are many more tasks than CPUs for 
#	this job and the walltime discrepancies among tasks are somewhat absorbed 
#	by the large number of tasks. The walltime is capped by the time to print 
#	chromosome 1, so the inputs are sorted by contig size so that the largest 
#	contigs are processed first. If no binomial grouping is desired, change 
#	group=true to group=false. Assumes all non-cancer samples have suffix '-N',
#	all other phenotype IDs are assigned to tumour.
#	Provide cohort name as argument. Sample info is read from <cohort>.config
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

group=false

t_input=./Inputs/bqsr_apply.inputs-tumour
n_input=./Inputs/bqsr_apply.inputs-normal
input=./Inputs/bqsr_apply.inputs

rm -f $t_input
rm -f $n_input
rm -f $input

dict=./Reference/hs38DH.dict
contigs=$(awk '$2~/^SN/ {print $2}' $dict | sed 's/SN\://')
contigs=($contigs)
contigs+=( "unmapped" )

awk 'NR>1' samples.config | while read CONFIG
do 
        labSampleID=`echo $CONFIG | cut -d ' ' -f 2`
	for ((i=0;i<${#contigs[@]};i++))	
	do
		if [[ $group = true ]]
		then
			if [[ $labSampleID = *-N ]]
       			then
                		printf "${labSampleID},${contigs[i]},$i\n" >> $n_input
        		else
                		printf "${labSampleID},${contigs[i]},$i\n" >> $t_input
        		fi
		else			
			printf "${labSampleID},${contigs[i]},$i\n" >> $input 
		fi      	
	done
done 


if [ -f $input ]
then
        tasks=`wc -l < $input`
        printf "Number of ApplyBQSR tasks to run: ${tasks}\n"
	sort -t ',' -n -k3 $input > ./Inputs/bqsr_apply_reordered.input
	mv ./Inputs/bqsr_apply_reordered.input $input
fi


if [ -f $n_input ]
then
        tasks=`wc -l < $n_input`
        printf "Number of ApplyBQSR normal sample tasks to run: ${tasks}\n"
	sort -t ',' -n -k3 $n_input > ./Inputs/bqsr_apply_reordered_normal.input
	mv ./Inputs/bqsr_apply_reordered_normal.input $n_input
fi

if [ -f $t_input ]
then
        tasks=`wc -l < $t_input`
        printf "Number of ApplyBQSR tumour sample tasks to run: ${tasks}\n"
	sort -t ',' -n -k3 $t_input > ./Inputs/bqsr_apply_reordered_tumour.input
	mv ./Inputs/bqsr_apply_reordered_tumour.input $t_input
fi
