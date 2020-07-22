#!/bin/bash

# Create input list for indexing
# Read sample info from <cohort>.config
# Default is to make one input list for all samples (split=false), 
# as indexing is fast so the SU cost saving for splitting into two 
# jobs is negligible
# If splitting into normal and tumour jobs is desired, update group=false 
# to group=true and in that case, script assumes all non-cancer samples 
# are designated 'N' (normal) or 'B' (blood), and are lower coverage. 
# All other phenotype IDs are assigned to 'tumour'
# Provide cohort name as argument

cohort=$1

group=false

t_input=./Inputs/index.inputs-tumour
n_input=./Inputs/index.inputs-normal
input=./Inputs/index.inputs

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
	printf "Number of tumour sample index tasks to run: ${tasks}\n"
fi

if [ -f $n_input ]
then
	tasks=`wc -l < $n_input`
	printf "Number of normal sample index tasks to run: ${tasks}\n"
fi

if [ -f $input ]
then
	tasks=`wc -l < $input`
	printf "Number of index tasks to run: ${tasks}\n"
fi
