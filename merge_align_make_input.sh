#!/bin/bash

# Create input list for merging parallel alignment output
# Read sample info from <cohort>.config
# Assumes all non-cancer samples are designated 'N' (normal)
# or 'B' (blood), and are lower coverage. All other phenotype 
# IDs are assigned to 'tumour'. Splitting in this way is 
# recommended when there is a binomial difference in coverage, 
# eg 30X/60X cancer projects. For cohorts where the input data 
# varies by a large degree and is not really binomial, the 
# grouping can be done on file size (see 'merge_align_make_input_by_fastq_size.sh')
# If no binomial grouping is desired, change group=true to group=false
# Provide cohort name as argument

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
