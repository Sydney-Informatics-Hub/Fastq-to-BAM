#!/bin/bash

# Create input list for parallel alignment with BWAkit
# Read sample info from <cohort>.config
# Unsplit fastq in Fastq directory, split fastq in Fastq_split
# Assumes flowcell ID is field 3 of ':' delim read ID
# Assumes flowcell lane number is field 4 of ':' delim read ID
# Assumes platform is Illumina
# Assumes library ID is either defined in samples.config, or 
# that the default of '1' is aceptable
# Assumes one library per sample - an example of 2 libraries 
# per sample is hashed out within the while loop, where these samples
# were given the value of '2' in the samples.config library field
# and were assigned 1a and 1b libary values based on '_a_' and '_b_'
# in their field name. 
# Provide cohort name as argument

cohort=$1

inputs=./Inputs/align.inputs

rm -f $inputs


awk 'NR>1' ${cohort}.config | while read LINE
do 
	sample=`echo $LINE | cut -d ' ' -f 1`
	labSampleID=`echo $LINE | cut -d ' ' -f 2`
	centre=`echo $LINE | cut -d ' ' -f 3`
	lib=`echo $LINE | cut -d ' ' -f 4`
	platform=illumina
	
	if [ ! "$lib" ]
	then
	    	lib=1
	fi
	
	#Find all fastq input pairs for $sample
	fqpairs=$(ls ./Fastq/${sample}_*f*q.gz | sed  's/_R1.*\|_R2.*\|_R1_*\|_R2_*\|.R1.*\|.R2.*//' | sort | uniq) #check regex per batch
	fqpairs=($fqpairs)
	
		
	#Get the flowcell and lane info from the original pairs:	
	for ((i=0; i<${#fqpairs[@]}; i++))
	do
		lib_diff=''
		fq1=$(ls ${fqpairs[i]}*R1.f*q.gz) #Must check regex for each batch.
		flowcell=$(zcat ${fq1} | head -1 | cut -d ':' -f 3)
		lane=$(zcat ${fq1} | head -1 | cut -d ':' -f 4)	
		
		#if [[ $lib == 2 ]] # this is a custom use case for OSCC batch 2, where some samples have 2 libraries ('_a_' and '_b_')
		#then 
			#if [[ $fq1 =~ _a_ ]]
			#then 
				#lib_diff=1a
			#else 
				#lib_diff=1b
			#fi	
		#fi		
		
		#Print each of the split chunks with flowcell and lane info to inputs file:
		set=$(basename ${fqpairs[i]})
		splitpairs=$(ls Fastq_split/*${set}* | sed 's/_R1.*\|_R2.*\|_R1_*\|_R2_*\|.R1.*\|.R2.*//' | uniq)
				
		splitpairs=($splitpairs)

		for ((c=0; c<${#splitpairs[@]}; c++))
		do
			if [ "$lib_diff" ]
			then
				printf "${splitpairs[c]},${labSampleID},${centre},${lib_diff},${platform},${flowcell},${lane}\n" >> align.input	
			else
				printf "${splitpairs[c]},${labSampleID},${centre},${lib},${platform},${flowcell},${lane}\n" >> align.input
			fi	
		done			
	done			
done	

tasks=`wc -l < $inputs`
printf "Number of alignment tasks to run: ${tasks}\n"
	
