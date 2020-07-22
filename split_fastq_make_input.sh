#!/bin/bash

# Create input text file for parallel fast splitting. Check the regex matches
# your data. Check the number of inputs printed matches your number of fastq pairs

inputs=./Inputs/split_fastq.inputs

rm -f $inputs

ls ./Fastq/*.f*q.gz | sed 's/_R1.*\|_R2.*\|_R1_*\|_R2_*\|.R1.*\|.R2.*//' | uniq > $inputs

tasks=`wc -l < $inputs`
printf "Number of fastq pairs to split: ${tasks}\n"
