#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: align split read pair files in parallel to non-human
# reference genome with BWA-mem 
# Usage: this script is executed by align_run_parallel.pbs
# Details:
# 	Parallel alignment of  split fastq pairs with BWA mem
# 	Output is read name sorted for compatibility with sambmaba merge
# 	Check that the regex at 'ls' matches your data, edit as required.
#	Install/edit paths to seqtk and bwa-postalt.js
# Author: Cali Willet
# cali.willet@sydney.edu.au
# Date last modified: 26/08/2020
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
 

set -e

ref=Reference/<ref>

module load samtools/1.10
module load bwa/0.7.17

fqpair=`echo $1 | cut -d ',' -f 1` 
fq1=$(ls ${fqpair}*R1_*.f*q.gz) #Must check regex for each batch
fq2=$(ls ${fqpair}*R2_*.f*q.gz)
sampleID=`echo $1 | cut -d ',' -f 2`
centre=`echo $1 | cut -d ',' -f 3`
lib=`echo $1 | cut -d ',' -f 4`
platform=`echo $1 | cut -d ',' -f 5`
flowcell=`echo $1 | cut -d ',' -f 6`
lane=`echo $1 | cut -d ',' -f 7` 	

outPrefix=$(basename $fqpair)
bam_out=./Align_split/${outPrefix}.bam
log=./Logs/BWA/${outPrefix}.log
err=./Error_capture/Align_split/${outPrefix}.err

bwa mem \
        -M -t $NCPUS $ref \
        -R "@RG\tID:${flowcell}.${lane}_${sampleID}_${lib}\tPL:${platform}\tPU:${flowcell}.${lane}\tSM:${sampleID}\tLB:${sampleID}_${lib}\tCN:${centre}" \
        $fq1 $fq2 2> $err \
        | samtools sort \
        -n -@ $NCPUS \
        -o $bam_out
 
# Error checking:
if ! samtools quickcheck ${bam_out}
then 
	printf "Corrupted or missing BAM\n" > $err  
fi

test=$(tail -1 ${log} | awk '$0~/^\[main\] Real time/') # This is correct for BWA version 0.7.15 and 0.7.17 - must check for new versions
if [[ ! $test ]]
then 
        printf "Error in BWA log\n" >> ${err}
fi

if ! grep -q "M::mem_process_seqs" ${log} # This is correct for BWA version 0.7.15 and 0.7.17 - must check for new versions
then 
        printf "Error in BWA log\n" >> ${err}
fi
