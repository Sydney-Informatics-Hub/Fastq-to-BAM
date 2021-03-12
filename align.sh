#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: align split read pair files in parallel to hs38+alt
# reference genome with BWA-mem including post-processing
# Usage: this script is executed by align_run_parallel.pbs
# Details:
# 	Parallel alignment of  split fastq pairs with BWAkit
# 	This is modified BWAkit to include '-m' flag and for custom sort
# 	Custom sort is required for this pipeline, because BWAkit default 
# 	sort order  is queryname, but does not include the SO tag in the 
#	 header which is required by sambamba (used to merge the 
# 	parallel-aligned BAMs). BWAkit includes an option to sort by 
# 	coordinate with samtools but does not allow the user to specify 
# 	sort by name (-n samtools flag) which would give the SO header tag.
# 	Below code is the teased out run-bwamem from the bwakit package,
# 	with added -M and sort -n flags.
# 	Check that the regex at 'ls' matches your data, edit as required.
#	Install/edit paths to seqtk and bwa-postalt.js
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
 

set -e

ref=Reference/hs38DH.fasta

module load samtools/1.10
module load bwa/0.7.17
module load k8/0.2.5
seqtk=/scratch/<project>/apps/seqtk #Compiled for Gadi
postalt=/scratch/<project>/apps/bwa-postalt.js

fqpair=`echo $1 | cut -d ',' -f 1` 
fq1=$(ls ${fqpair}*R1.f*q.gz) #Must check regex for each batch
fq2=$(ls ${fqpair}*R2.f*q.gz)
sampleID=`echo $1 | cut -d ',' -f 2`
centre=`echo $1 | cut -d ',' -f 3`
lib=`echo $1 | cut -d ',' -f 4`
platform=`echo $1 | cut -d ',' -f 5`
flowcell=`echo $1 | cut -d ',' -f 6`
lane=`echo $1 | cut -d ',' -f 7`
NCPUS=`echo $1 | cut -d ',' -f 8` 	

outPrefix=$(basename $fqpair)
bam_out=./Align_split/${outPrefix}.aln.bam
hla_out=./HLA_fastq/${outPrefix}
log=./Logs/BWA/${outPrefix}.log.bwamem
err=./Error_capture/Align_split/${outPrefix}.err

rm -rf $log $err ${bam_out}* $hla_out

#bwakit emits default sort order (queryname, but no SO tag in headers) or option to sort 
#by coordinate with samtools but does not allow the -n flag to specify sort order by name.
#Sambamba requires queryname sorted with SO tag. Below code does it by force. 


${seqtk} mergepe $fq1 $fq2 \
        | bwa mem -p -t $NCPUS \
        -R "@RG\tID:${flowcell}.${lane}_${sampleID}_${lib}\tPL:${platform}\tPU:${flowcell}.${lane}\tSM:${sampleID}\tLB:${sampleID}_${lib}\tCN:${centre}" \
        -M ${ref} - 2> ${log} \
        | k8 ${postalt} \
        -p ${hla_out} ${ref}.alt \
        | samtools sort -n -@ $NCPUS -o ${bam_out}

# Error checking:
if ! samtools quickcheck ${bam_out}
then 
	printf "Corrupted or missing BAM\n" > $err  
fi

test=$(awk '$0~/^\[main\] Real time/' ${log}) # This is correct for BWA version 0.7.15 and 0.7.17 - must check for new versions
if [[ ! $test ]]
then 
        printf "Error in BWA log\n" >> ${err}
fi

if ! grep -q "M::mem_process_seqs" ${log} # This is correct for BWA version 0.7.15 and 0.7.17 - must check for new versions
then 
        printf "Error in BWA log\n" >> ${err}
fi
