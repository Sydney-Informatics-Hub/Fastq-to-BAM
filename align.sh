#!/bin/bash

# Parallel alignment of  split fastq pairs with BWAkit
# This is modified BWAkit to include '-m' flag and for custom sort
# Custom sort is required for this pipeline, because BWAkit default 
# sort order  is queryname, but does not include the SO tag in the 
# header which is required by sambamba (used to merge the 
# parallel-aligned BAMs). BWAkit includes an option to sort by 
# coordinate with samtools but does not allow the user to specify 
# sort by name (-n samtools flag) which would give the SO header tag.
# Below code is the teased out run-bwamem from the bwakit package,
# with added -M and sort -n flags 

set -e

ref=Reference/hs38DH.fasta

module load samtools/1.10
module load bwa/0.7.17
module load k8/0.2.5
seqtk=/scratch/<project>/apps/seqtk #Compiled for Gadi
postalt=/scratch/<project>/apps/bwa-postalt.js

outdir=./Align_split
hladir=./HLA_fastq
logdir=./Logs/BWA
errdir=./Error_capture/Align_split

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
err=${errdir}/${outPrefix}.err
log=${logdir}/${outPrefix}.log.bwamem

#bwakit emits default sort order (queryname, but no SO tag in headers) or option to sort 
#by coordinate with samtools but does not allow the -n flag to specify sort order by name.
#Sambamba requires queryname sorted with SO tag. Below code does it by force. 


${seqtk} mergepe $fq1 $fq2 \
        | bwa mem -p -t $NCPUS \
        -R "@RG\tID:${flowcell}.${lane}_${sampleID}_${lib}\tPL:${platform}\tPU:${flowcell}.${lane}\tSM:${sampleID}\tLB:${sampleID}_${lib}\tCN:${centre}" \
        -M ${ref} - 2> ${log} \
        | k8 ${postalt} \
        -p ${hladir}/${outPrefix} ${ref}.alt \
        | samtools sort -n -@ $NCPUS -o ${outdir}/${outPrefix}.aln.bam

# Error checking:
if ! samtools quickcheck ${outdir}/${outPrefix}.aln.bam
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
