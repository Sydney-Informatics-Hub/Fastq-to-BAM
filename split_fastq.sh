#!/bin/bash 

#Split each fastq pair into 2,000,000 lines (500,000 reads)

fastp=/scratch/<project>/apps/fastp #or replace with module load if appropriate

fqpair=`echo $1 | cut -d ',' -f 1`
file=$(basename $fqpair)
NCPUS=`echo $1 | cut -d ',' -f 2`

fq1=$(ls ${fqpair}*R1.f*q.gz) #Must check regex for each batch.
fq2=$(ls ${fqpair}*R2.f*q.gz) 

log=./Logs/Fastp/${file}.log

$fastp -i ${fq1} \
	-I ${fq2} \
	-AGQL \
	-w $NCPUS \
	-S 2000000 \
	-d 0 \
	--out1 Fastq_split/${file}_R1.fastq.gz \
	--out2 Fastq_split/${file}_R2.fastq.gz 2>${log}
