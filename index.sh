#!/bin/bash

# Index the dedup/sort bam file
# Running index as a separate job re high charge rate for dedup/sort job

module load samtools/1.10

labSampleID=`echo $1 | cut -d ',' -f 1`
NCPUS=`echo $1 | cut -d ',' -f 2`

bam_out=./Dedup_sort/${labSampleID}.coordSorted.dedup.bam

samtools index -@ $NCPUS $bam_out

