#!/bin/bash

# Merge the split BAMs for each sample into sample-level, queryname sorted BAMs

module load samtools/1.10
module load sambamba/0.7.1

sample=`echo $1 | cut -d ',' -f 1` 
labSampleID=`echo $1 | cut -d ',' -f 2` 
NCPUS=`echo $1 | cut -d ',' -f 3`


cd ./Align_split


err=../Error_capture/Align_merge/${labSampleID}.err
bams=$(find . -name "*${sample}*.bam" | xargs echo) #Check regex per batch - for CSCC, this is OK, all sample names are unique with this simple regex


\rm -rf ../Align_merged/${labSampleID}.merged.nameSorted.bam   #attempt to set stripe on existingfile will cause fatal error
lfs setstripe -c 15 ../Align_merged/${labSampleID}.merged.nameSorted.bam


sambamba merge -t $NCPUS ../Align_merged/${labSampleID}.merged.nameSorted.bam $bams


if ! samtools quickcheck ../Align_merged/${labSampleID}.merged.nameSorted.bam
then 
	printf "Corrupted or missing BAM\n" > $err  
fi
