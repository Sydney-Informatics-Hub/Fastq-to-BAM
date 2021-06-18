#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: this script is run by merge_align_run_parallel.pbs
# Version: 2.0
#
# For more details see: https://github.com/Sydney-Informatics-Hub/Fastq-to-BAM
#
# If you use this script towards a publication, please acknowledge the
# Sydney Informatics Hub (or co-authorship, where appropriate).
#
# Suggested acknowledgement:
# The authors acknowledge the support provided by the Sydney Informatics Hub,
# a Core Research Facility of the University of Sydney. This research/project
# was undertaken with the assistance of resources and services from the National
# Computational Infrastructure (NCI), which is supported by the Australian
# Government, and the Australian BioCommons which is enabled by NCRIS via
# Bioplatforms Australia funding.
#
##########################################################################

sample=`echo $1 | cut -d ',' -f 1`
labSampleID=`echo $1 | cut -d ',' -f 2`
err=./Logs/Merge_align_error_capture/${labSampleID}.err

cd ../Align_split

#Sample names must be unique
bams=$(find . -name "*${sample}*.bam" | xargs echo) #check regex for each batch

\rm -rf ../Merge_align/${labSampleID}.merged.nameSorted.bam   #attempt to set stripe on existingfile will cause fatal error
lfs setstripe -c 15 ../Merge_align/${labSampleID}.merged.nameSorted.bam

sambamba merge -t $NCPUS ../Merge_align/${labSampleID}.merged.nameSorted.bam $bams

if ! samtools quickcheck ../Merge_align/${labSampleID}.merged.nameSorted.bam
then
        printf "Corrupted or missing BAM\n" > $err
fi
