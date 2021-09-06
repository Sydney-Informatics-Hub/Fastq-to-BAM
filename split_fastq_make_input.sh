#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: bash split_fastq_make_input.sh <config prefix>
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

cohort=../$1

inputs=./Inputs/split_fastq.inputs

rm -f $inputs

awk 'NR>1' ${cohort}.config | while read sample
        do
        sample=`echo $sample | cut -f 1 -d ' '`
        #echo $sample #use this to check $sample is correctly capturing your sampleID

        # find fastq pairs for each $sample. Check regex matches your samples
        fqpairs=$(find ../Fastq -name ${sample}*.f*q.gz | sed  's/.R1.*\|.R2.*\|_R1.*\|_R2.*\|_R1_*\|_R2_*//' $
        #echo $fqpairs #use this to check $fqpairs is capturing your fqs

        # print fq pair info to input file
        printf '%s\n' $fqpairs >> $inputs

        done
tasks=`wc -l < ${inputs}`
printf "Number of fastq pairs to split: ${tasks}\n"
