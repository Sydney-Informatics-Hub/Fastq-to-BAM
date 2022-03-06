#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: this script is run by bqsr_merge_run_parallel.pbs
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

set -e

labSampleID=$1

jvm=10 # Set this variable to 28 for hugemem queue or 10 for normal queue

log=./Logs/BQSR_merge/${labSampleID}.log
err=./Logs/BQSR_merge_error_capture/${labSampleID}.err
bam_out=../Final_bams/${labSampleID}.final.bam
ref=

\rm -rf $bam_out #attempt to set stripe on existing file will cause fatal error
lfs setstripe -c 15 $bam_out

gatk GatherBamFiles \
        --java-options "-Xmx${jvm}G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
        -I ./Inputs/BQSR_merge_lists/${labSampleID}.list \
        -O $bam_out \
        -R $ref \
        --CREATE_INDEX true >> $log 2>&1
        #\
        #--CREATE_MD5_FILE true #Do not create md5 in this job, as it will cost 3x SU

if ! samtools quickcheck $bam_out
then
        printf "Corrupted or missing BAM\n" > $err
fi

if grep -q -i error $log
then
        printf "Error in GATK log ${log}\n" >> $err
fi

if grep -q Exception $log
then
        printf "Exception in GATK log ${log}\n" >> $err
fi

