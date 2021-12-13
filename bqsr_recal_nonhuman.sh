#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: this script is run by bqsr_recal_run_parallel.pbs
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


labSampleID=`echo $1 | cut -d ',' -f 1`
intNum=`echo $1 | cut -d ',' -f 2`
interval=`echo $1 | cut -d ',' -f 3`

echo $labSampleID $intNum $interval

bam=../Dedup_sort/${labSampleID}.coordSorted.dedup.bam
log=./Logs/BQSR_recal/${labSampleID}.${intNum}.log
err=./Logs/BQSR_recal_error_capture/${labSampleID}.${intNum}.err
out=../BQSR_tables/${labSampleID}.${intNum}.recal_data.table

ref=
dbsnp= # edit this depending on your reference

gatk BaseRecalibrator \
        --java-options "-Xmx20g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
        -R $ref \
        -L $interval \
        -I $bam  \
        --known-sites $dbsnp \
        -O $out >> $log 2>&1


if grep -q -i error $log
then
        printf "Error in GATK log ${log}\n" >> $err
fi

if grep -q Exception $log
then
        printf "Exception in GATK log ${log}\n" >> $err
fi
