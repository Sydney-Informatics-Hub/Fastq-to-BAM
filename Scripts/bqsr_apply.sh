#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: this script is run by bqsr_apply_run_parallel.pbs
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

# special chracters in contig names requires this specific method to read in positional arguments:
IFS=','
params=($1)
labSampleID=${params[0]}
contig=${params[1]}
counter=${params[2]}

ref=../Reference/canFam4.fasta
table=../BQSR_tables/${labSampleID}.recal_data.table
bam_in=../Dedup_sort/${labSampleID}.coordSorted.dedup.bam
bam_out=../BQSR_apply/${labSampleID}.${counter}.recal.bam
log=./Logs/BQSR_apply/${labSampleID}.${counter}.log
err=./Logs/BQSR_apply_error_capture/${labSampleID}.${counter}.err

gatk ApplyBQSR \
        --java-options "-Xmx6G -Dsamjdk.compression_level=5 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
        -R $ref \
        -I $bam_in  \
        --bqsr-recal-file $table \
        -L $contig \
        -O $bam_out >> $log 2>&1


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

