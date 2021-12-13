#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: this script is run by split_fastq_run_parallel.pbs
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

fqpair=`echo $1 | cut -d ',' -f 1`
file=$(basename $fqpair)
fq1=$(ls ${fqpair}*R1*f*q.gz) #Must check regex for each batch.
fq2=$(ls ${fqpair}*R2*f*q.gz)
log=./Logs/Split_fastq/${file}.log


fastp -i ${fq1} \
        -I ${fq2} \
        -AGQL \
        -w $NCPUS \
        -S 2000000 \
        -d 0 \
        --out1 ../Split_fastq/${file}_R1.fastq.gz \
        --out2 ../Split_fastq/${file}_R2.fastq.gz 2>${log}

