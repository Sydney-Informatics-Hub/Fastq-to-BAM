#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: this script is run by reformat_split_disc_run_parallel.pbs
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

sample=$1

cd ../SplitDisc
temp=Temp_${sample}
mkdir -p $temp

in=( ".split" ".disc")

for (( i = 0; i < ${#in[@]}; i ++ ))
do
        samtools sort -@ $NCPUS \
                -o ${temp}/${sample}${in[$i]}.sort.bam \
                ${sample}${in[$i]}.sam

        samtools index -@ $NCPUS ${temp}/${sample}${in[$i]}.sort.bam

        mv ${temp}/* .

done

rmdir $temp
