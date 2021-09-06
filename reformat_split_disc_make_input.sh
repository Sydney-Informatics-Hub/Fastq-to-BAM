#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: reformat_split_disc_make_input.sh <cohort prefix>
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

group=false

t_input=./Inputs/reformat_split_disc.inputs-tumour
n_input=./Inputs/reformat_split_disc.inputs-normal
input=./Inputs/reformat_split_disc.inputs

rm -f $t_input
rm -f $n_input
rm -f $input

awk 'NR>1' ${cohort}.config | while read CONFIG
do
        labSampleID=`echo $CONFIG | cut -d ' ' -f 2`

        if [[ $group = true ]]
        then
                if [[ $labSampleID = *-N ]]
                then
                        printf "${labSampleID}\n" >> $n_input
                else
                        printf "${labSampleID}\n" >> $t_input
                fi
        else
                printf "${labSampleID}\n" >> $input
        fi
done


if [ -f $input ]
then
        tasks=`wc -l < $input`
        printf "Number of reformat split/disc tasks to run: ${tasks}\n"
fi

if [ -f $n_input ]
then
        tasks=`wc -l < $n_input`
        printf "Number of reformat split/disc normal sample tasks to run: ${tasks}\n"
fi

if [ -f $t_input ]
then
        tasks=`wc -l < $t_input`
        printf "Number of reformat split/disc tumour sample tasks to run: ${tasks}\n"
fi
