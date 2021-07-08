#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: bash merge_align_make_input.sh <config prefix>
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

t_input=./Inputs/merge_align.inputs-tumour
n_input=./Inputs/merge_align.inputs-normal
input=./Inputs/merge_align.inputs

rm -f $t_input
rm -f $n_input
rm -f $input

while read -r sampleid labid seqcentre lib
do
        # Ignore headers starting with #
        if [[ ! ${sampleid} =~ ^#.*$ ]]
        then
                if [[ $group = true ]]
                then
                        if [[ $labid = -B$ || $labid =~ -N$ ]]
                        then
                                printf "${sampleid},${labid}\n" >> $n_input
                        else
                                printf "${sampleid},${labid}\n" >> $t_input
                        fi
                else
                        printf "${sampleid},${labid}\n" >> $input
                fi
        fi
done < "${cohort}.config"


if [ -f $t_input ]
then
        tasks=`wc -l < $t_input`
        printf "Number of tumour sample merge tasks to run: ${tasks}\n"
fi

if [ -f $n_input ]
then
        tasks=`wc -l < $n_input`
        printf "Number of normal sample merge tasks to run: ${tasks}\n"
fi

if [ -f $input ]
then
        tasks=`wc -l < $input`
        printf "Number of merge tasks to run: ${tasks}\n"
fi
