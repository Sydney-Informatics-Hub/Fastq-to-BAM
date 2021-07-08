#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: bash dedup_sort_make_input.sh <config prefix>
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

t_input=./Inputs/dedup_sort.inputs-tumour
n_input=./Inputs/dedup_sort.inputs-normal
input=./Inputs/dedup_sort.inputs

rm -f $t_input
rm -f $n_input
rm -f $input

group=false

# Collect labids from samples.config. Assumes each labid is unique
# For cancer projects, assumes normal lab ids are appended with -N or -B.
# Everything else goes in Tumour group, but downstream scripts use -T, -M or -P appended to labid to indicate tumour
while read -r sampleid labid seq_center library; do
        if [[ ! ${sampleid} =~ ^#.*$ ]]; then
                if [[ $group =~ true ]]; then
                       if [[ $labid =~ *-N|B ]]; then
                                echo $labid >> $n_input
                        else
                                echo $labid >> $t_input
                        fi
                else
                        echo $labid >> $input
                        samples+=("${labid}")
                fi
        fi
done < "${cohort}.config"

if [ -f $t_input ]
then
        tasks=`wc -l < $t_input`
        printf "Number of tumour sample tasks to run: ${tasks}\n"
fi

if [ -f $n_input ]
then
        tasks=`wc -l < $n_input`
        printf "Number of normal sample tasks to run: ${tasks}\n"
fi

if [ -f $input ]
then
        tasks=`wc -l < $input`
        printf "Number of tasks to run: ${tasks}\n"
fi

