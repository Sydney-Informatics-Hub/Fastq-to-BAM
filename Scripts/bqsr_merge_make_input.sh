#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: bash bqsr_merge_make_input.sh <config prefix>
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
dict=
contigs=$(awk '$2~/^SN/ {print $2}' $dict | sed 's/SN\://')
contigs=($contigs)
contigs+=( "unmapped" )
num_recal=$((${#contigs[@]}-1)) # Count starts at 0

group=false

t_input=./Inputs/bqsr_merge.inputs-tumour
n_input=./Inputs/bqsr_merge.inputs-normal
input=./Inputs/bqsr_merge.inputs

rm -f $t_input
rm -f $n_input
rm -f $input

# Make BQSR merge list
# List of bam files to merge per sample
mkdir -p ./Inputs/BQSR_merge_lists

while read -r sampleid labid seq_center library; do
        if [[ ! ${sampleid} =~ ^#.*$ ]]; then
                list=./Inputs/BQSR_merge_lists/${labid}.list
                \rm -rf $list
                for i in $(seq 0 ${num_recal})
                do
                        printf "../BQSR_apply/${labid}.${i}.recal.bam\n" >> $list
                done
        fi
done < "${cohort}.config"

# Make input list
while read -r sampleid labid seq_center library; do
        if [[ ! ${sampleid} =~ ^#.*$ ]]
        then
                if [[ $group = true ]]
                then
                        if [[ ${labid} =~ -B$ || ${labid} =~ -N$ ]]
                        then
                                printf "${labid}\n" >> $n_input
                        else
                                printf "${labid}\n" >> $t_input
                        fi
                else
                        printf "${labid}\n" >> $input
                fi
        fi
done < "${cohort}.config"

if [ -f $input ]
then
        tasks=`wc -l < $input`
        printf "Number of BQSR merge tasks to run: ${tasks}\n"
fi

if [ -f $n_input ]
then
        tasks=`wc -l < $n_input`
        printf "Number of BQSR merge normal sample tasks to run: ${tasks}\n"
fi

if [ -f $t_input ]
then
        tasks=`wc -l < $t_input`
        printf "Number of BQSR merge tumour sample tasks to run: ${tasks}\n"
fi
