#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: bash bqsr_apply_make_input.sh <config prefix>
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

t_input=./Inputs/bqsr_apply.inputs-tumour
n_input=./Inputs/bqsr_apply.inputs-normal
input=./Inputs/bqsr_apply.inputs

rm -f $t_input
rm -f $n_input
rm -f $input

dict=../Reference/canFam4.dict
#contigs=$(awk '$2~/^SN/ {print $2}' $dict | sed 's/SN\://')
#contigs=($contigs)
#contigs+=( "unmapped" )

contigs=$(ls -1 ../Reference/BQSR_apply_intervals | cut -d '-' -f 1)
contigs=($contigs)
contigs+=( "unmapped" )

while read -r sampleid labid seq_center library; do
        if [[ ! ${sampleid} =~ ^#.*$ ]]
        then
                for ((i=0;i<${#contigs[@]};i++))
                do
                        if [[ $group = true ]]
                        then
                                if [[ ${labid} =~ -B$ || ${labid} =~ -N$ ]]
                                then
                                        printf "${labid},${contigs[i]},$i\n" >> $n_input
                                else
                                        printf "${labid},${contigs[i]},$i\n" >> $t_input
                                fi
                        else
                                printf "${labid},${contigs[i]},$i\n" >> $input
                        fi
                done
        fi
done < "${cohort}.config"


if [ -f $input ]
then
        tasks=`wc -l < $input`
        printf "Number of ApplyBQSR tasks to run: ${tasks}\n"
        sort -t ',' -n -k3 $input > ./Inputs/bqsr_apply_reordered.input
        mv ./Inputs/bqsr_apply_reordered.input $input
fi


if [ -f $n_input ]
then
        tasks=`wc -l < $n_input`
        printf "Number of ApplyBQSR normal sample tasks to run: ${tasks}\n"
        sort -t ',' -n -k3 $n_input > ./Inputs/bqsr_apply_reordered_normal.input
        mv ./Inputs/bqsr_apply_reordered_normal.input $n_input
fi

if [ -f $t_input ]
then
        tasks=`wc -l < $t_input`
        printf "Number of ApplyBQSR tumour sample tasks to run: ${tasks}\n"
        sort -t ',' -n -k3 $t_input > ./Inputs/bqsr_apply_reordered_tumour.input
        mv ./Inputs/bqsr_apply_reordered_tumour.input $t_input
