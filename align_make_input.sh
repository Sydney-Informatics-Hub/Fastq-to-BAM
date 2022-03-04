#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: bash align_make_input.sh <config prefix>
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
inputs=./Inputs/align.inputs
platform=illumina

rm -f $inputs

while read -r sampleid labid seqcentre lib
do
        # Ignore headers starting with #
        if [[ ! ${sampleid} =~ ^#.*$ ]]
        then
                # Set lib=1 if not defined
                if [ ! "$lib" ]
                then
                        lib=1
                fi

                #Find all fastq input pairs for $sample
                #Will collect multiplexed samples named with standard FASTQ naming convention
                #check regex per batch
                fqpairs=$(ls ../Fastq/*${sampleid}*f*q.gz | sed  's/_R1.*\|_R2.*\|_R1_*\|_R2_*\|.R1.*\|.R2.*//' | sort | uniq)
                fqpairs=($fqpairs)


                #Get the flowcell and lane info for each unique FASTQ pair
                for ((i=0; i<${#fqpairs[@]}; i++))
                do
                        lib_diff=''
                        fq1=$(ls ${fqpairs[i]}*R1*.f*q.gz) #Must check regex for each batch.
                        flowcell=$(zcat ${fq1} | head -1 | cut -d ':' -f 3)
                        lane=$(zcat ${fq1} | head -1 | cut -d ':' -f 4)

                        #Print each of the split chunks with flowcell and lane info to inputs file:
                        set=$(basename ${fqpairs[i]})
                        splitpairs=$(ls ../Split_fastq/*[0-9].${set}* | sed 's/_R1.*\|_R2.*\|_R1_*\|_R2_*\|.R1.*\|.R2.*//' | uniq)

                        splitpairs=($splitpairs)

                        for ((c=0; c<${#splitpairs[@]}; c++))
                        do
                                if [ "$lib_diff" ]
                                then
                                        printf "${splitpairs[c]},${labid},${seqcentre},${lib_diff},${platform},${flowcell},${lane}\n" >> $inputs
                                else
                                        printf "${splitpairs[c]},${labid},${seqcentre},${lib},${platform},${flowcell},${lane}\n" >> $inputs
                                fi
                        done
                done
        fi
done < "${cohort}.config"

tasks=`wc -l < $inputs`
printf "Number of alignment tasks to run: ${tasks}\n"
printf "Before running align_run_parallel.pbs, edit the script and unhash the correct task script for human/non-human samples\n"
