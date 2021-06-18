#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: bash merge_align_check_fastq_input_vs_bam_output_sizes.sh <config prefix>
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

output=./Logs/Merge_align/check_fastq_vs_bam_sizes.txt

if [ "$cohort" ]
then
        output=./Logs/Merge_align/${cohort}_check_fastq_vs_bam_sizes.txt
fi
rm -rf $output

printf "#ID\tFastq\tAlignMerged\tRatio\n" >> $output

while read -r sampleid labid seqcentre lib
do
        # Ignore headers starting with #
        if [[ ! ${sampleid} =~ ^#.*$ ]]
        then
                fq_size=$(find  ../Fastq -type f -name *${sampleid}* -print | xargs ls -lh | awk '{print $5}' | sed 's/G$//'| awk '{s+=$1} END {print s}')
                bam_size=$(ls -lh ../Merge_align/${labid}.merged.nameSorted.bam | awk '{print $5}' | sed 's/G$//')
                ratio=$(bc -l <<< "scale=2; ${bam_size}/${fq_size}")
                printf "${labid}\t${fq_size}\t${bam_size}\t${ratio}\n" >> $output
        fi
done < "${cohort}.config"

