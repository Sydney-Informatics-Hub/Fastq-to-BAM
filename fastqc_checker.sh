#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: bash fastqc_checker.sh
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

# This script pulls all .oe logs from Fastq-to-BAM/Logs/FastQC/ and checks for "Analysis complete for $fastq"
# Will print list of samples that failed to $cohort_FastQC_FAILED_$(date) to Logs, otherwise will print "FastQC completed for all fastqs successfully!" to CL
# Takes fastq files from inputs, same as fastqc.sh

#variables to collect
input=./Inputs/fastqc.inputs
logs=Logs/FastQC
when=$(date +"%d_%m_%Y")

rm -f ./Logs/fastqc_checker_${when}.log

while IFS=, read fastq; do

        fastqc=`echo $fastq | cut -d ',' -f 2`
        fastq=`echo $fastq | cut -d ',' -f 2 | cut -d '/' -f 3`
        if grep -q "Analysis complete" $logs/$fastqc.oe; then
                echo "FastQC was run successfully for $fastq"
        elif grep -q "Failed to process file" $logs/$fastqc.oe; then
                echo "$fastq failed FastQC, check ./Logs/FastQC"

        fi >> ./Logs/fastqc_checker_${when}.log

done < $input

echo -e "\nFastQC checker script run successfully"
echo -e "See ./Logs/fastqc_checker_${when}.log for summary of fastQC run"
