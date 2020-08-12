#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: calculate base and alignment metrics for final BAMs
# Usage: this script is executed with by bam_metrics_run_parallel.pbs
# Details:
#	Calculate the following metrics for final BAMs:
#	alignment_summary_metrics
#	insert_size_metrics
#	quality_distribution_metrics
#	quality_by_cycle_metrics
#	base_distribution_by_cycle_metrics
#	A combination of TSV and PDF are created, 9 files per sample 
# Author: Cali Willet
# cali.willet@sydney.edu.au
# Date last modified: 12/08/2020
#
# If you use this script towards a publication, please acknowledge the
# Sydney Informatics Hub (or co-authorship, where appropriate).
#
# Suggested acknowledgement:
# The authors acknowledge the scientific and technical assistance 
# <or e.g. bioinformatics assistance of <PERSON>> of Sydney Informatics
# Hub and resources and services from the National Computational 
# Infrastructure (NCI), which is supported by the Australian Government
# with access facilitated by the University of Sydney.
# 
#########################################################

module load gatk/4.1.2.0 #loads also java/jdk1.8.0_60
module load python3/3.7.4
module load R/3.6.1

sample=$1
	
bam_in=./Final_bams/${sample}.final.bam
out=./Final_bams_metrics/${sample}
ref=./Reference/ARS-UCD1.2_Btau5.0.1Y.fa
log=./Bam_metrics_logs/${sample}.log

gatk  CollectMultipleMetrics \
        --java-options "-Xmx120G " \
        -R $ref \
        -I $bam_in  \
        -O $out >> $log 2>&1
	

