#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: run GATK GatherBQSRReports over parallel tasks
# Usage: this script is executed by bqsr_gather_run_parallel.pbs
# Author: Cali Willet
# cali.willet@sydney.edu.au
# Date last modified: 24/07/2020
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

labSampleID=$1

tables=$(find ./BQSR_tables -name "${labSampleID}.*.recal_data.table")
tables=($tables)
tables_list=''

for file in ${tables[@]}
do 
	tables_list+=" -I $file"
done

gatk GatherBQSRReports \
	$tables_list \
	-O ./BQSR_tables/${labSampleID}.recal_data.table
