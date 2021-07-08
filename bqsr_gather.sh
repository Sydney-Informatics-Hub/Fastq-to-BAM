#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: this script is run by bqsr_gather_run_parallel.pbs
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

labSampleID=$1

tables=$(find ../BQSR_tables -name "${labSampleID}.*.recal_data.table")
tables=($tables)
tables_list=''

for file in ${tables[@]}
do
        tables_list+=" -I $file"
done

gatk GatherBQSRReports \
        $tables_list \
        -O ../BQSR_tables/${labSampleID}.recal_data.table

