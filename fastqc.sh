#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: run fastQC and unpack output
# Usage: this script is executed by fastqc_run_parallel.pbs
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

fastq=$1

fastqc --extract $fastq -o FastQC

# Optional: 
#    -a              Specifies a non-default file which contains the list of
#    --adapters      adapter sequences which will be explicity searched against
#                    the library. The file must contain sets of named adapters
#                    in the form name[tab]sequence.  Lines prefixed with a hash
#                    will be ignored.
