#!/bin/bash
 
#########################################################
# 
# Platform: NCI Gadi HPC
# Description: reformatting of split and disc files
# Usage: this script is executed by reformat_split_disc_run_parallel.pbs 
# Details:
# 	Convert split and disc files from sam to BAM, sort, index	
#	Not using --write-index feature of samtools/1.10 re known bug with 
#	threaded sort + index. Has been fixed but Gadi install not updated
# Author: Cali Willet
# cali.willet@sydney.edu.au
# Date last modified: 31/07/2020
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

set -e

module load samtools/1.10

sample=`echo $1 | cut -d ',' -f 1` 
NCPUS=`echo $1 | cut -d ',' -f 2`

cd SplitDisc
temp=Temp_${sample}
mkdir -p $temp

in=( ".split" ".disc")

for (( i = 0; i < ${#in[@]}; i ++ ))
do
	samtools sort -@ $NCPUS \
		-o ${temp}/${sample}${in[$i]}.sort.bam \
		${sample}${in[$i]}.sam
	
	samtools index -@ $NCPUS ${temp}/${sample}${in[$i]}.sort.bam
	
	mv ${temp}/* .	

done

rmdir $temp	
