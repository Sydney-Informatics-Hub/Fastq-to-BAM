#!/bin/bash

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: sort by coordiante with SAMtools and mark 
# duplicates with SAMblaster
# Usage: this script is executed by bqsr_recal_run_parallel.pbs
# Details:
#	A private temp dir is created per sample for sorting,
#	on /scratch. Testing found significant improvements
#	in performance when using scratch compared to jobfs,
#	and for using a private sort dir compared to relying
#	on the SAMtools '-T prefix' flag
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


module load samtools/1.10
module load samblaster/0.1.24

labSampleID=`echo $1 | cut -d ',' -f 1`

cd Dedup_sort

\rm -rf Dedup_sort_${labSampleID}
mkdir Dedup_sort_${labSampleID}

bam_in=../Align_merged/${labSampleID}.merged.nameSorted.bam
bam_out=./Dedup_sort_${labSampleID}/${labSampleID}.coordSorted.dedup.bam #private dir for sorting gives better I/O than jobfs or samtools tmp flag
split_out=../SplitDisc/${labSampleID}.split.sam
disc_out=../SplitDisc/${labSampleID}.disc.sam
dup_log=${labSampleID}.dedup.log
err=../Error_capture/Dedup_sort/${labSampleID}.err

\rm -rf $bam_out
lfs setstripe -c 15 $bam_out


samtools view -@ 12 -h $bam_in \
	| $samblaster -M -e -d $disc_out -s $split_out 2>$dup_log \
	| samtools sort -@ 24 -m 4G -o $bam_out -


mv ${bam_out}* . 
rmdir ./Dedup_sort_${labSampleID}

if ! samtools quickcheck ${labSampleID}.coordSorted.dedup.bam
then 
        printf "Corrupted or missing BAM\n" > $err  
fi

test=$(tail -1 ${dup_log} | awk '$0~/^samblaster\: Marked/') # any samblaster version updates need to be checked for last log line output. This is correct for 0.1.24
if [[ ! $test ]]
then 
        printf "Error in samblaster log\n" >> ${err}
fi

