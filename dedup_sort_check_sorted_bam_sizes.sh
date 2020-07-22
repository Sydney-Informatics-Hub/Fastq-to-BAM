#!/bin/bash

# Expected ratio of merged unsorted BAM to sorted BAM is ~ 0.55 (sorted BAM smaller than unsorted)
# Supply cohort name as first positional argument to bash if desired

cohort=$1 

output=./Logs/Dedup_sort/check_unsorted_vs_sorted_bam_output_sizes.txt

if [ "$cohort" ]
then
	output=./Logs/Dedup_sort/${cohort}_check_unsorted_vs_sorted_bam_output_sizes.txt
fi

rm -rf $output

printf "#ID\tAlignMerged\tDedupSort\n" >> $output

awk 'NR>1' samples.config | while read LINE
do 
	labSampleID=`echo $LINE | cut -d ' ' -f 2`
	merged_bam_size=$(ls -lh Align_merged/${labSampleID}.merged.nameSorted.bam | awk '{print $5}' | sed 's/G$//')
        sorted_bam_size=$(ls -lh Dedup_sort/${labSampleID}.coordSorted.dedup.bam | awk '{print $5}' | sed 's/G$//')	
	ratio=$(bc -l <<< "scale=2; ${sorted_bam_size}/${merged_bam_size}")
        printf "${labSampleID}\t${merged_bam_size}\t${sorted_bam_size}\t${ratio}\n" >> $output
done
