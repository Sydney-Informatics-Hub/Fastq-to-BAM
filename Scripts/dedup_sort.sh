#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: this script is run by dedup_sort_run_parallel.pbs
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

# Does not like being in run_parallel script
samblaster=/scratch/er01/apps/samblaster/0.1.24/samblaster

\rm -rf ../Dedup_sort/${labSampleID}
mkdir ../Dedup_sort/${labSampleID}

bam_in=../Merge_align/${labSampleID}.merged.nameSorted.bam
bam_out=../Dedup_sort/${labSampleID}/${labSampleID}.coordSorted.dedup.bam #private dir for sorting gives better I/O than jobfs or samtools tmp flag
split_out=../SplitDisc/${labSampleID}.split.sam
disc_out=../SplitDisc/${labSampleID}.disc.sam
dup_log=./Logs/Dedup_sort/${labSampleID}.dedup.log
err=./Logs/Dedup_sort_error_capture/${labSampleID}.err

\rm -rf $bam_out
lfs setstripe -c 15 $bam_out

echo bam_in:$bam_in bam_out:$bam_out split_out:$split_out disc_out:$disc_out dup_log:$dup_log err:$err NCPUS:$NCPUS

samtools view -@ 12 -h $bam_in \
        | $samblaster -M -e -d $disc_out -s $split_out 2>$dup_log \
        | samtools sort -@ ${NCPUS} -m 4G -o $bam_out -

mv ${bam_out}* ../Dedup_sort
rmdir ../Dedup_sort/${labSampleID}

if ! samtools quickcheck ../Dedup_sort/${labSampleID}.coordSorted.dedup.bam
then
        printf "Corrupted or missing BAM\n" > $err
fi

test=$(tail -1 ${dup_log} | awk '$0~/^samblaster\: Marked/') # any samblaster version updates need to be checked for last log line output. This is correct for 0.1.24
if [[ ! $test ]]
then
        printf "Error in samblaster log\n" >> ${err}
fi
