#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: this script is run by align_run_parallel.pbs
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

set -e

# align.sh for non-human samples

ref=
outdir=../Align_split
logdir=Logs/Align
errdir=Logs/Align_error_capture

fqpair=`echo $1 | cut -d ',' -f 1`
fq1=$(ls ${fqpair}*R1.f*q.gz) #Must check regex for each batch
fq2=$(ls ${fqpair}*R2.f*q.gz)
sampleID=`echo $1 | cut -d ',' -f 2`
centre=`echo $1 | cut -d ',' -f 3`
lib=`echo $1 | cut -d ',' -f 4`
platform=`echo $1 | cut -d ',' -f 5`
flowcell=`echo $1 | cut -d ',' -f 6`
lane=`echo $1 | cut -d ',' -f 7`

outPrefix=$(basename $fqpair)
err=${errdir}/${outPrefix}.err

#bwakit emits default sort order (queryname, but no SO tag in headers) or option to sort
#by coordinate with samtools but does not allow the -n flag to specify sort order by name.
#Sambamba requires queryname sorted with SO tag. Below code does it by force.

echo fqpair:$fqpair fq1:$fq1 fq2:$fq2 sampleID:$sampleID centre:$centre lib:$lib platform:$platform flowcell:$flowcell lane:$lane outPrefix:$outPrefix err:$err ref:$ref NCPUS:$NCPUS

seqtk mergepe $fq1 $fq2 | bwa mem -p -t $NCPUS \
        -R "@RG\tID:${flowcell}.${lane}_${sampleID}_${lib}\tPL:${platform}\tPU:${flowcell}.${lane}\tSM:${sampleID}\tLB:${sampleID}_${lib}\tCN:${centre}" \
        -M ${ref} - 2> ${logdir}/${outprefix}.log.bwamem \
        | samtools sort -n -@ $NCPUS -o ${outdir}/${outPrefix}.aln.bam

if ! samtools quickcheck ${outdir}/${outPrefix}.aln.bam
then
        printf "Corrupted or missing BAM\n" > $err
fi
