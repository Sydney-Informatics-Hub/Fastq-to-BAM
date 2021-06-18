##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: bash fastqc_make_input.sh </path/to/fastqdir>
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

if [ -z "$1" ]
then
        echo "Please provide the path to the directory with your FASTQ files, e.g. sh fastqc_make_input.sh ../Fastq"
        exit
fi

# INPUTS
Fastq=$1
outdir=../FastQC
logdir=./Logs/FastQC
INPUTS=./Inputs
unsorted=${INPUTS}/fastqc_unsorted.inputs
input_file=${INPUTS}/fastqc.inputs

mkdir -p ${INPUTS}
mkdir -p ${logdir}
mkdir -p ${outdir}

rm -rf ${input_file}

for fastq in ${Fastq}/*f*q.gz; do
        basename=$(basename "$fastq" | cut -d. -f1)
        out=${outdir}/${basename}
        logfile=${logdir}/${basename}.oe
        bytes=$(ls -s "${fastq}" | awk '{print $1}')
        mkdir -p ${out}
        echo "${fastq},${out},${logfile},${bytes}" >> ${unsorted}
done

# Reverse numeric sort bytes, comma delimited unsorted file
sort -rnk 4 -t ',' ${unsorted} > ${input_file}
rm -rf ${unsorted}
num_tasks=`wc -l ${input_file} | cut -d' ' -f 1`
echo "$(date): Number of FASTQ files (tasks): $num_tasks"

