#!/bin/bash

# After splitting, check that the number of output reads in split fastq files
# matches that detected by fastqc in unsplit input
# Update <project> and <cohort>


#1 - Get original number of pairs from fastQC data.txt
#2 - Get number of pairs processed from fastp logs
#3 - Get number of pairs output from line-counting the split fastq
#4 - print warning if mismatch detected

#Gather config:
files=$(ls -1 /scratch/<project>/<cohort>/Fastq/*gz |  sed 's/_R[1-2]\.fastq.gz//'  | cut -d '/' -f 6 | sort | uniq)
files=($files)

module load parallel

fastq=`echo $1 | cut -d ',' -f 1`
sample=`echo $1 | cut -d ',' -f 2`

for fastq in ${files[@]}
do

#1 - Get original number of pairs from fastQC data.txt
qc1=/scratch/<project>/<cohort>/FastQC/${fastq}_R1_fastqc/fastqc_data.txt
qc2=/scratch/<project>/<cohort>/FastQC/${fastq}_R2_fastqc/fastqc_data.txt

fastqc_1=$(grep "Total Sequences" $qc1 | awk '{print $3}')
fastqc_2=$(grep "Total Sequences" $qc2 | awk '{print $3}')

if [[ $fastqc_1 -ne $fastqc_2 ]]
then 
	printf "$sample error: fastqc R1 and R2 read counts do not match - $fastqc_1 and $fastqc_2 respectively\n"
fi

#2 - Get number of pairs processed from fastp logs
fastp_log=/scratch/<project>/<cohort>/Logs/Fastp/${fastq}.log

fastp_1=$(grep -A 1 "Read1 before filtering" $fastp_log | tail -1 | awk '{print $3}')
fastp_2=$(grep -A 1 "Read2 before filtering" $fastp_log | tail -1 | awk '{print $3}')
fastp_tot=$(grep -m 1 "reads passed filter" $fastp_log | awk '{print $4}')
fastp_sum=$(expr $fastp_1 + $fastp_2)

if [[ $fastp_sum -ne $fastp_tot ]]
then 
	printf "$sample error: fastp R1 and R2 input sum does not match fastp total reads read - $fastp_sum and $fastp_tot respectively\n"
fi

if [[ $fastp_1 -ne $fastp_2 ]]
then 
	printf "$sample error: fastp R1 and R2 read counts do not match - $fastp_1 and $fastp_2 respectively\n"
fi


#3 - Get number of pairs output from line-counting the split fastq
splits=/scratch/<project>/<cohort>/Fastq_split

split_1_lines=$(ls ${splits}/*${fastq}_R1.fastq.gz | parallel --will-cite "zcat {} | wc -l " | awk '{s+=$1} END {print s}')
split_2_lines=$(ls ${splits}/*${fastq}_R2.fastq.gz | parallel --will-cite "zcat {} | wc -l " | awk '{s+=$1} END {print s}')

if [[ $split_1_lines -ne $split_2_lines ]]
then 
	printf "$sample error: fastp split R1 and R2 line counts do not match - $split_1_lines and $split_1_lines respectively\n"
fi

split_pairs=$(expr $split_1_lines \/ 4)

#Check 3 sources
if [[ $fastqc_1 -ne $fastp_1 ]]
then 
	printf "$sample error: fastQC and fastp pair counts do not match - $fastqc_1 and $fastp_1 respectively\n"
fi

if [[ $fastp_1 -ne $split_pairs ]]
then 
	printf "$sample error: fastp input and fastp split output pair counts do not match - $fastp_1 and $split_pairs respectively\n"
fi

done
