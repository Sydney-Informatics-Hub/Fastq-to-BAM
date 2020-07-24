#!/usr/bin/env perl

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: summarise fastp logs into a TSV file after running 'split fastq'
# Usage: perl split_fastq_collect_fastp_logs.pl <cohort_name>
# Details:
#	If no cohort name specified, default output name is fastp_summary.txt
#	Output will be in ./Logs/Fastp
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


use warnings;
use strict; 

my $out = './Logs/Fastp/fastp_summary.txt'; 
if ($ARGV[0]) {
	my $prefix = $ARGV[0]; 
	chomp $prefix; 
	$out = "./Logs/Fastp/$prefix\_fastp_summary.txt";
	print "Writing output to $out\n";
}
else {
	print "No output prefix supplied. Please re-run with prefix if desired.\nWriting output to default outfile name: $out\n"; 
}

my @logs = (split ' ', `ls ./Logs/Fastp/*log`);
open (O, ">$out") || die "$! $out\n"; 
print O "#Fastq\tR1_reads\tR2_reads\tR1_bases\tR2_bases\tR1_Q30_bases(%)\tR2_Q30_bases(%)\tR1_Q20_bases(%)\tR2_Q20_bases(%)\tDuplication_rate(%)\tInsert_size_peak(bp)\n";


foreach my $log (@logs) {
	my $sample = $log;
	$sample=`basename $sample`; 
	chomp $sample; 
	$sample=~s/\.log$//; 

	#R1:
	my @data = split('\n', `grep -A 4 "Read1 before filtering:" $log | tail -4`);
	my @r1_a = split(' ', $data[0]); #take elem 2
	my @r1_b = split(' ', $data[1]); #take elem 2 
	my @r1_q20 = split(' ', $data[2]); 
	my $r1_q20 = $r1_q20[2];
	$r1_q20=~s/^\d+\(//;
	$r1_q20=~s/\%\)$//;
	my @r1_q30 = split(' ', $data[3]);
	my $r1_q30 = $r1_q30[2];
	$r1_q30=~s/^\d+\(//;
	$r1_q30=~s/\%\)$//;

	#R2:
	@data = split('\n', `grep -A 4 "Read2 before filtering:" $log | tail -4`);
	my @r2_a = split(' ', $data[0]); #take elem 2
	my @r2_b = split(' ', $data[1]); #take elem 2 
	my @r2_q20 = split(' ', $data[2]); 
	my $r2_q20 = $r2_q20[2];
	$r2_q20=~s/^\d+\(//;
	$r2_q20=~s/\%\)$//;
	my @r2_q30 = split(' ', $data[3]);
	my $r2_q30 = $r2_q30[2];
	$r2_q30=~s/^\d+\(//;
	$r2_q30=~s/\%\)$//;


	#Duplication:
	my @dup = split(' ', `grep "Duplication rate" $log`); 
	my $dup = $dup[2];
	$dup=~s/\%$//;
	
	#Insert size:Insert size peak (evaluated by paired-end reads): 269
	my @ins = split(' ', `grep "Insert size peak" $log`); 
	my $ins = $ins[-1];
		
	
	print O "$sample\t$r1_a[2]\t$r2_a[2]\t$r1_b[2]\t$r2_b[2]\t$r1_q30\t$r2_q30\t$r1_q20\t$r2_q20\t$dup\t$ins\n";
	
}close O; 
