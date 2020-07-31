#!/usr/bin/env perl

#########################################################
# 
# Platform: NCI Gadi HPC
# Description: Cat the HLA from split BWakit to one fastq per sample per gene
# Usage: this script is executed by cat_HLA_run_parallel.pbs 
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

use warnings;
use strict; 

my ($sample,$lab_id) = split(',', $ARGV[0]); 
chomp ($sample,$lab_id); 

my @HLA = ('A','B','C','DQA1','DQB1','DRB1');

my $dir = "./HLA_fastq"; #  check the regex at find depending on whether sample names have shorter non-unique matches

my $out_dir = './HLA_fastq_cat/';


foreach my $gene (@HLA) {
	my $out = "$out_dir\/$lab_id\_HLA\-$gene\.fastq";
	if (-e $out) {
    		`\rm -rf $out`;
	}
	my $split_hla = `find $dir -name "*$sample\*HLA\-$gene\.fq" -print`; #use this regex if all samples in same dir
	$split_hla =~ s/\s+/ /g; 
	#`cat $split_hla > $out`; #Can't exec "/bin/sh": Argument list too long .
	my @split_hla = split(' ', $split_hla);
	my $c = @split_hla; 
	for my $sh (@split_hla) {
		`cat $sh >> $out`; 
	}
}
