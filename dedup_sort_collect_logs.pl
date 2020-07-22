#!/usr/bin/perl

# Run as perl dedup_sort_collect_logs.pl <cohort> 
# If no cohort name specified, default output name is dedup_summary.txt

use warnings;
use strict; 

my $out = './Logs/Dedup_sort/dedup_summary.txt'; 
if ($ARGV[0]) {
	my $prefix = $ARGV[0]; 
	chomp $prefix; 
	$out = "./Logs/Dedup_sort/$prefix\_dedup_summary.txt";
	print "Writing output to $out\n";
}
else {
	print "No output prefix supplied. Please re-run with prefix if desired.\nWriting output to default outfile name: $out\n"; 
}

my @logs = (split ' ', `ls ./Logs/Dedup_sort/*log`);
open (O, ">$out") || die "$! dedup_log.txt\n"; 
print O "#Sample\tTotalReads\tDupReads\tDupPercent\tDiscordantReads\tSplitReads\n"; 

foreach my $log (@logs) {
	my $sample = $log;
	$sample=`basename $sample`; 
	chomp $sample; 
	$sample=~s/\.dedup\.log$//; 
	
	my $reads = `grep Marked $log`;
	chomp $reads; 
	my @reads = split(' ', $reads); #not paticularly robust - check the output each run, code a better method if it ever fails
	my $tot = $reads[4];
	my $dup = $reads[2];
	my $percent = $reads[5];
	$percent=~s/^\(//;
	$percent=~s/\%\)$//;
	
	my $disc = `grep discordant $log`; 
	chomp $disc;
	my @disc = split(' ', $disc); 
	
	my $split = `grep "split reads to" $log`; 
	chomp $split;
	my @split= split(' ', $split); 
	
	print  O "$sample\t$tot\t$dup\t$percent\t$disc[2]\t$split[2]\n"; 	
}close O; 
