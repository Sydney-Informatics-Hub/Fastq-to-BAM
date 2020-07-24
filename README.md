# Fastq-to-BAM
Optimised pipeline to process whole genome sequence data from fastq to BAM on NCI Gadi

24/7/20 - please note this pipeline has not been tested post upload. Edits have been made. I will remove this message once I have run sample data through the scripts. I will also create a proper user guide and --help flags in the near future. For now, please read the headers of each script as they contain guide and resource info
Also to come is a collection of software exectuables on Gadi. Some of the modules used are 'global' but some require the user to install.

<img src="https://user-images.githubusercontent.com/49257820/87630794-4e7ce680-c779-11ea-9b79-ff22fdb379af.png" width="100%" height="100%">

# Description

Paired fastq files were split into smaller files of approximately 500,000 read pairs with fastp v. 0.20.0 (Chen et al. 2018) for parallel alignment. Reads were aligned to hg38 + alt contigs as downloaded by bwakit v. 0.7.17 ‘run-gen-ref’. Alignment was performed with BWA-MEM v. 0.7.17 (Li 2013), using the ‘M’ flag to mark split hits as secondary. Post-processing of the reads was performed with bwa-postalt.js from bwakit v. 0.7.17 to improve mapping quality of alt hits. During alignment, reads mapping to any of the six HLA genes were extracted to fastq format for later analysis. Scattered BAM files were merged into sample-level BAMs with SAMbamba v. 0.7.1 (Tarasov et al. 2015). Duplicate reads were marked with SAMblaster v. 0.1.24 (Faust and Hall 2014), then sorted by genomic coordinate and indexed with SAMtools v. 1.10 (Li et al. 2009). During duplicate read marking, split and discordant reads were extracted to BAM format for later analysis. Base quality score recalibration was performed with GATK v 4.1.2.0 (Van der Auwera et al. 2013). GATK SplitIntervals was used to define 32 evenly-sized genomic intervals over which GATK BaseRecalibrator was run for each sample. The 32 recalibration tables per sample were merged into one table per sample with GATK GatherReports. Base-recalibrated BAM files were produced with GATK ApplyBQSR in parallel over each of the 3,366 contigs in the hg 38 + alt reference genome, plus the unmapped reads. These were merged into a final BAM per sample with SAMbamba merge or GATK GatherBamFiles.

All computation was performed on NCI ‘Gadi’ HPC, running Centos 8, PBS Pro v. 19, on Intel Xeon Cascade Lake 2 x 24 core nodes each with 192 GB ('normal' queue) or 1536 GB ('hugemem' queue) RAM. All stages of the analysis were run in parallel, either massively parallel using the scatter-gather method, or parallel by sample. Parallelisation across the cluster was achieved through either GNU parallel v. 20191022 (Tange 2018) or Open MPI v. 4.0.2 (Graham et al. 2005), along with NCI's utility script nci-parallel v 1.0.0. 

# Download the reference genome

This pipeline is written for alignment of human whole genome sequence data to the hg38 + alt contigs reference. The 'Reference' folder must be located within the main working directory for the cloned repo (or a symlink created), to ensure the relative paths that are used within the scripts are correct. The Reference directory is 73 GB unpacked. 

The complete set of required files for hg38 + alt can be downloaded with the following commands, which should be subitted as a job to Gadi's 'copyq':

wget https://cloudstor.aarnet.edu.au/plus/s/CHeAuEsBkHalDvI/download -O Reference.tar.gz
wget https://cloudstor.aarnet.edu.au/plus/s/0fN5X8tksH2PGso/download -O Reference.tar.gz.md5

Verify the completeness of the download by checking md5sums.

If you would like to use this pipeline with an alternate reference genome/species, you will need to prepare the required files (SAMtools index, GATK index, GATK split intervals for BQSR* and GATK split intervals  x 3200 for variant calling) as well as edit the scripts to the new reference name. *Split intervals for BQSR must be a minimum of 100 Mbp (eg for 3.2 Gbp hg38 + alt, this is 32 intervals). For ApplyBQSR adjust the number of contigs over which to print recalibrated BAMs. You will also need to source the known variants databases to be used for BQSR and VQSR and edit those scripts that use them accordingly. 

# Overview of steps

0) download, unpack and checksum the Reference
0) make <cohort_name>.config
1) split_fastq
2) align
3) merge_align
4) dedup_sort
5) index
6) bqsr_recal
7) bqsr_gather
8) bqsr_apply
9) bqsr_merge

<cohort_name>.config MUST match the following format:
#SampleID       LabSampleID     SeqCentre       Library(default=1)

where the 'SampleID' is the unique identifier enabling one to recognise which fastq belong to the sme sample, and the 'LabSampleID' is the desired name for the output files for that sample, ie the 'in-house ID'. Some metadata is collected in align_make_input.sh - please read the headers and check this is suitable for your data, and edit as required.
For samples with multiple fastq per sample, this format is fine assuming that all the fastq pairs for any one sample have some ID in common. For samples with more than one library prep, some modifications to 'align_make_input.sh' will be required - please view the headers and 'hashed out' example within that script and edit as required.

Each stage has (at least) 3 scripts: a 'make_input.sh', a 'run_parallel.sh' and a '.sh'. The make input script is to be run first, then the run parallel script submits parallel instances of the .sh task script.
Some stages have additional steps, eg collecting logs (split_fastq and dedup_sort), checking outputs (split_fastq, merge_align, dedup_sort) and bqsr_merge_GATK (bqsr_merge_make_bamLists.sh).
Error checking for each job MUST include:
	a) Checking the 'error capture' directory for each job. Should be empty.
	b) Checking the .e PBS logs. All successful tasks should have '...exited with status 0' but not all '...exited with status 0' tasks are successful
	c) Checking that the right number and size of outputs have been created
	d) For tasks with check scripts (split_fastq, merge_align, dedup_sort) verify the TSV output mathces expectations
	
# References

BWA-MEM: Li 2013 https://arxiv.org/abs/1303.3997
BWA-kit: Li 2014 https://github.com/lh3/bwa/tree/master/bwakit
Fastp: Chen et al 2018 https://academic.oup.com/bioinformatics/article/34/17/i884/5093234
GATK4: Van der Auwera et al. 2013 https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/0471250953.bi1110s43
GNU Parallel: Tange 2018 https://doi.org/10.5281/zenodo.1146014
OpenMPI: Graham et al. 2015 https://dl.acm.org/doi/10.1007/11752578_29
SAMbamba: Tarasov et al. 2015 https://academic.oup.com/bioinformatics/article/31/12/2032/214758
SAMblaster: Faust and Hall 2014 https://academic.oup.com/bioinformatics/article/30/17/2503/2748175
SAMtools: Li et al. 2009 https://www.ncbi.nlm.nih.gov/pubmed/19505943

	
