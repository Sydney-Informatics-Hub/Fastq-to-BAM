# FASTQ-TO-BAM User Guide

## Description  

This repository contains a genome alignment workflow that takes raw FASTQ files, aligns them to a reference genome and outputs analysis ready BAM files. This workflow is designed for the National Computational Infrastructure's (NCI) Gadi supercompter, leveraging multiple nodes on NCI Gadi to run all stages of the workflow in parallel, either massively parallel using the scatter-gather approach or parallel by sample. Parellelisation across the cluster was achieved with Open MPI v.4.0.2 [(Graham et al. 2005)](https://dl.acm.org/doi/10.1007/11752578_29), along with NCI's utility script nci-parallel v 1.0.0. It consists of a number of stages and follows the BORAD Institute's best practice recommendations. The stages of the pipeline are as follows:

#### Step 1: Split FASTQ files 
Paired FASTQ files are initially split into smaller files of approximately 500 K read pairs with fastp v.0.20.0 [(Chen et al. 2018)](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234) for parallel alignment. Reads are then aligned to Hg38 + alt contigs as downloaded by bwakit v. 0.7.17 ‘run-gen-ref’. 

#### Step 2: Align and merge
Alignment of reads to the reference assembly is then performed by BWA-MEM v.0.7.17 [(Li 2013)](https://arxiv.org/abs/1303.3997), using the ‘M’ flag to mark split hits as secondary. Post-processing of the reads is performed with bwa-postalt.js from bwakit v. 0.7.17 to improve mapping quality of alt hits. Scattered BAM files are merged into sample-level BAMs with SAMbamba v.0.7.1 [(Tarasov et al. 2015)](https://academic.oup.com/bioinformatics/article/31/12/2032/214758). 

#### Step 3: Mark duplicate reads, extract split and discordant reads
Duplicate reads are marked with SAMblaster v. 0.1.24 [(Faust and Hall 2014)](https://academic.oup.com/bioinformatics/article/30/17/2503/2748175), then reads are sorted by genomic coordinate and indexed with SAMtools v. 1.10 [(Li et al. 2009)](https://www.ncbi.nlm.nih.gov/pubmed/19505943). During duplicate read marking, split and discordant reads can be extracted to BAM format for later analysis.  

#### Step 4: Perform and apply BQSR (optional)
Base quality score recalibration is performed with GATK v 4.1.2.0 [(Van der Auwera et al. 2013)](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/0471250953.bi1110s43). GATK SplitIntervals is used to define a number of evenly-sized genomic intervals over which GATK BaseRecalibrator is run for each sample. Eah of the recalibration tables per sample are merged into one table per sample with GATK GatherReports. Base-recalibrated BAM files are produced with GATK ApplyBQSR in parallel over each of the contigs in the reference genome assembly, plus the unmapped reads. These are then merged into a final BAM per sample with SAMbamba merge or GATK GatherBamFiles.

<img src="https://user-images.githubusercontent.com/49257820/87630794-4e7ce680-c779-11ea-9b79-ff22fdb379af.png" width="100%" height="100%">  

## Set up 

To run this pipeline you will need to clone the GitHub repository, prepare your input files including a config file and reference assembly.  

After completing set up, `Fastq-to-BAM` directory structure should resemble the following:

```
|-- Reference
    |-- BQSR_intervals 
    |-- BQSR_apply_intervals
|-- Scripts 
    |-- Inputs  
    |-- Logs  
|-- <cohort>.config
|-- create_project.bash
|-- Fastq
```

### Installation 

Clone the Fastq-to-bam repository: 

```
module load git
git clone https://github.com/Sydney-Informatics-Hub/Fastq-to-BAM.git
```  

### Software   

All software required to run the fastq-to-bam pipeline are globally installed on Gadi: 

* bwa/0.7.17
* fastp/0.20.0
* fastqc/0.11.7
* gatk/4.1.8.1
* java/jdk1.8.0_60
* k8/0.2.5
* nci-parallel/1.0.0
* openmpi/4.0.2
* sambamba/0.7.1
* samblaster/0.1.24
* samtools/1.10

### Inputs

The fastq-to-bam pipeline requires users supply the following inputs: 
* A config file listing all samples  
* Short read sequences in FASTQ format  
* An indexed reference genome assembly in FASTA format  
* A set of high-confidence variants in VCF format (if performing BQSR)  

#### 1. Prepare the config file 

The config file must have the suffix '.config' (i.e. cohortname.config) and one row per sample, matching the format SampleID\tLabSampleID\tSeqCentre\tLibrary(default=1) where: 

   - SampleID is the unique identifier enabling one to recognise which FASTQs belong to the same sample 
   - LabSampleID is the desired name for the output files for that sample i.e. in-house ID 
   - SeqCentre 
   - Library (default=1)  

Save the config file to the base `Fastq-to-BAM` directory. 
     
#### 2. Prepare the reference genome

Provided with this pipeline is the (+ alt contigs) reference assembly. The complete set of Hg38 required files can be downloaded with the following commands, which should be submitted to Gadi's 'copyq':

```
wget https://cloudstor.aarnet.edu.au/plus/s/CHeAuEsBkHalDvI/download -O Reference.tar.gz
wget https://cloudstor.aarnet.edu.au/plus/s/0fN5X8tksH2PGso/download -O Reference.tar.gz.md5
```

Next, unpack the Reference tar file: 

```
tar -xvf Reference.gz.tar
````

The `Reference` directory must be located within the main working directory for the cloned repository (or a symlink created), to ensure the relative paths that are used within the scripts are correct. 
  
If you would like to use this pipeline with a different reference genome, you will need to prepare the required files (SAMtools index, GATK index, GATK split intervals for BQSR and GATK split intervals for variant calling). Split intervals for BQSR must be a minimum of 100 Mbp (e.g. for 3.2 Gbp hg38 + alt, this is 32 intervals). For ApplyBQSR adjust the number of contigs over which to print recalibrated BAMs. You will also need to source the known variants databases to be used for BQSR and VQSR and edit those scripts that use them accordingly. We have provided a script to prepare reference intervals. Use it by editing and running: 
   
```
qsub create_gatk_ref_intervals.pbs
```

#### 3. Run `bash create_project.bash` from the Fastq-to-BAM base directory 

 This script will edit all pipeline scripts and create the required directory set up. Users will be prompted to provide the following:

* NCI project e.g. ab00  
* Name of scripts directory (relative or absolute) e.g. ./Scripts or /scratch/ab00/Fastq-to-BAM/Scripts  
* Basename of config file e.g. cohortname for cohortname.config
* Full path to the reference directory e.g. /scratch/ab00/Fastq-to-BAM/Reference  
* Name of your reference genome sequence (include suffix) e.g. Hg38.fasta   
* Required read/write access to any Gadi storage other than the supplied NCI project e.g. er01 

## Usage 

Users will run a series of scripts for each stage of the fastq-to-bam pipeline. After the Create Project stage, each subsequent stage consists of a `make_input.sh` script which creates a list of input files and a `run_parallel.pbs` script which submits the `<task>.sh` in parallel to the PBS job scheduler. Some stages have additional 'checker' scripts that can be run to confirm the integrity of a process' output. Each stage is detailed below.
   
All scripts must be run from the `Scripts` directory. For how to calculate the resource requirements for each `_run_parallel.pbs` script, please see [Resource Usage](#resource-usage) and [Benchmarking metrics](#benchmarking-metrics) sections. 

### FastQC 

This step will produce quality reports for all fastq files in a list. Each fastq file will be run as a separate task and each task is processed in parallel. For an explanation of reports, see the [FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/). We recommend creating aggregate reports of FastQC results for all FASTQ files using MultiQC. See the [MultiQC documentation](https://multiqc.info/docs/) for more information. 

Steps for running fastqc scripts are: 

**1. Run `bash fastqc_make_input.sh /path/to/fastqs`**

Creates a list of inputs to be submitted as separate tasks. This script assumes your fastq files are stored in a directory called `Fastq` and will output a list of fastq files to the `Inputs` directory. From this list, FASTQ files will be organised by size and run from largest to smallest in the next step to improve overall job efficiency. 

**2. Run `fastqc_run_parallel.pbs`**  

This script runs `fastqc.sh` for each FASTQ file as an independent task.  

### Split FASTQ

In this step paired FASTQ files are split into smaller files of approximately 500 K read pairs for parallel alignment to the reference genome assembly.  

Steps for running the split_fastq scripts are: 

**1.  Run `bash split_fastq_make_input.sh <config prefix>`**

This script will search the `Fastq` directory and search for all FASTQ pairs belonging to the samples listed in the config file. It will then print a list of unique FASTQ identifiers to the .inputs file. This script assumes FASTQ files are paired. For samples with multiple pairs, please do not concatenate them into one massive pair per sample. That is general advice that extends beyond the scope of this analysis. 

**2. Run `qsub split_fastq_run_parallel.pbs`**  

This script will run `split_fastq.sh`, which runs a splitting task for each of the input files created previously. Before running, confirm the fq1 and fq2 regex in `split_fastq.sh` match your FASTQ files and edit as required. 

**3. Run `qsub split_fastq_check_fastq_input_vs_split_output_run_parallel.pbs`**  

Check that the number of output reads in split FASTQ files matches that detected by fastQC in the unsplit FASTQ files. This step assumes fastQC has previously been run, with unzipped results in `./FastQC`. Before running, check the regex in `split_fastq_check_fastq_input_vs_split_output.sh` match your FASTQ files and edit as required.  

### Align

In this step parallel alignment of 500 K FASTQ pairs is performed with bwakit.  

Steps for running the align scripts are: 

**1. Run `bash align_make_inputs.sh <config prefix>`**

This script will report the total number of split FASTQ pairs to be aligned. If > 20 K input files then users will need to manually batch the list of inputs to make approximately evenly sized input files of ~20 K inputs per input file. Name these files `align.input-chunk<N>` and modify `run_parallel.pbs` script accordingly between submissions. Use the number of tasks printed to stout to determine resource requirements. 

**2. Run `qsub align_run_parallel.pbs`**

This script will run `align.sh` or `align_nonhuman.sh`. These scripts run BWA-mem and post-processing over parallel tasks for each chunk of ~20,000 inputs. Before running, confirm the correct `SCRIPT=` is unhashed in `align_run_parallel.pbs` and edit accordingly. 

### Initial merge

In this step all parallel alignment output files will be merged to form a single BAM file per sample. 

Steps for running the merge_align scripts are: 

**1. Run `bash merge_align_make_input.sh <config prefix>`**  

This script will make a list of input files for parallel merging of alignment outputs. This script can be adjusted for cohorts with a binomial difference in coverage (e.g. 30x/60x normal/tumour samples). If a binomial grouping is desired, change `group=false` to `group=true`. 

**2. Run `qsub merge_align_run_parallel.pbs`**  

This script will run `merge_align.sh` which merges the chunked small BAM files into a single BAM file for each sample. 

**3. Run `bash merge_align_check_fastq_input_vs_bam_output_sizes.sh <config prefix>`**  

This script will check that merged alignments meet the expected size. We expect merged BAM files will be ~1.4-1.5x large than FASTQ files. If a preprocessing step has output FASTQ files in a different compression level, this ratio can deviate from expected. 

### Mark duplicates and sort BAM

In this step, duplicate reads are marked, then sorted by genomic coordinate and indexed. During duplicate read marking, split and discordant reads are extracted to BAM format for later analysis. 

Steps for running dedup_sort scripts are: 

**1. Run `bash dedup_sort_make_input.sh <config prefix>`**  

This script reads samples from the config file and creates an inputs list for duplicate read marking and sorting. This script can be adjusted for cohorts with a binomial difference in coverage (e.g. 30x/60x normal/tumour samples). If a binomial grouping is desired, change `group=false` to `group=true`.  

**2. Run `qsub dedup_sort_run_parallel.pbs`**  

This script will run `dedup_sort.sh` which sorts reads by coordinates and mark duplicate reads, for each sample. It will also extract split and discordant reads to separate files. 

### Index BAM

In this step, BAI index files are created for the dedup/sorted BAM files. 

Steps for running the index scripts are: 

**1. Run `bash index_make_input.sh <config prefix>`**  

This script creates an inputs list for indexing of BAM files. This script can be adjusted for cohorts with a binomial difference in coverage (e.g. 30x/60x normal/tumour samples). If a binomial grouping is desired, change `group=false` to `group=true`. 

**2. Run `qsub index_run_parallel.pbs`**

This script will run `index.sh` for all dedup/sorted BAM files in parallel. 

### BQSR 

This step will create base quality score recalubration (BQSR) tables based on high-confidence, population level variant call sets. For more information on this process see GATK's [BaseRecalibrator tool documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-). 

The steps for running the BQSR_recal scripts are:  

**1. Run `bash bqsr_recal_make_input.sh <config prefix>`**  

This script creates an inputs list for parallel execution of GATK's BaseRecalibrator tool. This script can be adjusted for cohorts with a binomial difference in coverage (e.g. 30x/60x normal/tumour samples). If a binomial grouping is desired, change `group=false` to `group=true`. The comma-delimited list of genomic intervals should consist of 100 Mb chunks per sample i.e. for Hg38, this is 32 entries per sample. GATK recommends BQSR operate over min 100 Mb per chunk. This script uses the human reference assmbly files by default, if using a different reference assembly edit the script accordingly. 

**2. Run `qsub bqsr_recal_run_parallel.pbs`**  

This script will run `bqsr_recal.sh` or `bqsr_recal_nonhuman.sh` in parallel for the number of interval chunks for the supplied reference assembly, per sample. Before running, confirm the correct `SCRIPT=` is unhashed in `bqsr_recal_run_parallel.pbs` and edit accordingly. 

### Gather BQSR reports 

This step gathers the BQSR reports into a single file for each sample. 

The steps for running the BQSR_gather scripts are:  

**1. Run `bash bqsr_gather_make_input.sh <config prefix>`**  

This script makes an inputs list for parallel or serial execution of GATK's GatherBQSRReports. It will create a sample list for gathering chunked BQSR recal tables to one table per sample. 

**2. Run `qsub bqsr_gather_run_parallel.pbs` OR `qsub bqsr_gather_serial.pbs`** 

Running `bqsr_gather_run_parallel.pbs` will run `bqsr_gather.sh` for each BAM in parallel. Because of resource requirements for this process, this script is best suited to larger bathces of samples. BQSR BaseRecalibrator produces one table file per BQSR task (number of samples * number of 100 Mbp intervals). These tables need to be gathered into a merged table for each sample. 

Running `bqsr_gather_serial.pbs` will run `bqsr_gather.sh` for each BAM serially. Because of resource requirements for this process, this script is best suited to small bathces of samples. BQSR BaseRecalibrator produces one table file per BQSR task (number of samples * number of 100 Mbp intervals). These tables need to be gathered into a merged table for each sample.  

### Apply BQSR

This process runs GATK's ApplyBQSR tool over each of the contigs for each sample in parallel to adjust base quality scores.   

The steps for running the BQSR_apply scripts are:  

**1. Run `bash bqsr_apply_make_input.sh <config prefix>`**  

This script creates an inputs list for applying BQSR reports to each BAM file. This script can be adjusted for cohorts with a binomial difference in coverage (e.g. 30x/60x normal/tumour samples). If a binomial grouping is desired, change `group=false` to `group=true`. This script uses the human reference assmbly files by default, if using a different reference assembly edit the script accordingly. 

**2. Run `qsub bqsr_apply_run_parallel.pbs`**  

This script will run `bqsr_apply.sh` for each contig for each sample in parallel. 

### Final merge

This is the final step which merges recalibrated split BAM files into a final BAM file for each sample using GATK's GatherBamFiles tool. 

**1. Run `bash bqsr_merge_make_input.sh <config prefix>`**  

This script creates an inputs list of split BAMs to be merged for each sample. This script can be adjusted for cohorts with a binomial difference in coverage (e.g. 30x/60x normal/tumour samples). If a binomial grouping is desired, change `group=false` to `group=true`. This script uses the human reference assmbly files by default, if using a different reference assembly edit the script accordingly. 

**2. Run `qsub bqsr_merge_GATK_run_parallel.pbs`**  

This script will run `bqsr_merge_GATK.sh` for each sample in parallel. 

## Resource usage

These resource request recommendations are based on benchmarking performed using [6x Platinum Genomes samples](#benchmarking-metrics). All computation was performed on NCI ‘Gadi’ HPC, running Centos 8, PBS Pro v. 19, on Intel Xeon Cascade Lake 2 x 24 core nodes each with 192 GB ('normal' queue), 128/256 GB ('normalbw' queue) or 1536 GB ('hugemem' queue) RAM. 

Compute performance is dependant on genome size, composition, and complexity, therefore we recommend those working with non-human datasets perform one round of benchmarking before running the pipeline on their full dataset.  


| Job          | Queue    | Job distribution           | NCPU/task | Mem/task   | Walltime | Notes                       |
|--------------|----------|----------------------------|-----------|------------|----------|-----------------------------|
| FastQC       | Normalbw | 1 task/FASTQ               | 1         | 4-5 GB     | ~30min   |                             |
| Split_fastq  | Normal   | 1 task/FASTQ pair          | 4         | see note   | ~40min   | Request whole nodes only    |
| Align        | Normal   | 1 task/Split FASTQ         | 2         | see note   | ~30min   | Request 120nodes/1hr for 20k chunk|
| Merge_align  | Normalbw | 1 task/Sample              | 28        | 128 GB     | ~1hr     | Batch samples by size if binomial |
| Dedup_sort   | Normalbw | 1 task/Sample              | 28        | 256 GB     | ~2hr     | Batch samples by size if binomial |
| Index        | Normal   | 1 task/Sample              | 24        | half node  | 3-6.5min |                             |
| BQSR_recal   | Normal   | BQSR intervals/sample      | 2         | see note   | 10-20min |                             |
| BQSR_gather  | Normal   | 1 task/Sample              | 1         | 4 GB       | 10min    |                             |
| BQSR_apply   | Normal   | Contigs/sample             | 1.5-2.5 node/sample  | see NCPU/task | 1.5hr    | Avg time to print a contig is 1min. max is chr1 (40-60min) |
| BQSR_merge   | Normal   | 1 task/Sample              | 3         | 12 GB      | 1hr    | 1 CPU/task on hugemem       |

  
## Benchmarking metrics 

**6x Platinum Genomes human samples (~100x coverage)**

| #JobName          | CPUs_requested | CPUs_used | Mem_requested(GB) | Mem_used(GB) | CPUtime   | CPUtime_mins | Walltime_req | Walltime_used | Walltime_mins | JobFS_req  | JobFS_used | Efficiency | Service_units(CPU_hours) |
|-------------------|----------------|-----------|-------------------|--------------|-----------|--------------|--------------|---------------|---------------|------------|------------|------------|--------------------------|
| fastqc            | 28             | 28        | 128.0GB           | 128.0GB      | 17:56:01  | 1076.02      | 2:00:00      | 0:45:20       | 45.33         | 100.0MB    | 8.11MB     | 0.85       | 26.44                    |
| split_fastq       | 96             | 96        | 380               | 198.2        | 29:19:58  | 1759.97      | 1:00:00      | 0:22:48       | 22.8          | 209715200b | 16714kb    | 0.80       | 72.96                    |
| check_split_fastq | 96             | 96        | 380               | 189.87       | 6:27:50   | 387.83       | 0:45:00      | 0:05:36       | 5.6           | 209715200b | 16714kb    | 0.72       | 17.92                    |
| align_5           | 5760           | 5760      | 22.27TB           | 18.47TB      | 968:51:47 | 58131.78     | 2:00:00      | 0:12:49       | 12.82         | 11.72GB    | 9.02MB     | 0.79       | 2460.8                   |
| merge_align       | 168            | 168       | 768.0GB           | 687.34GB     | 57:51:12  | 3471.2       | 1:30:00      | 0:55:31       | 55.52         | 600.0MB    | 8.11MB     | 0.37       | 194.31                   |
| dedup_sort_3      | 168            | 168       | 1.5TB             | 1.3TB        | 76:09:56  | 4569.93      | 6:00:00      | 2:06:53       | 126.88        | 600.0MB    | 8.11MB     | 0.21       | 451.14                   |
| index             | 144            | 144       | 570.0GB           | 281.75GB     | 4:48:21   | 288.35       | 0:20:00      | 0:05:49       | 5.82          | 300.0MB    | 8.16MB     | 0.34       | 27.92                    |
| bqsr_recal        | 384            | 384       | 1.48TB            | 1.37TB       | 43:45:07  | 2625.12      | 0:40:00      | 0:16:28       | 16.47         | 800.0MB    | 8.18MB     | 0.42       | 210.77                   |
| bqsr_gather       | 6              | 6         | 24.0GB            | 22.44GB      | 0:04:06   | 4.1          | 0:10:00      | 0:00:49       | 0.82          | 100.0MB    | 150.0B     | 0.83       | 0.16                     |
| bqsr_apply        | 720            | 720       | 2.78TB            | 1.08TB       | 137:59:57 | 8279.95      | 1:30:00      | 0:27:58       | 27.97         | 1.46GB     | 9.17MB     | 0.41       | 671.2                    |
| bqsr_merge_2      | 6              | 6         | 72.0GB            | 72.0GB       | 4:32:52   | 272.87       | 2:00:00      | 0:55:59       | 55.98         | 100.0MB    | 174.0B     | 0.81       | 33.59                    |

**16x Canine samples**

| #JobName        | CPUs_requested | CPUs_used | Mem_requested | Mem_used | CPUtime    | CPUtime_mins | Walltime_req | Walltime_used | Walltime_mins | JobFS_req | JobFS_used | Efficiency | Service_units(CPU_hours) |
|-----------------|----------------|-----------|---------------|----------|------------|--------------|--------------|---------------|---------------|-----------|------------|------------|--------------------------|
| fastqc   | 28             | 28        | 128.0GB       | 128.0GB  | 10:42:42   | 642.7        | 2:00:00      | 0:30:44       | 30.73         | 100.0MB   | 8.11MB     | 0.75       | 17.93                    |
| split_fastq   | 96             | 96        | 200.0GB       | 197.53GB | 18:44:06   | 1124.1       | 0:40:00      | 0:27:07       | 27.12         | 200.0MB   | 8.16MB     | 0.43       | 86.77                    |
| align | 5760           | 5760      | 22.27TB       | 14.53TB  | 524:08:31  | 31448.52     | 1:00:00      | 0:15:18       | 15.3          | 11.72GB   | 8.74MB     | 0.36       | 2937.6                   |
| merge_align   | 448            | 448       | 2.0TB         | 1.28TB   | 31:43:02   | 1903.03      | 2:00:00      | 0:12:34       | 12.57         | 1.56GB    | 8.11MB     | 0.34       | 117.29                   |
| index         | 384            | 384       | 1.48TB        | 467.44GB | 2:51:33    | 171.55       | 0:20:00      | 0:02:47       | 2.78          | 800.0MB   | 8.16MB     | 0.16       | 35.63                    |
| dedup_sort    | 448            | 448       | 4.0TB         | 2.88TB   | 37:50:32   | 2270.53      | 1:00:00      | 0:25:31       | 25.52         | 1.56GB    | 8.11MB     | 0.2        | 241.94                   |
| bqsr_recal    | 768            | 768       | 2.97TB        | 2.23TB   | 25:25:12   | 1525.2       | 0:40:00      | 0:05:26       | 5.43          | 1.56GB    | 8.19MB     | 0.37       | 139.09                   |
| bqsr_gather   | 16             | 16        | 64.0GB        | 16.2GB   | 0:05:28    | 5.47         | 0:10:00      | 0:00:40       | 0.67          | 100.0MB   | 388.0B     | 0.51       | 0.36                     |
| bqsr_apply  | 1920           | 1920      | 7.42TB        | 51.51GB  | 1:09:03    | 69.05        | 3:00:00      | 0:00:13       | 0.22          | 800.0GB   | 115.18KB   | 0.16       | 13.87                    |
| bqsr_merge    | 48             | 48        | 190.0GB       | 130.35GB | 3:13:45    | 193.75       | 2:00:00      | 0:19:44       | 19.73         | 100.0MB   | 8.16MB     | 0.2        | 31.57                    |                                        |                                        |

## Additional notes
  
* Users with non-human datasets will need to create the required reference interval files using `create_gatk_ref_intervals.pbs` and edit the align_run_parallel.pbs and bqsr_recal_run_parallel.pbs scripts to run the non-human version of these steps.  
* Given the impact genome complexity and size can have on compute resources requirements, we recommend users with non-human datasets perform a round of benchmarking using their own data before running this pipeline.  
* Users can adjust scripts for cohorts with a binomial difference in coverage (e.g. 30x/60x normal/tumour samples). If a binomial grouping is desired, change group=false to group=true in scripts.  

## Acknowledgements
   
Acknowledgements (and co-authorship, where appropriate) are an important way for us to demonstrate the value we bring to your research. Your research outcomes are vital for ongoing funding of the Sydney Informatics Hub and national compute facilities.  
   
### Authors 

- Cali Willet (Sydney Informatics Hub, University of Sydney)  
- Tracy Chew (Sydney Informatics Hub, University of Sydney)
- Georgina Samaha (Sydney Informatics Hub, University of Sydney)
- Rosemarie Sadsad (Sydney Informatics Hub, University of Sydney)
- Andrey Bliznyuk (National Computational Infrastructure)
- Roger Edberg (National Computational Infrastructure)

### Suggested acknowledgements

The authors acknowledge the technical assistance provided by the Sydney Informatics Hub, a Core Research Facility of the University of Sydney. This research/project was undertaken with the assistance of resources and services from the National Computational Infrastructure (NCI), which is supported by the Australian Government and Australian BioCommons which is enabled by NCRIS via Bioplatforms Australia funding.

## Cite us to support us! 
   
If you use our pipeline, please cite us:

Sydney Informatics Hub, Core Research Facilities, University of Sydney, 2021, The Sydney Informatics Hub Bioinformatics Repository, \<date accessed\>, https://github.com/Sydney-Informatics-Hub/Bioinformatics   

## References 

Chen, S., Zhou Y., Chen Y., Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 01 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560  

Faust G.G., Hall I.M., SAMBLASTER: fast duplicate marking and structural variant read extraction, Bioinformatics, Volume 30, Issue 17, 1 September 2014, Pages 2503–2505, https://doi.org/10.1093/bioinformatics/btu314  

Graham, R.L., Woodall, T.S., Squyres, J.M. 2005. Open MPI: a flexible high performance MPI. In Proceedings of the 6th international conference on Parallel Processing and Applied Mathematics (PPAM'05). Springer-Verlag, Berlin, Heidelberg, 228–239. DOI:https://doi.org/10.1007/11752578_29  

Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R., & 1000 Genome Project Data Processing Subgroup. 2009. The Sequence Alignment/Map format and SAMtools. Bioinformatics (Oxford, England), 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352  

Li, H. 2013. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997  

Tarasov A., Vilella A.J., Cuppen E., Nijman I.J., Prins P., Sambamba: fast processing of NGS alignment formats, Bioinformatics, Volume 31, Issue 12, 15 June 2015, Pages 2032–2034, https://doi.org/10.1093/bioinformatics/btv098  

Van der Auwera G.A., Carneiro M.O., Hartl C., Poplin R., del Angel G., Levy-Moonshine A., Jordan T., Shakir K., Roazen D., Thibault J., Banks E., Garimella K.V., Altshuler D., Gabriel S. and DePristo M.A. (2013), From FastQ Data to High-Confidence Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline. Current Protocols in Bioinformatics, 43: 11.10.1-11.10.33. https://doi.org/10.1002/0471250953.bi1110s43  
