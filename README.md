# Metagenomic Assembly Pipeline and Notes

- [Metagenomic Assembly Pipeline and Notes](#metagenomic-assembly-pipeline-and-notes)
  * [Overview:](#overview)
  * [Directory Structure:](#directory-structure)
  * [Installation:](#installation)
  * [Usage:](#usage)
    + [Method Wrappers:](#method-wrappers)
    + [Using the Config file & Command Line Options:](#using-the-config-file---command-line-options)
    + [Command line input:](#command-line-input)
  * [3rd Party Docs:](#3rd-party-docs)
    + [Assembly](#assembly)
    + [Mapping](#mapping)
    + [Annotation](#annotation)
  * [Useful shell/msub/Quest commands:](#useful-shell-msub-quest-commands)

## Overview:
* The code in this repo is designed to quickly write and submit jobs to Northwestern's computing cluster to run 3rd party bioinformatics tools on multiple samples.
* The general workflow is to manage methods and options by editing the config.yaml file and then run the various method wrapper scripts. Paths are automatically handled by the config file to make running the pipeline relatively straight forward.
* I also included links to the documentation for all of the 3rd party programs I used and a few links I found helpful when choosing programs/paramaters.
* Lastly there are some useful shell/msub commands I didn't know about for small things, maybe useful for other Wells people later '¯\\_(ツ)_/¯
* The overall workflow I am using for my MEC project is:
<img src="/workflow.png" width="500">




## Directory Structure:
The tree below is taken directly from what my project folder on quest looks like:
* assembled is created by assemble_reads.py and modified by cleanup_intermediates.py
* annotations is created by annotate_prokka.py
* mapped_reads is created by map_reads.py

```bash
├───reads     # Raw *demultiplexed* reads go here
    ├───{Samples}  
        ├───{Sample}*R1*.gzip  
        ├───{Sample}*R2*.gzip  
├───clean_reads  # After trimmomatic to remove adapters and QC
    ├───{Samples}  
        ├───{Sample}*R1*.gzip  
        ├───{Sample}*R2*.gzip  
├───assembled # subdirectories made by assemble_reads.py
    ├───idba  
        ├───{Samples}  
            ├───contig.fa  
            ├───scaffolds.fa  
            ├───log  
    ├───megahit  
        ├───{Samples}  
            ├───final.contigs.fa  # Renamed to contig.fa by cleanup_int
├───quast_stats_full
├───annotations  
    ├───{Samples}  
├───mapped_reads
    ├───idba_bwa
├───scripts  
```


## Installation:
Clone this repo while logged into Quest:
`Git clone https://github.com/jimbopants/Metagenomics_assembly`

---
## Usage:

### Method Wrappers:
1. **assemble_reads.py**: This script writes (and submits) assembly jobs. Currently IDBA and megahit are implemented, although other methods should be easy to add by following the method templates in the script.
    * Input: looks for sample folders in paths/clean_reads
    * Output: assembled reads in /assembled/{METHOD}/{SAMPLE}

2. **map_reads.py**: This script wraps the mapping programs BWA and bowtie2 as specified by the config file.
    * Input: contig files defined by [paths]["reads_to_map"] + ["contig_names"] in the config file
    * Output: BAM alignment files in /mapped_reads in /assembled/{METHOD}/{SAMPLE}

3. **quast_assembly_stats.py**: Calculates assembly statistics for all of the assemblies to compare performance and QC the assembly step.
    * Input: assembled contigs found in "assembled"/{idba,megahit} based on active assembly methods in the config file.

---
### Using the Config file & Command Line Options:
By default, all of the method wrappers will load options from config.yaml. Parameters are divided into 3 sections:
1. **Run and header options**:
    * Controls msub options used in the job submission scripts (allocation, resources, email, etc.) as well as whether to submit jobs after they are written.
2. **Assembly and mapping methods**:
    * Methods to use for assembly and mapping as well as method specific parameters are stored here.
    * The defaults come from Gao Han and Yubo's work.
    * Commented out methods are ignored.
3. **Paths**:
    * The default directory structure for input/output from various methods is listed here for organization and to save time when scaling analyses for many samples.
    * **base_dir** needs to be updated to wherever your analyses are running.

---
### Command line input:
For debugging or quickly updating parameters without editing the config file, several common paths and options can be overwritten using command line arguments. All of the method wrappers can be called with the following options which will override the config file values:  
```bash
--directory, -d     Directory containing per-sample folders to process.

--method, -m        idba,megahit,bwa,bowtie (comma separated list)

--sample, -s        Path to single sample directory. Forces script to operate on a single sample only.

--out, -o           Output folder

--sub SUB           (Boolean: True/False) (overrides config)
                    Writes shell scripts but doesnt submit jobs
```

---
## 3rd Party Docs:
### Assembly
* [idba](https://github.com/loneknightpy/idba)
    * idba requires interleaved fasta files (Not allowed: fastq, compressed, or forward/reverse paired end reads)
    * Has a built in function **fq2fa** for formatting demultiplexed fastq reads: `fq2fa --merge --filter FWD REV FA_OUT`
* [megahit](https://github.com/voutcn/megahit)
    * In most of our experience so far, this underperforms IDBA & metaspades
* [(meta)-spades](http://cab.spbu.ru/software/spades/)
    * I don't have a working wrapper for this yet. =(
---
### Mapping
* [Mapping comparison by the amazing Meren](http://merenlab.org/2015/06/23/comparing-different-mapping-software/)
* [BWA](https://github.com/lh3/bwa)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
---
### Annotation
* [prokka](https://github.com/tseemann/prokka)
    * Basically a wrapper for prodigal/hmmer with some automatic database construction. Seems easier to use than setting up each component individually
---
### Useful Links
* [Uppsala Metagenomics Workshop rtd](https://metagenomics-workshop.readthedocs.io/en/latest/annotation/quantification.html)
* [Metagenomics for microbial ecologists rtd](https://metagenomic-methods-for-microbial-ecologists.readthedocs.io/en/latest/#)

---
## Useful shell/msub/Quest commands:
* [Quest user guide](https://kb.northwestern.edu/page.php?id=72406)
* `showq -u NETID`: Shows your q'd jobs
* `checkjob JOBID`: Show details about a job
* `module load utilities; grouplist <allocationID>`: Check whether you are part of an allocation
* `gunzip *.fq.gz` Decompress every .fq file in the current directory. Useful for idba which doesn't handle compressed files
* `du -h --max-depth=0 DIR`: Calculates directory size in GB for DIR (useful to see who is killing our quest storage. No storage makes me sad, delete your intermediate files.)
