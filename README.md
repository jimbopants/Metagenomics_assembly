## Metagenomic Assembly Pipeline and Notes

* [Metagenomic Assembly Pipeline and Notes](#metagenomic-assembly-pipeline-and-notes)  
         * [Contents:](#contents)  
         * [Quest Directory Structure:](#input-directory-structure)
         * [File descriptions:](#file-descriptions)  
         * [Assembly docs:](#assembly-docs)  
         * [Useful shell/msub commands:](#useful-shellmsub-commands)  

### Contents:
* This repo contains mostly python scripts for job submission scripts for metagenomic assembly on Northwestern's Quest computing cluster. See the diagram below for an overview of what steps were used.  
* I also included links to the documentation for all of the 3rd party programs I used and a bit of description of what the programs do and what intermediate files are created and if they have any use.  
* Lastly there are some useful shell commands I didn't know about for small things, maybe useful for other Wells people later ¯\_(ツ)_/¯


### Input Directory Structure:
Root = /projects/b1042/Wells/User/ProjectName (or /projects/b1052/User/ProjectName)
* int_assembled directory is created by assembly_job_writer.py
* assembled is created by cleanup_intermediates.py
* annotations is created by prokka_job_writer.py

```bash
├───raw     # Raw *demultiplexed* reads go here
    ├───{Samples}  
        ├───{Sample}*R1*.gzip  
        ├───{Sample}*R2*.gzip  
├───clean_reads  # After trimmomatic to remove adapters and QC
    ├───{Samples}  
        ├───{Sample}*R1*.gzip  
        ├───{Sample}*R2*.gzip  
├───int_assembled #Deleted by cleanup_intermediates.py  
    ├───idba  
        ├───{Samples}  
            ├───{aligned#, other useless stuff}  
            ├───contigs.fa  
            ├───scaffolds.fa  
            ├───log  
    ├───megahit  
        ├───{Samples}  
            ├───contigs.fa  
            ├───{other stuff}  
├───assembled  # Created after cleanup_intermediates.py
    ├───{Methods}  
        ├───{Sample}_contigs.fa  
├───annotations  
    ├───{Samples}  
```

### File descriptions:
* `assembly_script_writer.py`:
  * Creates a valid msub job submission script for every PE fastq sample in a clean_reads directory, then submits w/ subprocess
  * Default parameters for the assembly algorithms are based on Gao Han's/Yubo's settings. You can edit these in the `config.json` file and pass it to the `--params` argument when you run the code on Quest.
  * command line options are:
```bash
usage: assembly_job_writer.py [-h] [--assembler ASSEMBLER] [--raw RAW]
                              [--out OUT] [--params PARAMS] [--check_modules]
                              [--nosub]

optional arguments:
  -h, --help            show this help message and exit
  --assembler ASSEMBLER, -a ASSEMBLER
                        Comma separated assemblers to use
                        ('idba','spades','megahit')
  --raw RAW, -r RAW     Input directory containing subdirectories with
                        demultiplexed paired end reads
  --out OUT, -o OUT     Directory where assembled reads will go i.e.:
                        out/METHOD/SAMPLE
  --params PARAMS, -p PARAMS
                        A file of tab separated key-val parameters for each
                        method
  --check_modules       Runs validate_modules to check module paths
  --nosub               writes shell scripts but doesn't submit jobs
```

### Assembly docs:
* [idba](https://github.com/loneknightpy/idba)
* idba requires interleaved fasta files (Not allowed: fastq, compressed, or forward/reverse paired end reads)
* Has a built in function **fq2fa** for fixing our reads: `fq2fa --merge --filter FWD REV FA_OUT`

### Useful shell/msub commands:
* `gunzip *.fq.gz` Decompress every .fq file in the current directory. Useful for idba which doesn't currently handle compressed files
* `showq -u NETID`: Shows your q'd jobs
* `du -h --max-depth=0 DIR`: Calculates directory size in GB for DIR (useful to see who is killing our quest storage. No storage makes me sad, delete your intermediate files.)
