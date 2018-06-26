## Metagenomic Assembly Pipeline and Notes

* [Metagenomic Assembly Pipeline and Notes](#metagenomic-assembly-pipeline-and-notes)  
         * [Contents:](#contents)  
         * [File descriptions:](#file-descriptions)  
         * [Assembly docs:](#assembly-docs)  
         * [Useful shell/msub commands:](#useful-shellmsub-commands)  

### Contents:
* Code for job submission scripts for metagenomic assembly on Northwestern's Quest computing cluster
* Links to documentation for the assemblers used
* Some useful shell commands I didn't know about for small things, maybe useful for other Wells people later


Input Directory Structure:
/projects/b1042/Wells/User/ProjectName (or /projects/b1052/User/ProjectName)
* int_assembled directory is created by assembly_job_writer.py
* assembled is created by cleanup_intermediates.py
* annotations is created by prokka_job_writer.py

```bash
|--raw  
    ├───{Samples}  
        ├───{Sample}*R1*.gzip  
        ├───{Sample}*R2*.gzip  
├───clean_reads  
    ├───{Samples}  
        ├───{Sample}*R1*.gzip  
        ├───{Sample}*R2*.gzip  
├───int_assembled #This gets deleted by cleanup_intermediates.py  
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
├───assembled  
    ├───{Methods}  
        ├───{Sample}_contigs.fa  
├───annotations  
    ├───{Samples}  
```


### File descriptions:
* `idba_assembly_script_writer.py`:
  * Creates a valid msub job submission script for every PE fastq sample in a directory, then submits w/ subprocess
  * Currently this is pretty inflexible. (hard-coded paths abound!)
  * TODO: argparse?, run locally with openSSH/paramiko?
  * Default parameters for this are what GH used: `--mink 20 --maxk 80 --step 10 --min_contig 500 --num_threads 24`

### Assembly docs:
* [idba](https://github.com/loneknightpy/idba)
* idba requires interleaved fasta files (Not allowed: fastq, compressed, or forward/reverse paired end reads)
* Has a built in function **fq2fa** for fixing our reads: `fq2fa --merge --filter FWD REV FA_OUT`

### Useful shell/msub commands:
* `gunzip *.fq.gz` Decompress every .fq file in the current directory. Useful for idba which doesn't currently handle compressed files
* `showq -u NETID`: Shows your q'd jobs
* `du -h --max-depth=0 DIR`: Calculates directory size in GB for DIR
