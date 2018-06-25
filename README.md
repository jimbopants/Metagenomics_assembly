## Metagenomic Assembly Pipeline and Notes

### Contents:
* Code for job submission scripts for metagenomic assembly on Northwestern's Quest computing cluster
* Links to documentation for the assemblers used
* Some useful shell commands I didn't know about for small things, maybe useful for other Wells people later


### Why not Snakemake?
Snakemake is a bioinformatics pipeline building tool that uses rules to automatically create directory structures as needed and then run a series of steps. I had planned on making a snakemake pipeline for my metagenomic analysis steps but decided to stick with python scripts whose output is batch job files to save time. 


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
