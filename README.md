## Metagenomic Assembly Pipeline and Notes

### Contents:
* Code for job submission scripts for metagenomic assembly on Northwestern's Quest computing cluster
* Links to documentation for the assemblers used
* Some useful shell commands I didn't know about for small things, maybe useful for other Wells people later

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
