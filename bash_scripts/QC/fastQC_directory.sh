#!/bin/bash
#MSUB -N fastQC_directory.sh
#MSUB -A b1042
#MSUB -q genomics
#MSUB -m ae
#MSUB -M jamesgriffin2013@u.northwestern.edu
#MSUB -l nodes=1:ppn=24
#MSUB -l walltime=4:00:00
#MSUB -o /projects/b1042/Wells/Jim/Wells01/fastQC.o
#MSUB -e /projects/b1042/Wells/Jim/Wells01/fastQC.e

# Edit the directories and sample wildcard to run on different sequences:
outdir=/projects/b1042/Wells/Jim/Wells01/fastQC/
readsdir=/projects/b1042/Wells/Jim/Wells01/clean_reads

module load fastqc/0.11.5
mkdir $out_dir
cd $reads_dir
# This expands all files that match $reads_dir/C[anything]/[anything].fq
fastqc C*/*.fq -o $out_dir -t 24
