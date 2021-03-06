#!/bin/bash
#MSUB -N metawrap-binning-maxbin.sh
#MSUB -A b1042
#MSUB -q genomics
#MSUB -m ae
#MSUB -M jamesgriffin2013@u.northwestern.edu
#MSUB -l nodes=1:ppn=24
#MSUB -l walltime=48:00:00
#MSUB -o /projects/b1042/Wells/Jim/Wells01/binning/metawrap-binning-maxbin.o
#MSUB -e /projects/b1042/Wells/Jim/Wells01/binning/metawrap-binning-maxbin.e

cd /projects/b1042/Wells/Jim/Wells01/
source activate metawrap-env
metawrap binning -o binning_maxbin -t 24 -a assembled/co_assembly_stored/contig.fa --maxbin2 clean_reads/all_reads_1.fastq clean_reads/all_reads_2.fastq
