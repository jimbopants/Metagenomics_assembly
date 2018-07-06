#!/bin/bash
#MSUB -N C1B_S2.bwa.sh
#MSUB -A b1042
#MSUB -q genomics
#MSUB -m ae
#MSUB -M jamesgriffin2013@u.northwestern.edu
#MSUB -l nodes=1:ppn=24
#MSUB -l walltime=48:00:00
#MSUB -o /projects/b1042/Wells/Jim/Wells01/mapped_reads/idba_bwa/C1B_S2.bwa2018-07-05-21-27.o
#MSUB -e /projects/b1042/Wells/Jim/Wells01/mapped_reads/idba_bwa/C1B_S2.bwa2018-07-05-21-27.e

module load bwa/0.7.12
bwa index /projects/b1042/Wells/Jim/Wells01/assembled/idba/C1B_S2/contig.fa
bwa mem /projects/b1042/Wells/Jim/Wells01/assembled/idba/C1B_S2/contig.fa /projects/b1042/Wells/Jim/Wells01/clean_reads/C1B/C1B_S2_R2_trimmed.fq /projects/b1042/Wells/Jim/Wells01/clean_reads/C1B/C1B_S2_R1_trimmed.fq > /projects/b1042/Wells/Jim/Wells01/mapped_reads/idba_bwa/C1B_S2/aln-pe.sam
