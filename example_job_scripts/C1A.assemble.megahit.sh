#!/bin/bash
#MSUB -N C1A.assemble.megahit.sh
#MSUB -A b1042
#MSUB -q genomics
#MSUB -m ae
#MSUB -M jamesgriffin2013@u.northwestern.edu
#MSUB -l nodes=1:ppn=24
#MSUB -l walltime=48:00:00
#MSUB -o /projects/b1042/Wells/Jim/Wells01/assembled/job_logs/C1A.assemble.megahit2018-07-02-01-55.o
#MSUB -e /projects/b1042/Wells/Jim/Wells01/assembled/job_logs/C1A.assemble.megahit2018-07-02-01-55.e

module load megahit/1.0.6.1
megahit --k-min 21 --k-max 121 --k-step 10 -t 12 -m 0.9     -1 /projects/b1042/Wells/Jim/Wells01/clean_reads/C1A/C1A_S1_R2_trimmed.fq -2 /projects/b1042/Wells/Jim/Wells01/clean_reads/C1A/C1A_S1_R1_trimmed.fq -o /projects/b1042/Wells/Jim/Wells01/assembled/megahit/C1A