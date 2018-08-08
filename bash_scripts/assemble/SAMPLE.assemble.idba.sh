#!/bin/bash
#MSUB -N C3B.assemble.idba.sh
#MSUB -A b1042
#MSUB -q genomics
#MSUB -m ae
#MSUB -M jamesgriffin2013@u.northwestern.edu
#MSUB -l nodes=1:ppn=24
#MSUB -l walltime=48:00:00
#MSUB -o /projects/b1042/Wells/Jim/Wells01/assembled/job_logs/C3B.assemble.idba2018-07-04-15-11.o
#MSUB -e /projects/b1042/Wells/Jim/Wells01/assembled/job_logs/C3B.assemble.idba2018-07-04-15-11.e

module load idba/2016_12_longread

idba idba_ud -r /projects/b1042/Wells/Jim/Wells01/fa_interleaved/C3B_trimmed.fa --mink 21 --maxk 81 --step 10 --min_contig 500 --num_threads 24 -o /projects/b1042/Wells/Jim/Wells01/assembled/idba/C3B/
