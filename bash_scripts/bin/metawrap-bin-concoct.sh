#!/bin/bash
#MSUB -N checkm_test.sh
#MSUB -A b1042
#MSUB -q genomics
#MSUB -m ae
#MSUB -M jamesgriffin2013@u.northwestern.edu
#MSUB -l nodes=1:ppn=24
#MSUB -l walltime=2:00:00
#MSUB -o /projects/b1042/Wells/Jim/Wells01/checkm_test.o
#MSUB -e /projects/b1042/Wells/Jim/Wells01/checkm_test.e

cd /projects/b1042/Wells/Jim/Wells01/
source activate metawrap-env
checkm test ./checkm_test_results
