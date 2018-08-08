#!/bin/bash
#MSUB -N Kraken_setup.sh
#MSUB -A b1042
#MSUB -q genomics
#MSUB -m ae
#MSUB -M jamesgriffin2013@u.northwestern.edu
#MSUB -l nodes=1:ppn=24
#MSUB -l walltime=48:00:00
#MSUB -o /projects/b1052/shared/kraken_setup.o
#MSUB -e /projects/b1052/shared/kraken_setup.e

cd /projects/b1052/shared/
mkdir kraken2_DB
kraken2-build --standard --db kraken2_DB --threads 24

# This doesn't work because I'm on a proxy that can't rsync without setting a variable.
