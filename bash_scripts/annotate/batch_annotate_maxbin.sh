#!/bin/bash

#MSUB -N annotate_maxbin.sh
#MSUB -A b1042
#MSUB -q genomics
#MSUB -m ae
#MSUB -M jamesgriffin2013@u.northwestern.edu
#MSUB -l nodes=1:ppn=24
#MSUB -l walltime=24:00:00
#MSUB -o /projects/b1042/Wells/Jim/Wells01/logs/annotate_maxbin.o
#MSUB -e /projects/b1042/Wells/Jim/Wells01/logs/annotate_maxbin.e

# PATHS
#METABAT
#BINS="/projects/b1042/Wells/Jim/Wells01/binning/metabat2_bins/final_bins"
#BASE_OUT="/projects/b1042/Wells/Jim/Wells01/binning/metabat2_bins/annotations/"
#EXT=".fa" # metabat
#MAXBIN:
BINS="/projects/b1042/Wells/Jim/Wells01/binning/maxbin_bins/final_bins"
BASE_OUT="/projects/b1042/Wells/Jim/Wells01/binning/maxbin_bins/annotations/"
EXT=".fasta"

mkdir -p $BASE_OUT
source activate env_prokka
cd $BINS

for f in *$EXT; do
  g=${BASE_OUT}$(basename -s .fasta $f);
  echo $g
  prokka --outdir $g --metagenome --cpus=0 $f;
done
