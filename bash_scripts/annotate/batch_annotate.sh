#!/bin/bash

#MSUB -N annotate_BINS.sh
#MSUB -A b1042
#MSUB -q genomics
#MSUB -m ae
#MSUB -M jamesgriffin2013@u.northwestern.edu
#MSUB -l nodes=1:ppn=24
#MSUB -l walltime=48:00:00
#MSUB -o /projects/b1042/Wells/Jim/Wells01/logs/annotate_BINS.o
#MSUB -e /projects/b1042/Wells/Jim/Wells01/logs/annotate_BINS.e

# PATHS
#maxbin_bins = "/projects/b1042/Wells/Jim/Wells01/binning/maxbin_bins/final_bins"
BINS="/projects/b1042/Wells/Jim/Wells01/binning/metabat2_bins/final_bins"
BASE_OUT="/projects/b1042/Wells/Jim/Wells01/binning/metabat2_bins/annotations/"
EXT=".fa" # change when running maxbin

mkdir -p $BASE_DIR
source activate env_prokka
cd $BINS

for f in *$EXT; do
  g=${BASE_OUT}$(basename -s .fasta $f);
  echo $g
  mkdir -p $g;
  prokka --outdir $g --metagenome --cpus=0 $f;
done
