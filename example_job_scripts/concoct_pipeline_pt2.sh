#!/bin/bash
#MSUB -N concoct_preprocess.sh
#MSUB -A b1042
#MSUB -q genomics
#MSUB -m ae
#MSUB -M jamesgriffin2013@u.northwestern.edu
#MSUB -l nodes=1:ppn=24
#MSUB -l walltime=48:00:00
#MSUB -o /projects/b1042/Wells/Jim/Wells01/logs/concoct_preprocess.o
#MSUB -e /projects/b1042/Wells/Jim/Wells01/logs/concoct_preprocess.e

# Written by JG 7/25/18
# This could easily be merged with part 1 for future runs.
# TODO: Apply initial contig cutting to other binning tools before refinement.

cd /projects/b1042/Wells/Jim/Wells01/
coassembly=/projects/b1042/Wells/Jim/Wells01/assembled/co_assembly_stored
concoct=/projects/b1052/shared/CONCOCT-0.4.0
samples=/projects/b1042/Wells/Jim/Wells01/clean_reads
basedir=/projects/b1042/Wells/Jim/Wells01/CONCOCT_input/
maxbins=200

# Source CONCOCT programs
source activate env_concoct
source source_concoct_paths.sh

# Subset coverage to mean for each contig:
cd $basedir
cut -f1,3- concoct_inputtable.tsv > concoct_inputtableR.tsv

# Run Concoct on 150 bp reads with up to a max of $maxbin bins
concoct -c $maxbins -r 150 --coverage_file concoct_inputtableR.tsv --composition_file $coassembly/contig_c10K.fa -b concoct-output/
