#!/bin/bash
#MSUB -N bin-concoct.sh
#MSUB -A b1042
#MSUB -q genomics
#MSUB -m ae
#MSUB -M jamesgriffin2013@u.northwestern.edu
#MSUB -l nodes=1:ppn=24
#MSUB -l walltime=48:00:00
#MSUB -o /projects/b1042/Wells/Jim/Wells01/binning/binning-concoct.o
#MSUB -e /projects/b1042/Wells/Jim/Wells01/binning/binning-concoct.e

cd /projects/b1042/Wells/Jim/Wells01/
source activate env_concoct

concoct [--coverage_file COVERAGE_FILE]
         [--composition_file COMPOSITION_FILE] [-c CLUSTERS] [-k KMER_LENGTH]
         [-l LENGTH_THRESHOLD] [-r READ_LENGTH]
         [--total_percentage_pca TOTAL_PERCENTAGE_PCA] [-b BASENAME] [-s SEED]
         [-i ITERATIONS] [-e EPSILON] [--no_cov_normalization]
         [--no_total_coverage] [--no_original_data] [-o] [-d] [-v]
metawrap binning -o binning_concoct -t 24 -a assembled/co_assembly_stored/contig.fa --concoct clean_reads/all_reads_1.fastq clean_reads/all_reads_2.fastq
