#!/bin/bash
# Minpath wrapper bash script
# Written by JG 8/6/18

#Command line arguments:
#    -d,--directory      Directory containing bins.
#    -o,--output         Output filepath for pathway predictions
#    -j,--jobname        Jobname for submission and errorlogs

# USAGE:
# msub SCRIPT.sh -d=my_directory -o=output

# Notes:
# This script will iterate through annotation folders created by prokka,
# and then run minpath for pathway prediction on the enzyme file.
# It's based on the steps in this tutorial:
# https://metagenomics-workshop.readthedocs.io/en/latest/annotation/functional_annotation.html


################## msub options ##################
#MSUB -N minpath_pathwys_job.sh
#MSUB -A b1042
#MSUB -q genomics
#MSUB -m ae
#MSUB -M jamesgriffin2013@u.northwestern.edu
#MSUB -l nodes=1:ppn=24
#MSUB -l walltime=48:00:00
#MSUB -o /projects/b1042/Wells/Jim/Wells01/logs/pathways_minpath_1.o
#MSUB -e /projects/b1042/Wells/Jim/Wells01/logs/pathways_minpath_1.e

# Paths:
module load python/anaconda
MINPY="/projects/b1052/shared/MinPath/MinPath1.4.py"
COGMAP="/projects/b1052/shared/MinPath/data/ec2path"
KOMAP="/projects/b1052/shared/MinPath/data/ko"

# Argument Parsing:
DIRECTORY = $1
cd $DIRECTORY

# Run Minpath on the enzyme annotations:
for f in *; do
    if [ -d ${f} ]; then
      python $MINPY -any $f/ECs.txt -map $COGMAP -report $f/ec.report.txt -details $f/ec.detailed.txt
      python $MINPY -any $f/ECs.txt -map $KOMAP -report $f/ko.report.txt -details $f/ec.detailed.txt
    fi
done
