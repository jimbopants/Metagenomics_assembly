#!/bin/bash
#MSUB -N bmtagger_clean_reads.sh
#MSUB -A b1042
#MSUB -q genomics
#MSUB -m ae
#MSUB -M jamesgriffin2013@u.northwestern.edu
#MSUB -l nodes=1:ppn=24
#MSUB -l walltime=4:00:00
#MSUB -o /projects/b1042/Wells/Jim/Wells01/.o
#MSUB -e /projects/b1042/Wells/Jim/Wells01/fastQC.e

source activate env_bmtagger

# Set for your files
#base_dir=/projects/b1042/Wells/Jim/Wells01/clean_reads/
#temp_dir= "${base_dir}bmtagger_temp"
#file_out="${base_dir}all_reads_no_human.fastq"
#fwd_reads="${base_dir}all_reads_1.fastq"
#rev_reads="${base_dir}all_reads_2.fastq"

# Human Genome Indexed Database Files, Don't change!
#hg38_bitmask = /projects/b1052/shared/DB/hg19/hg38.reference.bitmask
#hg38_srprism = /projects/b1052/shared/DB/hg19/hg38.reference.srprism
#bmtagger.sh -b $hg38_bitmask -x $hg38_srprism -T $temp_dir -q 1 -1 $fwd_reads -2 $rev_reads -o $file_out

cd /projects/b1042/Wells/Jim/Wells01/clean_reads

bmtagger.sh -b /projects/b1052/shared/DB/hg19/hg38.reference.bitmask -x /projects/b1052/shared/DB/srprism-db -q 1 -1 all_reads_1.fastq -2 all_reads_2.fastq -o all_reads_no_human
