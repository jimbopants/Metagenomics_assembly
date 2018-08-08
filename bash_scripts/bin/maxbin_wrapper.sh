#!/bin/bash
#MSUB -N maxbin_wrapper_7.24.18.sh
#MSUB -A b1042
#MSUB -q genomics
#MSUB -m ae
#MSUB -M jamesgriffin2013@u.northwestern.edu
#MSUB -l nodes=1:ppn=24
#MSUB -l walltime=48:00:00
#MSUB -o /projects/b1042/Wells/Jim/Wells01/binning/binning-maxbin_7.24.18.o
#MSUB -e /projects/b1042/Wells/Jim/Wells01/binning/binning-maxbin_7.24.18.e

# Lazy but <2 months til I graduate so fuck it hardcode everything!
contigs=/projects/b1042/Wells/Jim/Wells01/assembled/co_assembly_stored/contig.fa
outdir=/projects/b1042/Wells/Jim/Wells01/binning_maxbin2/
readspath=/projects/b1042/Wells/Jim/Wells01/fa_interleaved
readfiles=(C1A_S1_trimmed.fa C1B_S2_trimmed.fa C2A_S3_trimmed.fa C2B_S4_trimmed.fa C3A_S5_trimmed.fa C3B_trimmed.fa)

cd /projects/b1052/shared/MaxBin-2.2.5/
perl run_MaxBin.pl -contig $contigs -out $outdir -thread 24 -reads $readspath/${readfiles[0]} -reads2 $readspath/${readfiles[1]} -reads3 $readspath/${readfiles[2]} -reads4 $readspath/${readfiles[3]} -reads5 $readspath/${readfiles[4]} -reads6 $readspath/${readfiles[5]}
