#!/bin/bash
#MSUB -N concoct_preprocess2.sh
#MSUB -A b1042
#MSUB -q genomics
#MSUB -m ae
#MSUB -M jamesgriffin2013@u.northwestern.edu
#MSUB -l nodes=1:ppn=24
#MSUB -l walltime=48:00:00
#MSUB -o /projects/b1042/Wells/Jim/Wells01/logs/concoct_preprocess2.o
#MSUB -e /projects/b1042/Wells/Jim/Wells01/logs/concoct_preprocess2.e

# Written by JG 7/25/18
# updaed 7/31/18

# filenames:
cd /projects/b1042/Wells/Jim/Wells01/
coassembly=/projects/b1042/Wells/Jim/Wells01/assembled/co_assembly_stored
concoct=/projects/b1052/shared/CONCOCT-0.4.0
samples=/projects/b1042/Wells/Jim/Wells01/clean_reads
basedir=/projects/b1042/Wells/Jim/Wells01/CONCOCT_input/
maxbins=200



# Source CONCOCT programs
# 7/31: Bash tip of the day:
# Need to export (not with $) in order for subprocesses to see the variable.
source activate env_concoct
source source_concoct_paths.sh
MRKDUP="/software/picard/2.6.0/picard-tools-2.6.0/picard.jar MarkDuplicates"
export MRKDUP

# cut up contigs into 10-20kbp length
#python $concoct/scripts/cut_up_fasta.py -c 10000 -o 0 -m $coassembly/contig.fa > $coassembly/contig_c10K.fa

# map w/ bowtie2:
#bowtie2-build $coassembly/contig_c10K.fa $coassembly/contig_c10K.fa

# for each sample, run CONCOCT's mapping pipeline:
# 1. Maps given paired library to given reference with bowtie2
# 2. Use picard to remove duplicates.
# 3. Requires:
#       MRKDUP to point to picard/MarkDuplicates
#       bowtie2
#       BEDtools
# 4. options:
# -t      Number of threads for bowtie2 and the java garbage collector
# -c      Calculate coverage with BEDTools
# 5. usage: bash `basename $0` [options] <reads1> <reads2> <qname> <ref> <rname> <outdir>
# bash tips: man to view command help. basename strips pre/suffix. mkdir -p does int. Paths
# sed s/r1/r2 replaces r1 with r2 in s.
cd $samples
# 7/31: This didn't run because I didn't export MRKDUP to environment:
for f in $samples/C*/*R1_trimmed.fq; do
    echo $f
    mkdir -p map/$(basename -s _R1_trimmed.fq $f);
    cd map/$(basename -s _R1_trimmed.fq $f);
    bash $concoct/scripts/map-bowtie2-markduplicates.sh -ct 24 -p '-q' $f $(echo $f | sed s/R1/R2/) pair $coassembly/contig_c10K.fa asm bowtie2;
    cd ../..;
done

# Generate an abundance/coverage input table for CONCOCT:
cd ${samples}/map
python ${concoct}/scripts/gen_input_table.py --isbedfiles \
    --samplenames <(for s in C*; do echo $s | cut -d'_' -f1; done) \
    $coassembly/contig_c10K.fa */bowtie2/asm_pair-smds.coverage \
> concoct_inputtable.tsv
mkdir $basedir
mv concoct_inputtable.tsv $basedir

# Generate a linkage table:
cd ${samples}/map
python $concoct/scripts/bam_to_linkage.py -m 24 \
    --regionlength 500 --fullsearch \
    --samplenames <(for s in C*; do echo $s | cut -d'_' -f1; done) \
    $coassembly/contig_c10K.fa C*/bowtie2/asm_pair-smds.bam \
> concoct_linkage.tsv
mv concoct_linkage.tsv $basedir


# Subset coverage to mean for each contig:
cd $basedir
cut -f1,3- concoct_inputtable.tsv > concoct_inputtableR.tsv

# Run Concoct on 150 bp reads with up to a max of $maxbin bins
concoct -c $maxbins -r 150 --coverage_file concoct_inputtableR.tsv --composition_file $coassembly/contig_c10K.fa -b concoct-output/
