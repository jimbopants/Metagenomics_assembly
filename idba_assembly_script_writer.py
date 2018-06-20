# Written by Gao Han, updated by Jim Griffin
# 6/20/18

# Imports:
import subprocess
import os

# hard-coded filepaths for now:
sample_id = ['C1A_S1', 'C1B_S2', 'C2A_S3', 'C2B_S4', 'C3A_S5', 'C3A_S6']
trim_reads = "/projects/b1042/Wells/Jim/Wells01/clean_reads/"

#output directories:
assembled_dir = "/projects/b1042/Wells/Jim/Wells01/assembled/"
fa_dir =        "/projects/b1042/Wells/Jim/Wells01/fa_interleaved/"

# Make directories:
for dir in [assembled_dir, fa_dir]:
    if not os.path.exists(dir):
        os.makedirs(dir)

def fill_header(sample, assembled_dir):
    """Creates an msub header for job submission using the filename specified."""
    return [
    "#!/bin/bash\n",
    "#MSUB -N {}_assemble_idba\n".format(sample),
    "#MSUB -A b1042\n",
    "#MSUB -q genomics\n",
    "#MSUB -m abe\n",
    "#MSUB -l nodes=1:ppn=24\n",
    "#MSUB -l walltime=48:00:00\n\n",
    "module load idba/2016_12_longread\n\n",
    "cd {}\n\n".format(assembled_dir)
    ]

for sample in sample_id:
    print("Writing MSUB script for sample ID: {}".format(sample))
    filename = "{}.run.assemble.idba.sh".format(sample)
    with open(filename, "w") as f:
        header_info = fill_header(sample, assembled_dir)
        for line in header_info:
            f.write(line)
    # Read names:
        forward = "{0}{1}_R1_trimmed.fq".format(trim_reads, sample)
        reverse = "{0}{1}_R2_trimmed.fq".format(trim_reads, sample)
        fa = "{0}{1}_trimmed.fa".format(fa_dir, sample)
    # Commands:
        fq2fa_cmd = "fq2fa --merge --filter {0} {1} {2}\n\n".format(forward, reverse, fa)
        idba_cmd = "idba idba_ud -r {0} --mink 20 --maxk 80 --step 10 --min_contig 500 --num_threads 24 -o {1}{2}/\n".format(fa, assembled_dir, sample)
        f.write(fq2fa_cmd)
        f.write(idba_cmd)

    shell_script = "msub {}".format(filename)
    subprocess.check_call([shell_script], shell=True)
    print('MSUB script submitted for sample ID: {}'.format(sample))
