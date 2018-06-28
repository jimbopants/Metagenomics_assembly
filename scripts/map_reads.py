"""
BWA/bowtie wrapper
This script writes and submits jobs to map reads using bwa or bowtie as specified by the config file

User Input:
    User input is optional for this script.
    If 2 command line arguments are passed:
        arg1: a directory w/ per sample subdirectories of assembled contigs
        arg2: a directory w/ per sample subdirectories of mapped contigs
    If command line arguments are not present it takes the assembled and mapped directories from the config file

Mapping to a co-assembly:
    I've often seen co-assembly performed, and then per-sample reads to the co-assembly
    I think the idea is you get the biggest assembly and then get per sample abundances
    To handle this, any sample named co-assembly should be ignored until the end and then used as an index to map all samples to
"""
__author__ = "Jim Griffin"
__author_email__ = "jamesgriffin2013@u.northwestern.edu"

# Imports:
import subprocess
import sys
import glob
import shared_utilities as ll
# unused:
import argparse
import json
import os

def main():
    """
    Map samples with the selected method(s)
    """
    # Load parameters and Optional command line input to overwrite config.yaml
    params = ll.read_config()
    if len(sys.argv) > 1:
        params["paths"]["assembly"] = sys.argv[1]
        params["paths"]["map"] = sys.argv[2]

    # Get sample paths + full_names
    samples_to_map = ll.paths_and_samples(params["paths"]["assembly"])

    # Map all samples with all active mapping methods from config file:
    for method in params["methods"]["map"].split():
        for sample in samples_to_map:
            outdir = ll.make_output_dir(params["paths"]["map"], method, sample[1])
            job_script = "{0}{1}.map.{2}.sh".format(outdir, sample_id, method)

            # Write job script:
            with open(job_script, "w") as file_out:
                header_info = ll.fill_header(job_script, params, filename[:-3]+".log")

                if method == "bwa":
                    cmd_lines = bwa(params, sample, outdir)
                if method == "bowtie":
                    cmd_lines = bowtie()
                # Merge header and command specific text:
                msub_lines = header_info + cmd_lines
                for line in msub_lines:
                    file_out.write(line)

            # Run job if turned on in config:
            if params["run"]:
                ll.submit_job(job_script)


def bwa(params, sample, outdir):
    """Writes lines to job submission script to map 1 sample using bwa
    """
    index = sample[0]+"contig.fa"
    trimmed_reads = "{0}{1}/*".format(params["paths"]["trim"], sample[1])
    fwd, rev = glob.glob(trimmed_reads)[:2]
    return [
            "module load {}".format(params["modules"]["bwa"]),
            "bwa index {}".format(index),
            "bwa mem {0} {1} {2} > {3}/aln-pe.sam".format(index, fwd, rev, outdir)
            ]

def bowtie():
    """Writes map script with bowtie"""
    pass

if __name__ == "__main__":
    main()
