"""
assembly_job_writer:
    This script writes (and submits) batch job submission scripts for several
    different assemblers. This code and other workflow scripts in this directory
    are designed to take user input and a separate parameter file and perform
    repeatable analyses.

Parameter File:
    The parameter file is a json file with key-value pairs that holds options for
    assembly as well as other processes. Make sure to edit the parameter file and
    save it on Quest if you want to change any of the default options.

User Input:
    --assembler ASSEMBLER, -a ASSEMBLER
                        Comma separated assemblers to use
                        ('idba','spades','megahit')
  --raw RAW, -r RAW
                        Input directory containing subdirectories with
                        demultiplexed paired end reads
  --out OUT, -o OUT     Directory where assembled reads will go i.e.:
                        out/METHOD/SAMPLE
  --params PARAMS, -p PARAMS
                        A file of tab separated key-val parameters for each
                        method
  --check_modules       Runs validate_modules to check module paths
  --nosub               writes shell scripts but doesn't submit jobs

New Methods:
    This script can be easily extended to new assembly methods by adding a valid
    method to load a module and run it on a given sample.
    See the megahit and idba functions for examples.

TODO:
Found out about snakemake, I really should have just done that lol.
"""
__author__ = "Jim Griffin"
__author_email__ = "jimbogrif@gmail.com"

# Imports:
import argparse
import json
import os
import subprocess
import sys


def main():
    """
    Assembles samples with the selected method(s)
    """
    # Parse Arguments, Assembly Methods, and Load Parameters
    args = parse_arguments()
    params = read_json_config(args.params)
    methods = args.assembler.split(',')
    if args.check_modules:
        subprocess.call(["python /projects/b1052/shared/validate_modules.py"])

    # Get Sample Directory Paths found in input argument
    sample_dirs = [f.path for f in os.scandir(args.input) if f.is_dir()]

    # Create job scripts for each method-sample pair
    for method in methods:
        for sample in sample_dirs:
            # Generate names and paths
            sample_id = sample.split('/')[-1]
            outdir = make_output_dir(args.out, method, sample_id)
            # Write job submission script
            job_script = write_assembly_script(args.raw, sample_id, method,
                                               outdir, params)
        if args.nosub is False:
            submit_job(job_script)


def write_assembly_script(input_dir, sample_id, method, outdir, params):
    """Writes a single MOAB job submission script """
    print("Writing MSUB submission script\n\
    Method: {0}\nSample ID: {1}".format(method, sample_id))

    filename = "{0}{1}.run.{2}.sh".format(outdir, sample_id, method)
    with open(filename, "w") as file_out:
        header_info = fill_header(sample_id, outdir, params["email"])
        # Makes a lot of assumptions about paths:
        forward = "{0}{1}_R1_trimmed.fq".format(input_dir, sample_id)
        reverse = "{0}{1}_R2_trimmed.fq".format(input_dir, sample_id)

        if method == 'idba':
            cmd_lines = idba(forward, reverse, params, sample_id, input_dir, outdir)
        if method == 'megahit':
            cmd_lines = megahit(forward, reverse, params, outdir)

        msub_lines = header_info + cmd_lines
        for line in msub_lines:
            file_out.write(line)

    return filename

def megahit(forward, reverse, params, outdir):
    """
    Writes lines to job submission script to assemble 1 sample using megahit.
    """
    load_megahit_cmd = "module load megahit/1.0.6.1"
    megahit_cmd = "megahit --k-min {0} --k-max {1} --k-step {2} -t {3} -m {4} \
    -1 {5} -2 {6} -o {7}".format(params['megahit']['k-min'],
                                 params['megahit']['k-max'],
                                 params['megahit']['k-step'],
                                 params['megahit']['t'],
                                 params['megahit']['m'], forward,
                                 reverse, outdir)
    return load_megahit_cmd, megahit_cmd

def idba(forward, reverse, params, sample_id, input_dir, outdir):
    """
    Writes lines to job submission script to assemble 1 sample using idba.
    """
    # Read names:
    fa_dir = input_dir+"fa_interleaved/"
    fa_reads = "{0}{1}_trimmed.fa".format(fa_dir, sample_id)
    # Commands:
    fq2fa_cmd = "fq2fa --merge --filter {0} {1} {2}\n\n".format(forward, reverse,
                                                                fa_reads)
    idba_cmd = ("idba idba_ud -r {0} --mink {1} --maxk {2} --step {3}\
            --min_contig {4} --num_threads {5} -o {6}/"
                .format(fa_reads, params['idba']['mink'], params['idba']['maxk'],
                        params['idba']['step'], params['idba']['min_contig'],
                        params['idba']['num_threads'], outdir
                       ))
    idba_load_cmd = "module load idba/2016_12_longread\n\n"
    return idba_load_cmd, fq2fa_cmd, idba_cmd

def fill_header(filename, assembled_dir, email):
    """
    Creates an msub header for job submission using the filename specified.
    """
    return [
        "#!/bin/bash\n",
        "#MSUB -N {}\n".format(filename),
        "#MSUB -A b1042\n",
        "#MSUB -q genomics\n",
        "#MSUB -m ae\n",
        "#MSUB -M {}".format(email),
        "#MSUB -l nodes=1:ppn=24\n",
        "#MSUB -l walltime=48:00:00\n\n",
        "cd {}\n\n".format(assembled_dir)
    ]

def submit_job(job_file):
    """
    Submits a job in a new shell using msub and subprocess.
    Remember that modules need to be loaded in that shell.
    """
    shell_script = "msub {}".format(job_file)
    subprocess.check_call([shell_script], shell=True)

def read_params(params):
    """
    Reads parameters in from a tab separated text file and returns a
    dictionary of key-value pairs
    """
    p_vals = {}
    with open(params) as param_f:
        for line in param_f:
            if not line.startswith('#'):
                cols = line.split()
                p_vals[cols[0]] = cols[1]
    return p_vals

def read_json_config(config):
    """
    Fixed my json object handling!
    """
    with open(config) as config_f:
        params = json.load(config_f)
    return params

def parse_arguments():
    """Parse and return command line arguments"""
    parser = argparse.ArgumentParser()
    #Default options:
    parser.add_argument("--assembler", "-a", type=str,
                        help="Comma separated assemblers to use\
                        ('idba','spades','megahit')")
    parser.add_argument("--raw", "-r", type=str,
                        help="Input directory containing subdirectories with \
                        demultiplexed paired end reads")
    parser.add_argument("--out", "-o", type=str,
                        help="Directory where assembled reads will go \
                        i.e.: out/METHOD/SAMPLE")
    parser.add_argument("--params", "-p", type=str,
                        help="A file of tab separated key-val parameters \
                        for each method")
    parser.add_argument("--check_modules", action="store_true",
                        help="Runs validate_modules to check module paths")
    parser.add_argument("--nosub", action="store_true",
                        help="writes shell scripts but doesn't submit jobs")

    # Print help if no options given.
    if len(sys.argv) == 1:
        parser.print_help()
        print("\n\nNeed command line input\n\n")
        sys.exit(1)

    # Parse Command Line Arguments:
    try:
        result = parser.parse_args()
        return result
    except Exception as e:
        parser.print_help()
        print(e)
        sys.exit(0)

def make_output_dir(out, method, sample):
    """
    Creates a sample and method specific output directory if it does not exist
    and returns the filepath
    """
    full_output_path = "{0}{1}/{2}".format(out, method, sample)
    if not os.path.exists(full_output_path):
        os.makedirs(full_output_path)
    return full_output_path


if __name__ == "__main__":
    main()
