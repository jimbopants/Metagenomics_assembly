"""
shared_utilities:
    This script contains functions shared across the job_writing scripts.
    I usually import it as ll for low level.

## See individual funcion docstrings for usage
"""
__author__ = "Jim Griffin"
__author_email__ = "jamesgriffin2013@u.northwestern.edu"

# Imports:
import argparse
import datetime
import os
import subprocess
import sys
import yaml

CONFIGFILE = "./config.yaml"

################
# TOC
# 1. Config values and user Input
# 2. Path Management & Expansion
# 3. Job header and submission

################
# 1. Config values and user Input
def read_config(config=CONFIGFILE):
    """
    Reads parameters in from a yaml file and returns a
    dictionary of key-value pairs
    """
    params = yaml.load(open(config))
    params = make_relative_paths(params)
    for key in params["methods"]:
        if params["methods"][key]:
            params["methods"][key] = params["methods"][key].split()
    return params

def parse_arguments(process):
    """Parse and return command line arguments as an argparse args object"""
    parser = argparse.ArgumentParser()
    process_choices = {"assembly" : "[megahit,idba]",
                     "map" : "[bwa,bowtie]",
                     None: "Ignored for this script"
                     }
    parser.add_argument("--directory", "-d", type=str,
                        help="Directory containing per-sample folders to process. \
                        (overrides config)")
    parser.add_argument("--methods", "-m", type=str,
                        help="comma separated list of methods to use:\
                        Valid choices are {}".format(process_choices[process]))
    parser.add_argument("--out", "-o", type=str,
                        help="Output filepath (overrides config)")
    parser.add_argument("--sub",
                        help="writes shell scripts but doesn't submit jobs \
                        (overrides config)")
    parser.add_argument("--sample", "-s", nargs='?', const=True,
                        help="Path to single file. \
                        Forces script to operate on a single file only")
    parser.add_argument("--ext", "-x", default='fna',
                        help="Bin file extension. Only used by checkm_wrapper")

    try:
        result = parser.parse_args()
        return result
    except Exception as e:
        parser.print_help()
        print(e)
        sys.exit(0)

def load_params_and_input(process):
    """
    Overrides config.yaml values with command line input.
    Returns both args (command line arguments) & params (config.yaml with updated values)
    """
    params = read_config()
    args = parse_arguments(process)
    if args.sub:
        params["sub"] = args.sub
    if args.methods:
        params["methods"][process] = args.methods.split(',')
    if args.directory:
        params["paths"]["clean"] = args.directory # This is a horrible idea.
        params["paths"]["directory"] = args.directory # This is less bad
    if args.out:
        params["paths"][process] = args.out
        params["paths"]["out"] = args.out
    return params, args

################
# 2. Path Management & Expansion
def make_relative_paths(params):
    """Update config paths by prepending base directory"""
    for path in params["paths"]:
        params["paths"][path] = "{0}{1}".format(params["base_dir"], params["paths"][path])
    return params

def paths_and_samples(directory):
    """
    Returns a zipped list of [(path, samplename), (path2, name2), ...]
    based on the subdirectories in DIRECTORY. Mostly used to iterate over all samples
    when applying a processing step to all samples
    ie: /assembled/{C1A,C2A} -> [["/assembled/C1A", "C1A"], ["/assembled/C2A/", "C2A"]]
    """
    sample_dirs = [f.path for f in os.scandir(directory) if f.is_dir()]
    sample_ids = [sample.split('/')[-1] for sample in sample_dirs]
    zipped_items = list(zip(sample_dirs, sample_ids))
    return zipped_items

def make_dir(path_out):
    """
    Creates a directory if it does not exist and returns the filepath.
    """
    os.makedirs(path_out, exist_ok=True)
    return path_out

################
# 3. Job header and submission
def fill_header(filename, params):
    """
    Creates a list of msub header strings that can be pasted into a job submission script.
    Input:
        filename: Name of quest job
        params: dictionary of params, typically loaded from config.yaml although
        these values can be overwritten by a job submission script before passing to fill_header.
    """
    jobname = filename.rsplit("/",1)
    timestamp = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d-%H-%M')
    output_log = "{0}{1}.o".format(filename[:-3], timestamp)
    err_log = "{0}{1}.e".format(filename[:-3], timestamp)
    return [
        "#!/bin/bash\n",
        "#MSUB -N {}\n".format(jobname[-1]),
        "#MSUB -A {}\n".format(params["allocation"]),
        "#MSUB -q {}\n".format(params["queue"]),
        "#MSUB -m ae\n",
        "#MSUB -M {}\n".format(params["email"]),
        "#MSUB -l {}\n".format(params["nodes"]),
        "#MSUB -l walltime={}\n".format(params["walltime"]),
        "#MSUB -o {}\n".format(output_log),
        "#MSUB -e {}\n\n".format(err_log)
    ]

def submit_job(job_file):
    """
    Submits a job in a new shell using msub and subprocess.
    Remember that modules need to be loaded in that shell.
    """
    print("Submitting MSUB submission script\n\
    Jobfile: {}".format(job_file))
    shell_script = "msub {}".format(job_file)
    subprocess.check_call([shell_script], shell=True)
