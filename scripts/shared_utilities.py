"""
shared_utilities:
    This script contains functions shared across the job_writing scripts.
    I usually import it as ll for low level.

## Functions:
read_config:
    Reads a config.yaml and saves the parameters as a dictionary

fill_header:
    Writes all of the generic properties of an msub submission script

submit_job:
    Submits a job using msub

TODO:
    other stuff as needed
"""
__author__ = "Jim Griffin"
__author_email__ = "jamesgriffin2013@u.northwestern.edu"

# Imports:
import os
import subprocess
import sys
import yaml

CONFIGFILE = "../parameters/config.yaml"

def read_config(config=CONFIGFILE):
    """
    Reads parameters in from a yaml file and returns a
    dictionary of key-value pairs
    """
    params = yaml.load(open(config))
    return params

def fill_header(filename, params, output):
    """
    Creates a list of msub header strings that can be pasted into a job submission script.
    Input:
        filename: Name of quest job
        params: dictionary of params, typically loaded from config.yaml although
        these values can be overwritten by a job submission script before passing to fill_header.
        output: Filename where stdout and stderr is written
    """
    return [
        "#!/bin/bash\n",
        "#MSUB -N {}\n".format(filename),
        "#MSUB -A {}\n".format(params["allocation"]),
        "#MSUB -q {}\n".format(params["queue"]),
        "#MSUB -m ae\n",
        "#MSUB -M {}".format(params["email"]),
        "#MSUB -l {}\n".format(params["nodes"]),
        "#MSUB -l walltime={}".format(params["walltime"]),
        "#MSUB -j {}\n\n".format(output)
    ]

def submit_job(job_file):
    """
    Submits a job in a new shell using msub and subprocess.
    Remember that modules need to be loaded in that shell.
    """
    shell_script = "msub {}".format(job_file)
    subprocess.check_call([shell_script], shell=True)

def paths_and_samples(directory):
    """
    Returns a zipped list of [(path, samplename), (path2, name2), ...]
    based on the subdirectories in DIRECTORY. Mostly used to iterate over all samples
    when applying a processing step to all samples
    """
    sample_dirs = [f.path for f in os.scandir(directory) if f.is_dir()]
    sample_ids = [sample.split('/')[-1] for sample in sample_dirs]
    zipped_items = zip(sample_dirs, sample_ids)
    return zipped_items

def make_output_dir(out, method, sample):
    """
    Creates a sample and method specific output directory if it does not exist
    and returns the filepath. Does not have an ending / slash
    """
    full_output_path = "{0}{1}/{2}".format(out, method, sample)
    if not os.path.exists(full_output_path):
        os.makedirs(full_output_path)
    return full_output_path
