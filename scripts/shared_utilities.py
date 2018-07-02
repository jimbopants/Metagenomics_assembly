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
import datetime
import os
import subprocess
import sys
import yaml

CONFIGFILE = "./config.yaml"

def read_config(config=CONFIGFILE):
    """
    Reads parameters in from a yaml file and returns a
    dictionary of key-value pairs
    """
    params = yaml.load(open(config))
    params = make_relative_paths(params)
    params["methods"]["assembly"] = params["methods"]["assembly"].split()
    params["methods"]["map"] = params["methods"]["map"].split()
    return params

def make_relative_paths(params):
    for path in params["paths"]:
        params["paths"][path] = "{0}{1}".format(params["base_dir"], params["paths"][path])
    return params

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

def make_output_dir(out, method, sample):
    """
    Creates a sample and method specific output directory if it does not exist
    and returns the filepath. Does not have an ending / slash
    """
    if sample == None:
        full_output_path = "{0}{1}".format(out, method)
    else:
        full_output_path = "{0}{1}/{2}".format(out, method, sample)
    if not os.path.exists(full_output_path):
        os.makedirs(full_output_path)
    return full_output_path
