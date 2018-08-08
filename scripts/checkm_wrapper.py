"""
CheckM Wrapper
Runs the checkM lineage_wf workflow command:
https://github.com/Ecogenomics/CheckM/wiki/Quick-Start

Command line arguments:
    --directory, -d       Directory containing bins.

    --ext, -x        Bin file extension (default fna)

    --out, -o        Output filepath for bin statistics

    --sub SUB        (Boolean: True/False) (overrides config)
                     Submit jobs

"""
__author__ = "Jim Griffin"
__author_email__ = "jamesgriffin2013@u.northwestern.edu"
# Imports:

import os
import datetime
import shared_utilities as ll


def main():
    """ Make directories, parse input
        Create a job submission script,
        Write standard header
        Write checkM lineage workflow commands
        Submit job if sub parameter is True
    """
    params, args = ll.load_params_and_input(None)
    os.makedirs(params["paths"]["checkM"], exist_ok=True)
    os.makedirs(params["paths"]["jobs"], exist_ok=True)

    now = datetime.datetime.strftime(datetime.datetime.now(), '%d_%m_%Y-%H-%M')
    job_script = "{0}/checkM_job_{1}.sh".format(params["paths"]["jobs"], now)

    with open(job_script, "w") as file_out:
        msub_lines = ll.fill_header(job_script, params)
        msub_lines += lineage_wf_wrapper(args, params)
        for line in msub_lines:
            file_out.write(line)
        if params["sub"] in (True, "True", "T", "true", 1):
            ll.submit_job(job_script)


def lineage_wf_wrapper(args, params):
    """
    Returns lines to run checkM's lineage_wf (workflow) command with several additional flags
    Type checkm lineage_wf -h to see a full list of options
    """
    return [
        "source activate metawrap-env\n",
        "cd {}\n".format(params["base_dir"]),
        "checkm lineage_wf --tab_table -x {0} -t {1} --pplacer_threads {1} {2} {3}\n".format(
            args.ext,
            params['num_threads'],
            args.directory,
            args.out)]

if __name__ == "__main__":
    main()
