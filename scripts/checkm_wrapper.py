"""
CheckM Wrapper
*** Update: Never do this it's a huge waste of time.
Run something once before writing anything like this again!

Submits a job script to check bin completeness/contamination with checkM.
Paths: Give relative paths to the params[base_dir] for -d and -o

By default, splits bins into groups of 20 and submits batch jobs for each group.

Runs the lineage_wf workflow command right now:
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
import glob
import os
import shutil
import shared_utilities as ll

def main():
    """ Handles parameters and splits job batches
    """
    params, args = ll.load_params_and_input(None)
    os.makedirs(params["paths"]["out"], exist_ok=True)
    bins = glob.glob("{0}/*{1}".format(args.directory, args.ext))
    bin_index = 0

    while len(bins) > 20:
        bin_subset = bins[:20]
        del bins[:20]
        sub_bin_dir, sub_out_dir = move_bin_subset(bin_subset, bin_index, params)
        checkm_wrap(params, args, sub_bin_dir, sub_out_dir, bin_index)
        bin_index += 1
    if len(bins) <= 20:
        sub_bin_dir, sub_out_dir = move_bin_subset(bins, bin_index, params)
        checkm_wrap(params, args, sub_bin_dir, sub_out_dir, bin_index)
    return

def checkm_wrap(params, args, bin_dir, out_dir, bin_index):
    """ Create a job submission script in OUT_DIR,
        Write standard header
        Write checkM lineage workflow commands
        Submit job if sub parameter is True
    """
    job_script = "{0}/checkM_batch_{1}.sh".format(params["paths"]["directory"], bin_index)
    with open(job_script, "w") as file_out:
        msub_lines = ll.fill_header(job_script, params)
        msub_lines += lineage_wf_wrapper(args.ext, params, bin_dir, out_dir)
        for line in msub_lines:
            file_out.write(line)
        if params["sub"] in (True, "True", "T", "true", 1):
            ll.submit_job(job_script)

def move_bin_subset(bin_subset, bin_index, params):
    """
    Make a new subdirectory and move the bin subset to it.
    Make a subdirectory for results for that section
    Return the new bin directory.
    """
    new_bin_dir = "{0}bin_subset_{1}".format(params["paths"]['directory'], bin_index)
    new_out_dir = "{0}/checkM_subset_{1}".format(params["paths"]['out'], bin_index)
    os.makedirs(new_bin_dir, exist_ok=True)
    for bins in bin_subset:
        name = bins.rsplit('/',1)[1]
        shutil.move(bins, "{0}/{1}".format(new_bin_dir, name))
    return new_bin_dir, new_out_dir

def lineage_wf_wrapper(file_ext, params, bin_dir, out_dir):
    """
    Returns lines to run checkM's lineage_wf (workflow) command with several additional flags
    Type checkm lineage_wf -h to see a full list of options
    """
    return [
        "source activate metawrap-env\n",
        "cd {}\n".format(params["base_dir"]),
        "checkm lineage_wf --tab_table -x {0} -t {1} --pplacer_threads {2} {3} {4}\n".format(
            file_ext,
            params['num_threads'],
            params['num_threads'],
            bin_dir,
            out_dir)]

if __name__ == "__main__":
    main()
