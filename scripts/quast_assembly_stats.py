"""
quast_assembly_stats:
    This script calculates assembly stats using quast for all assemblies.
    It doesn't use metaquast as it seems useless.
    Quast documentation: http://quast.bioinf.spbau.ru/manual.html
    By default, it loads assembled contigs from the paths listed below.
    See the readme on github for more help.

Usage:
    1. Set the following parameters in the config.yaml file before running:
    paths
        assembled
            idbda,megahit
        quast_results
    2. If running from command line you can optionally pass user input to override the config file
    (useful for running 1 sample at a time or overriding sub without needing a text editor)

Arguments (Optional):
    --directory, -d     Directory containing per-sample folders to process.
                        (overrides config)

    --sample, -s        Path to single file. Forces script to operate on a
                        single file only

    --out OUT, -o OUT   Output filepath (overrides config)

    --sub SUB           Writes shell scripts but doesn't submit jobs
                        (overrides config)
"""
__author__ = "Jim Griffin"
__author_email__ = "jimbogrif@gmail.com"

# Imports:
import datetime
import os
import shared_utilities as ll

def main():
    params, args = ll.load_params_and_input("assembly")
    contig_names = {"idba" : "contig.fa",
                    "megahit" : "final.contigs.fa"}
    # If -sample is passed, just check assembly stats for that file,
    # otherwise expand the directory based on config or -directory
    if args.sample is None:
        quast_list = []
        for method in params["methods"]["assembly"]:
            samples = ll.paths_and_samples("{0}{1}/".format(params["paths"]["assembly"], method))
            full_paths =  ["{}/{}".format(x[0], contig_names[method]) for x in samples]
            quast_list.extend(full_paths)
        files_to_quast = " ".join(quast_list)
    else:
        file_to_quast = args.sample

    if not os.path.exists(params["paths"]["quast_out"]):
        os.makedirs(params["paths"]["quast_out"])

    timestamp = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d')
    job_script = "{0}{1}.quast.sh".format(params["paths"]["quast_out"], timestamp)
    #Actually it's easier to just run on a head node.
    header_info = ll.fill_header(job_script, params)

    # Could add support/option for metaquast later, I mostly found it useless.
    quast_lines = ["{0} {1} -o {2} --threads {3} -L".format(
                                                        params["binaries"]["quast"]
                                                        files_to_quast,
                                                        params["paths"]["quast_out"],
                                                        params["num_threads"])]
    # Write job script
    cmd_lines = header_info + quast_lines
    with open(job_script, 'w') as file_out:
        for line in cmd_lines:
            file_out.write(line)

    if params["sub"] is True:
        ll.submit_job(job_script)

if __name__ == "__main__":
    main()
