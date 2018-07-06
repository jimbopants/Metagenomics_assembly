"""
assembly_job_writer:
    This script writes (and submits) batch job submission scripts for several
    different assemblers. See the readme on github for more help.

By default, assumes you have a directory of demultiplexed samples like:
    PROJECTDIR/
        clean_reads/
            S1/{forward, reverse}
            S2/

Usage:
    1. Set the following parameters in the config.yaml file before running:
    Including:
        sub (submit jobs Yes/No?)
        co-assembly / sample-assemblys
        Job header info: email, allocation, etc.
        base_dir -> Your project directory
        paths/:
            clean_reads
            assembled
            clean_suffix
            fa_interleaved
        parameters/:
            idba
            megahit
        Methods/assembly:
            Comment out any you don't want to use

    2. Run the script at the command line or as part of a workflow script.

Arguments (Optional, will override config file, mostly for convenience/speed):
    --directory, -d     Directory containing per-sample folders to process.

    --method, -m        idba,megahit (comma separated list)

    --sample, -s        Path to single sample directory. Forces script to operate on a
                        single sample only

    --out OUT, -o OUT   Output filepath

    --sub SUB           (Boolean: True/False) (overrides config)
                        Writes shell scripts but doesn't submit jobs

Config File:
    The config file is a yaml file with key-value pairs that holds options for
    assembly as well as other processes. Make sure to edit the parameter file and
    save it on Quest if you want to change any of the default options.

Adding New Methods:
    This script can be easily extended to new assembly methods by doing 2 things:
    1. Add an assembly method function:
        ASSSEMBLER(forward, reverse, params, outdir, sample):
    2. Add the assembler name to the config.yaml
    3. Add {"assembler":assembler} to the assemblers dictionary at the beginning of main
    4. See the megahit and idba functions for examples.

"""
__author__ = "Jim Griffin"
__author_email__ = "jimbogrif@gmail.com"

# Imports:
import glob
import shared_utilities as ll

def main():
    """
    Assembles samples with the selected method(s)
    """
    # Load Parameters
    params, args = ll.load_params_and_input('assembly')
    # Maps assemblers to functions. Comment out assemblers in config.yaml or override with -m at runtime
    assemblers = {"megahit" : megahit, "idba" : idba}

    # Returns a list of [path, sampleID] for each sample to assemble
    # From command line input (-s or -d or config file)
    if args.sample:
        samples_to_assemble = [[args.sample, args.sample.rsplit('/', 1)[1]]]
    else:
        samples_to_assemble = ll.paths_and_samples(params["paths"]["clean"])

    # Create job scripts for each method-sample pair
    for method in params["methods"]["assembly"]:
        for sample in samples_to_assemble:
            # Generate job script name and file out paths
            # megahit can't recursively make paths or write to existing paths...
            outdir = "{0}{1}/{2}".format(params["paths"]["assembly"], method, sample[1])
            outdir = ll.make_dir("{}{}/".format(params["paths"]["assembly"], method)
            job_out = ll.make_dir(params["paths"]["assembly"], "job_logs", None)
            job_script = "{0}/{1}.assemble.{2}.sh".format(job_out, sample[1], method)

            # Open job script for writing and write header, then add method-specific lines
            with open(job_script, "w") as file_out:
                print("Writing MSUB submission script\n\
                Jobfile: {}\n".format(job_script))
                header_info = ll.fill_header(job_script, params)
                # todo: this might be a relative path. update to absolute.
                forward, reverse = sorted(glob.glob(sample[0]+"/*.fq")[:2])
                cmds = assemblers[method](forward, reverse, params, outdir, sample)
                # Write lines to job:
                msub_lines = header_info + cmds
                for line in msub_lines:
                    file_out.write(line)
            if params["sub"] is True:
                ll.submit_job(job_script)

def megahit(forward, reverse, params, outdir, sample):
    """
    Writes lines to job submission script to assemble 1 sample using megahit.
    """
    load_megahit_cmd = "module load {0}\n".format(params["modules"]["megahit"])
    megahit_cmd = "megahit --k-min {0} --k-max {1} --k-step {2} -t {3} -m {4} \
    -1 {5} -2 {6} -o {7}".format(params['megahit']['k-min'],
                                 params['megahit']['k-max'],
                                 params['megahit']['k-step'],
                                 params['megahit']['t'],
                                 params['megahit']['m'],
                                 forward,
                                 reverse,
                                 outdir)
    return [load_megahit_cmd, megahit_cmd]

def idba(forward, reverse, params, outdir, sample):
    """
    Writes lines to job submission script to assemble 1 sample using idba.
    """
    fa_reads = "{0}{1}_trimmed.fa".format(params["paths"]["raw_interleaved"], sample[1])
    idba_load_cmd = "module load {}\n\n".format(params["modules"]["idba"])
    fq2fa_cmd = "fq2fa --merge --filter {0} {1} {2}\n\n".format(forward,
                                                                reverse,
                                                                fa_reads)
    idba_cmd = ("idba idba_ud -r {0} --mink {1} --maxk {2} --step {3}\
            --min_contig {4} --num_threads {5} -o {6}/"
                .format(fa_reads,
                        params['idba']['mink'],
                        params['idba']['maxk'],
                        params['idba']['step'],
                        params['idba']['min_contig'],
                        params['num_threads'],
                        outdir
                       ))
    return [idba_load_cmd, fq2fa_cmd, idba_cmd]

if __name__ == "__main__":
    main()
