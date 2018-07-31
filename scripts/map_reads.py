"""
TODO:
    Add samtools view -> bam, remove sam file.
    Check whether I have header info or not from bwa.


BWA/bowtie wrapper:
    This script writes and submits jobs to map reads using bwa or bowtie
    as specified by the config file. See the readme on github for more help.

Usage:
    1. Check or set the following parameters in the config.yaml file before running:
    Including:
        base_dir -> Your project directory
        paths/:
            assembled
            map
            reads_to_map
            contig_names
        Methods/map: Comment out any methods you don't want to use
            BWA
            bowtie

    2. Run the script at the command line or as part of a workflow script.
    3. Mapping to a co-assembly:
        Manually override with --sample

Arguments (Optional, will override config file, mostly for convenience/speed):
    --directory, -d     Directory containing per-sample folders to process.

    --method, -m        idba,megahit (comma separated list)

    --sample, -s        Path to single sample directory. Forces script to operate on a
                        single sample only

    --out OUT, -o OUT   Output filepath

    --sub SUB           (Boolean: True/False) (overrides config)
                        Writes shell scripts but doesn't submit jobs
    If command line arguments are not present it takes the assembled and mapped directories from the config file
"""
__author__ = "Jim Griffin"
__author_email__ = "jamesgriffin2013@u.northwestern.edu"

# Imports:
import glob
import os
import shared_utilities as ll

def main():
    """
    Map samples with the selected method(s)
    """
    # Load parameters and optional command line input to overwrite config.yaml
    params, args = ll.load_params_and_input("map")
    mappers = {"bwa" : bwa, "bowtie" : bowtie}
    os.makedirs(params["paths"]["map"], exist_ok=True)

    # Returns a list of [path, sampleID] for each sample to assemble
    # From command line input (-s or -d or config file)
    if args.sample:
        samples_to_assemble = [[args.sample, args.sample.rsplit('/', 1)[1]]]
    else:
        samples_to_assemble = ll.paths_and_samples(params["paths"]["reads_to_map"])

    # Map all samples with all active mapping methods from config file:
    for method in params["methods"]["map"]:
        for sample in samples_to_map:
            outdir = ll.make_dir("{0}{1}".format(params["paths"]["map"], sample[1]))
            job_script = "{0}/{1}.{2}.sh".format(outdir.rsplit('/',1)[0], sample[1], method)
            # Write job script:
            with open(job_script, "w") as file_out:
                msub_lines = ll.fill_header(job_script, params)
                msub_lines += mappers[method](params, sample, outdir)
                for line in msub_lines:
                    file_out.write(line)
            if params["sub"] in (True, "True", "T", "true", 1): #Move str->bool conversion to ll.parse_argument
                ll.submit_job(job_script)

def bwa(params, sample, outdir):
    """Writes lines to job submission script to map 1 sample using bwa"""
    index = "{0}/{1}".format(sample[0], params["contig_names"])
    clean_reads = "{0}{1}/*".format(params["paths"]["clean"], sample[1].split('_')[0])
    try:
        fwd, rev = glob.glob(clean_reads)[:2]
    except ValueError:
        print("Can't find sequence data at {0}. Edit params[paths][clean] & contig_names".format(clean_reads))
    return [
        "module load {0}\n".format(params["modules"]["bwa"]),
        "bwa index {0}\n".format(index),
        "bwa mem {0} {1} {2} > {3}/aln-pe.sam".format(index, fwd, rev, outdir)
        # Need to add samtools view -> BAM, rm SAM, sort + index BAM, to mapping or add a mapping cleanup function to do this

           ]

def bowtie(params, sample, outdir):
    """Writes map script with bowtie, not finished as of 7/4/18,
    The structure below is based on Meren's call here:
    http://merenlab.org/2015/06/23/comparing-different-mapping-software/
    I'm not 100% sure how Meren defines things, so I didn't implement yet."""
    index = sample[0]+"contig.fa"
    pass
    #return [
#            "module load {}\n".format(params["modules"]["samtools"],
#            "module load {}\n".format(params["modules"]["bowtie"],
#            "bowtie2-build {0} {1}\n".format(contigs, samplename),
#            "bowtie2 -x {0} -f $query -S {}.bt2.sam\n".format(x=ref, query, sample),
#            "samtools view -bS {0}.bt2.sam > {0}.bt2.bam\n".format(sample),
#            "rm {0}.bt2.sam\n".format(sample)
#            ]

if __name__ == "__main__":
    main()
