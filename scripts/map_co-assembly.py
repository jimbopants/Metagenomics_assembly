"""
Map to the file specified in config.yaml as the co-assembly.
Similar options & usage to map_reads.py otherwise

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
    params, args = ll.load_params_and_input("map")
    os.makedirs(params["path_out"]["co_assembly_map"], exist_ok=True)
    samples_to_map = ll.paths_and_samples(params["paths"]["reads_to_map"])

    for sample in samples_to_map:
        outdir = ll.make_dir("{0}{1}".format(params["path_out"]["co_assembly_map"], sample[1]))
        job_script = "{0}/{1}.{2}.sh".format(outdir.rsplit('/',1)[0], sample[1], "BWA_co_assembly")

        # Write job script:
        with open(job_script, "w") as file_out:
            msub_lines = ll.fill_header(job_script, params)
            msub_lines += BWA_co(params, sample, outdir, params["paths"]["co_index"])
            for line in msub_lines:
                file_out.write(line)
        if params["sub"] in (True, "True", "T", "true", 1): #Move str->bool conversion to ll.parse_argument
            ll.submit_job(job_script)


def BWA_co(params, sample, outdir, co_index):
    """All individual -> single co-assembly BWA mapping and convert to binary"""
    clean_reads = "{0}{1}/*".format(params["paths"]["clean"], sample[1].split('_')[0])
    try:
        fwd, rev = glob.glob(clean_reads)[:2]
    except ValueError:
        print("Can't find sequence data at {0}. Edit params[paths][clean] & contig_names".format(clean_reads))
    return [
        "module load {0}\n".format(params["modules"]["bwa"]),
        "bwa index {0}\n".format(co_index),
        "bwa mem {0} {1} {2} > {3}/aln-pe.sam\n".format(co_index, fwd, rev, outdir),
        "samtools view -bS {0}/aln-pe.sam > {0}/aln-pe.bam\n".format(outdir),
        "rm {0}/aln-pe.sam\n".format(outdir)
            ]

if __name__ == "__main__":
    main()
