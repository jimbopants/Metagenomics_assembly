"""
Works similarily to map_reads, assumes you have a directory with
per sample directory of contigs to annotate.
"""

# Imports:
import glob
import os
import shared_utilities as ll


def main():
    """Anotate sample contigs with prokka using default settings + metagenome option"""
    params, args = ll.load_params_and_input(None)
    os.makedirs(params["paths"]["annotations"], exist_ok=True)

    # Annotate the same assemblies we mapped to earlier in the pipeline:
    samples_to_annotate = ll.paths_and_samples(params["paths"]["reads_to_map"])

    for sample in samples_to_annotate:
        outdir = ll.make_dir("{0}{1}".format(params["paths"]["annotations"], sample[1]))
        job_script = "{0}/{1}.prokka.sh".format(outdir.rsplit('/',1)[0], sample[1])

        # Write job script:
        with open(job_script, "w") as file_out:
            msub_lines = ll.fill_header(job_script, params)
            msub_lines += prokka(params, sample, outdir)
            for line in msub_lines:
                file_out.write(line)
        if params["sub"] in (True, "True", "T", "true", 1): #TODO: Move str->bool conversion to ll.parse_argument
            ll.submit_job(job_script)

def prokka(params, sample, outdir):
    """ Reasonable looking metagenome settings with prokka"""
    full_sample_path = "{0}/contig.fa".format(sample[0])
    lines = [ "source activate {0}\n".format(params["python_envs"]["prokka"]),
              "prokka {0} --outdir {1}/ --metagenome\
               --cpus {2} --force".format(full_sample_path,
                                          outdir,
                                          params["num_threads"])
            ]
    return lines


if __name__ == "__main__":
    main()
