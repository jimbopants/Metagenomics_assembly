"""
Cleanup intermediate files:
    This deletes unused files created by programs like idba and megahit that
    aren't used and take up a lot of space.
    Check your filenames to make sure they match the directory structure below:

Directory structures:
    assembled /
        megahit /
            S1 /
                final_contigs.fa, log, [intermediate_contigs]
        idba /
            S1 /
                contigs, scaffolds, log, [OTHERFILES]

**Files in [] will be deleted.

Arguments (Mandatory):
  --directory DIRECTORY, -d DIRECTORY
                        Directory with samples to cleanup
  --methods METHOD, -m METHOD
                        Method that was used to generate the files to be
                        cleaned up
Update date: 7/2/18
"""
__author__ = "Jim Griffin"
__author_email__ = "jimbogrif@gmail.com"

# Imports:
import shutil
import os
import shared_utilities as ll

def main():
    # Parse arguments
    params, args = ll.load_params_and_input(None)

    # Map method argument to cleanup function
    run_cleanup = {"megahit" : megahit, "idba" : idba}

    # Get all samples to cleanup:
    cleanup_folders = ll.paths_and_samples(args.directory)

    # Run cleanup scripts:
    for sample in cleanup_folders:
        run_cleanup[args.methods](sample[0].rstrip('/'))

def idba(path_in):
    """ Removes IDBA intermediates in place."""
    files_to_keep = ["contig.fa", "scaffold.fa", "log"]
    for file in os.scandir(path_in):
        if file.name not in files_to_keep:
            os.remove(file.path)
    return

def megahit(path_in):
    """ Removes megahit intermediates in place. Renames contigs to contig.fa"""
    shutil.rmtree("{0}/{1}".format(path_in, "intermediate_contigs/")) #
    os.rename("{0}/final.contigs.fa".format(path_in), "{0}/contig.fa".format(path_in))
    os.listdir(path_in)
    return

if __name__ == "__main__":
    main()
