"""
Annotation processing
Written by JG 8/6/18

Follows outline here: https://metagenomics-workshop.readthedocs.io/en/latest/annotation/functional_annotation.html
After this, a bash script calling minpath is run on EC annotations.

Command line arguments:
    -d,--directory      Directory containing bins.
"""
__author__ = "Jim Griffin"
__author_email__ = "jimbogrif@gmail.com"

# Imports:
import argparse
import os
import sys
import pandas as pd


def main():
    counts_df = pd.DataFrame(columns=["Annotations", "Enzymes", "COGs"])
    args = parse_arguments()

    for root, dirs, files in os.walk(args.directory):
        for file in files:
            if '.gff' in file:
                # Initialize new file locations based on this root directory:
                bin_name = root.rsplit('/', 1)[-1]
                file_in = os.path.join(root, file)
                ec_out = os.path.join(root, 'ECs.txt')
                cog_out = os.path.join(root, 'COGs.txt')

                # Add new row to dataframe and initialize count data:
                counts_df.loc[bin_name] = [0, 0, 0]
                enzymes, COGs, total_annotations = 0, 0, 0
                with open(file_in, 'r') as f, open(ec_out, 'w') as ec_file, open(cog_out, 'w') as cog_file:
                    for line in f:
                        total_annotations += 1
                        if 'eC_number=' in line:
                            ID, EC = line.split()[8].split(';')[:2]
                            ID = ID.replace('ID=', '')
                            EC = EC.replace('eC_number=', '')
                            ec_file.write("{0}\t{1}\n".format(ID, EC))
                            enzymes += 1
                        if 'COG' in line:
                            COG_annotation = ' '.join(line.split()[8:])
                            ID = COG_annotation.split(';')[0].replace('ID=', '')
                            COG_details = COG_annotation.split(';', 1)[1]
                            cog_file.write("{0}\t{1}\n".format(ID, COG_details))
                            COGs += 1

                    counts_df.loc[bin_name] = [total_annotations, enzymes, COGs]
    counts_df.to_csv(args.directory+'annotation_counts.txt', sep='\t')

def parse_arguments():
    """Parse and return command line arguments as an argparse args object"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", "-d", type=str,
                        help="Directory containing per-bin annotation subdirectories.")
    try:
        result = parser.parse_args()
        return result
    except Exception as e:
        parser.print_help()
        print(e)
        sys.exit(0)

if __name__ == "__main__":
    main()
