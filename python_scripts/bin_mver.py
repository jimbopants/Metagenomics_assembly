"""
Move high quality bins into a new annotations folder

"""
__author__ = "Jim Griffin"
__author_email__ = "jimbogrif@gmail.com"

import os
import sys
import shutil


def main():
    # Hard code for now cuz doing in 10 minutes
    root = '/projects/b1042/Wells/Jim/Wells01/binning/'
    bins = [root+'annotations/maxbin_annotations/',
            root + 'annotations/metabat_annotations/']
    bin_lists = [root +  'bin_subsets/HQ_bins.txt',
                 root + 'bin_subsets/HQ_bins.txt']

    dest = root + 'bin_subsets/'
    os.makedirs(dest+'maxbin', exist_ok=True)
    os.makedirs(dest+'metabat', exist_ok=True)

    # Check both bins against both bin list
    for bin_list in bin_lists:
        bins_to_move = read_bin_list(bin_list)
        move_bins(bins_to_move, bins[0], dest+'maxbin/')
        move_bins(bins_to_move, bins[1], dest+'metabat/')


def read_bin_list(file):
    bin_list = []
    with open(file, 'r') as f:
        for line in f:
            bin_list.append(line.strip())
    return bin_list

def move_bins(bin_list, bin_dir, dest):
    dirs = [name for name in os.listdir(bin_dir)
            if os.path.isdir(os.path.join(bin_dir, name))]
    found_bins = [bin_name for bin_name in dirs if bin_name in bin_list]
    for bins in found_bins:
        full_path = bin_dir+bins
        shutil.move(full_path, dest)

if __name__ == "__main__":
    main()
