#!/usr/bin/env python

import argparse
import logging
import os
import sys
import gzip
from itertools import groupby
from operator import itemgetter



logger = logging.getLogger("{}".format(__file__))

def get_arguments():

    parser = argparse.ArgumentParser(prog="calculate methylation level per fragment", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", type=str, help="input file from Bismark Extractor")
    parser.add_argument("output", type=str, help="output file")
    args = parser.parse_args()
    return args

def absolute_path(path):
    """convert relative path to absolute path"""
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(os.getcwd(), path)
def process_group(key, group):
    total = 0
    z_count = 0
    chrom = None
    left = sys.maxsize
    right = 0

    for item in group:
        total += 1
        if not chrom:
            chrom = item[2]
        left = min(left, int(item[3])-1)
        right = max(right, int(item[3]))
        if item[-1] == 'Z':
            z_count += 1

    # Calculate fraction of 'Z'
    return chrom, left, right, key, total, z_count, z_count / total if total > 0 else 0

def process(options):
    split_data = None
    with gzip.open(options.input, 'rt') as f: # Open the gzipped file
        with open(options.output, 'w') as fo:
            next(f)
            split_data = (line.split() for line in f)
            for key, group in groupby(split_data, key=itemgetter(0)):
                chrom, left, right, key,total, Z_count, fraction = process_group(key, group)
                # Print out the fractions
                print(f"{chrom}\t{left}\t{right}\t{key}\t{total}\t{Z_count}\t{fraction:.4f}", file=fo)

if __name__ == '__main__':
    sys.exit(process(get_arguments()))
# Convert input data into a pandas DataFrame
