#!/usr/bin/env python

import argparse
import logging
import os
import sys
import pysam

logger = logging.getLogger("{}".format(__file__))

def get_arguments():

    parser = argparse.ArgumentParser(prog="filter converted strand for methylation extraction", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_bam", type=str, help="input bam file")
    parser.add_argument("output_bam", type=str, help="output bam file")
    parser.add_argument("--pstrand_min_mapq", "-a", type=int, default=0, help="min mapping quality for protected strand")
    parser.add_argument("--cstrand_min_mapq", "-b", type=int, default=0, help="min mapping quality for converted strand")
    parser.add_argument("--cstrand_max_softclip", "-s", type=int, default=sys.maxsize, help="max softclip len for converted strand")
    args = parser.parse_args()
    return args

def absolute_path(path):
    """convert relative path to absolute path"""
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(os.getcwd(), path)


def get_soft_clip_len(read):
    cigar_tuples = read.cigartuples

    soft_clipping_length = 0
    if cigar_tuples:
        # Initialize soft clipping length

        # CIGAR code 4 corresponds to soft clipping
        for operation, length in cigar_tuples:
            if operation == 4:  # Soft clipping
                soft_clipping_length += length
    return soft_clipping_length
def process(options):
    with pysam.AlignmentFile(options.input_bam, "rb") as infile:
        # Create an output BAM file with the same header as the input
        with pysam.AlignmentFile(options.output_bam, "wb", header=infile.header) as outfile:
            # Iterate over each read in the input BAM file
            for read in infile:
                # Check if the read and its mate satisfy the mapping quality thresholds
                pstrand_mapq = read.get_tag("MQ")
                if read.mapping_quality >= options.cstrand_min_mapq and \
                    pstrand_mapq >= options.pstrand_min_mapq and \
                    get_soft_clip_len(read) < options.cstrand_max_softclip:
                    outfile.write(read)


if __name__ == '__main__':
    sys.exit(process(get_arguments()))