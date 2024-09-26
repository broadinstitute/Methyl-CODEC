#!/usr/bin/env python

import argparse
import logging
import os
import sys

logger = logging.getLogger("{}".format(__file__))

def get_arguments():

    parser = argparse.ArgumentParser(prog="foo", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inbam", type=str, help="input bam file")
    parser.add_argument("outbam", type=str, help="output bam file")
    args = parser.parse_args()
    return args

def absolute_path(path):
    """convert relative path to absolute path"""
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(os.getcwd(), path)

import pysam

def convert_c_to_t_in_bam(input_bam, output_bam):
    # Open the input BAM file for reading
    with pysam.AlignmentFile(input_bam, "rb") as in_bam:
        # Open the output BAM file for writing
        with pysam.AlignmentFile(output_bam, "wb", header=in_bam.header) as out_bam:
            for read in in_bam:
                if read.is_duplicate:
                    continue
                # Convert the sequence from bytes to a mutable list
                if not read.has_tag("XM"):
                    sequence = list(read.query_sequence)
                    
                    # Iterate through the sequence and change all 'C' to 'T'
                    if read.is_reverse:
                        for i, base in enumerate(sequence):
                            if base == 'G':
                                sequence[i] = 'A'
                    else:
                        for i, base in enumerate(sequence):
                            if base == 'C':
                                sequence[i] = 'T'
                    
                    # Update the read with the new sequence
                    read.query_sequence = ''.join(sequence)
                
                # Write the modified read to the output BAM file
                out_bam.write(read)



def process(options):
    convert_c_to_t_in_bam(options.inbam, options.outbam)

if __name__ == '__main__':
    sys.exit(process(get_arguments()))

