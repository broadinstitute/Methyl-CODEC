#!/usr/bin/env python

import argparse
import logging
import os
import sys

logger = logging.getLogger("{}".format(__file__))

def get_arguments():

    parser = argparse.ArgumentParser(prog="foo", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bam", type=str, help="input bam file")
    parser.add_argument("output", type=str, help="output bam file")
    args = parser.parse_args()
    return args

def absolute_path(path):
    """convert relative path to absolute path"""
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(os.getcwd(), path)

import pysam

def filter_reads_by_tag(input_bam, output_bam, tag_name):
    # Open the input BAM file
    with pysam.AlignmentFile(input_bam, 'rb') as infile:
        # Create an output BAM file
        with pysam.AlignmentFile(output_bam, 'wb', header=infile.header) as outfile:
            # Iterate through reads
            for read in infile:
                # Check if the specified tag is present in the read
                if read.has_tag(tag_name):
                    # Check if the tag value matches the desired value
                    if read.has_tag(tag_name) and not read.is_duplicate:
                        # Write the read to the output BAM file
                        if read.is_reverse:
                            read.flag = 16
                        else:
                            read.flag = 0
                        if read.is_secondary:
                            read.flag += 256
                        if read.is_supplementary:
                            read.flag += 2048
                        read.next_reference_id = -1
                        read.next_reference_start = -1
                        read.template_length = 0
                        outfile.write(read)



def process(options):
  # Example usage
  tag_to_filter = 'XM'  # Replace with your desired tag
  
  filter_reads_by_tag(options.bam, options.output, tag_to_filter)

if __name__ == '__main__':
    sys.exit(process(get_arguments()))

