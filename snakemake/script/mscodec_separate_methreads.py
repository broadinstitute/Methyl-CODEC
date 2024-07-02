#!/usr/bin/env python

import argparse
import logging
import os
import sys
import pysam

logger = logging.getLogger("{}".format(__file__))

def get_arguments():

    parser = argparse.ArgumentParser(prog="foo", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bam", type=str, help="input bam file")
    parser.add_argument("methyl", type=str, help="methyl reads bam output file")
    parser.add_argument("regular", type=str, help="regular reads bam output file")
    args = parser.parse_args()
    return args

def absolute_path(path):
    """convert relative path to absolute path"""
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(os.getcwd(), path)


def filter_reads_by_tag(input_bam, methyl_bam, regular_bam, tag_name):
    # Open the input BAM file
    with pysam.AlignmentFile(input_bam, 'rb') as infile:
        # Create an output BAM file
        with pysam.AlignmentFile(methyl_bam, 'wb', header=infile.header) as outfile, pysam.AlignmentFile(regular_bam, 'wb', header=infile.header) as outfile2:
            # Iterate through reads
            for read in infile:
                # Check if the specified tag is present in the read
                if read.is_duplicate:
                    continue
                #make read single end
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
                #output read
                if read.has_tag(tag_name):
                    outfile.write(read)
                else:
                    outfile2.write(read)
        # Write the read to the output BAM file



def process(options):
  # Example usage
  tag_to_filter = 'XM'  # Replace with your desired tag
  
  filter_reads_by_tag(options.bam, options.methyl, options.regular, tag_to_filter)

if __name__ == '__main__':
    sys.exit(process(get_arguments()))

