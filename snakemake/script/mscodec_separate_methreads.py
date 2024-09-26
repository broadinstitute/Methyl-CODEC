#!/usr/bin/env python

import argparse
import logging
import os
import sys
import pysam
from collections import defaultdict

logger = logging.getLogger("{}".format(__file__))

def get_arguments():

    parser = argparse.ArgumentParser(prog="filter and separate strands", formatter_class=argparse.RawDescriptionHelpFormatter)
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

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

def filter_reads_by_tag_2(input_bam, methyl_bam, regular_bam, tag_name):
    # Open the input BAM file
    with pysam.AlignmentFile(input_bam, 'rb') as infile:
        # Create an output BAM file
        with pysam.AlignmentFile(methyl_bam, 'wb', header=infile.header) as outfile, pysam.AlignmentFile(regular_bam, 'wb', header=infile.header) as outfile2:
            # Iterate through reads

            for read1, read2 in read_pair_generator(input_bam):
                # Check if the specified tag is present in the read
                if read1.is_duplicate or read2.is_duplicate:
                    continue
                #make read single end
                assert(read1.is_reverse != read2.is_reverse)
                if read1.is_reverse:
                    read1.flag = 16
                    read2.flag = 0
                else:
                    read1.flag = 0
                    read2.flag = 16
                read1.next_reference_id = -1
                read1.next_reference_start = -1
                read1.template_length = 0
                read2.next_reference_id = -1
                read2.next_reference_start = -1
                read2.template_length = 0
                #output read1
                if read1.has_tag(tag_name):
                    outfile.write(read1)
                    outfile2.write(read2)
                else:
                    outfile2.write(read1)
                    outfile.write(read2)
        # Write the read to the output BAM file

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


def process(options):
  # Example usage
  tag_to_filter = 'XM'  # Replace with your desired tag
  
  filter_reads_by_tag(options.bam, options.methyl, options.regular, tag_to_filter)

if __name__ == '__main__':
    sys.exit(process(get_arguments()))

