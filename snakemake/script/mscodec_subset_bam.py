#!/usr/bin/env python

import argparse
import logging
import os
import sys
import pysam

logger = logging.getLogger("{}".format(__file__))

def get_arguments():

    parser = argparse.ArgumentParser(prog="foo", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bam_file", type=str, help="input bam file")
    parser.add_argument("readname_file", type=str, help="input readname file")
    parser.add_argument("bam_output", type=str, help="output bam")
    parser.add_argument("--readname_index", "-n", default = 0, type=int, help="column idx for read name")
    parser.add_argument("--condition_index", "-i", default = None, type=int, nargs="+", help="column idx which needs to satisfy condition")
    parser.add_argument("--condition", "-c", default = None, type=str, nargs="+", help="condition to satisfy")
    args = parser.parse_args()
    return args

def absolute_path(path):
    """convert relative path to absolute path"""
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(os.getcwd(), path)

def check_if_number(input_arg):
    try:
        float(input_arg)
        return True
    except ValueError:
        return False


def subset_bam_by_read_names(options):
    # Read the read names into a set for quick lookup
    read_names = set()
    with open(options.readname_file, 'r') as f:
        if options.condition_index is not None and options.condition is not None:
            print(options.condition_index, options.condition)
            for line in f:
                columns = line.split('\t')  # Split the line into columns

                # Ensure there are enough columns in the line
                if len(columns) > 3:
                    # Iterate through each condition index
                    for ii, cond in zip(options.condition_index, options.condition):
                        value_to_check = columns[ii].strip()  # The value in the relevant column

                        # Check condition based on whether cond is numeric or not
                        if check_if_number(cond):
                            condition_met = float(value_to_check) == float(cond)
                        else:
                            condition_met = value_to_check == cond

                        # If the condition is met, add the read name to the set
                        if condition_met:
                            read_name = columns[options.readname_index].strip()
                            read_names.add(read_name)

            # Print the size of read_names set after processing each line
            print("size of read_names:", len(read_names))
        else:
            read_names = set(line.split('\t')[options.readname_index].strip() for line in f if len(line.split('\t')) > 3)

    #print(read_names)
    # Open the input BAM file
    bam = pysam.AlignmentFile(options.bam_file, "rb")
    outbam = pysam.AlignmentFile(options.bam_output, "wb", header=bam.header)

    # Iterate through each read in the BAM file
    for read in bam:
        # Check if the read name is in the set of read names
        if read.query_name in read_names:
            # Write the read to the output BAM file
            outbam.write(read)

    # Close the input BAM file
    bam.close()
    outbam.close()

def process(options):
    # Example usage

    subset_bam_by_read_names(options)

if __name__ == '__main__':
    sys.exit(process(get_arguments()))