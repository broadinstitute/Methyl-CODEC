#!/usr/bin/env python

import argparse
import logging
import os
import sys
from pybedtools import BedTool

logger = logging.getLogger("{}".format(__file__))

def get_arguments():

    parser = argparse.ArgumentParser(prog="foo", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("interval_bed", type=str, help="the bed file which is to be queried")
    parser.add_argument("methylation_levels", type=str, help="methylation levels bed")
    args = parser.parse_args()
    return args

def absolute_path(path):
    """convert relative path to absolute path"""
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(os.getcwd(), path)

def group_overlapping_entries(bed1_path, bed2_path):
    # Load BED files using pybedtools
    bed1 = BedTool(bed1_path)
    bed2 = BedTool(bed2_path)

    # Group entries from bed2 based on overlap with bed1
    grouped_bed2 = bed1.intersect(bed2, wa=True, wb=True, loj=True)

    # Create a dictionary to store grouped entries
    grouped_entries = {}
    # Process the grouped_bed2 and organize entries in the dictionary
    for entry in grouped_bed2:
        bed1_entry = str(entry).strip('\n').split('\t')[:len(bed1[0].fields)]
        bed2_entry = str(entry).strip('\n').split('\t')[len(bed1[0].fields):]

        # Check if bed2 entry overlaps with bed1 entry
        if bed1_entry != ['.'] * len(bed1[0]):
            key = tuple(bed1_entry)
            value = tuple(bed2_entry)

            if key in grouped_entries:
                grouped_entries[key].append(value)
            else:
                grouped_entries[key] = [value]

    return grouped_entries

def process(options):
    result = group_overlapping_entries(options.interval_bed, options.methylation_levels)
    #header="\t".join(['chrom', 'start', 'end', 'name', 'length', 'cpgNum', 'gcNum', 'perCpG', 'perGC', 'obsExp', 'met_count', 'unmet_count', 'total_count', 'per_met'])
    #print(header)
    for key, values in result.items():
        met_count = 0
        unmet_count = 0
        for value in values:
            try:
                met_count += int(value[-2])
                unmet_count += int(value[-1])
            except ValueError as ve:
                print(f"Error: {ve}", file=sys.stderr)
        try:
            perc_met = met_count/(met_count + unmet_count)
        except ZeroDivisionError:
            perc_met = "NA"
        keystr = "\t".join(key)
        print(f"{keystr}\t{met_count}\t{unmet_count}\t{met_count+unmet_count}\t{perc_met}")

if __name__ == '__main__':
    sys.exit(process(get_arguments()))

