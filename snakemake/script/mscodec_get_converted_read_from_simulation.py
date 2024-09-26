#!/usr/bin/env python
import sys
import pysam
import argparse
from collections import defaultdict


def get_arguments():
    parser = argparse.ArgumentParser(prog="foo", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputbam", type=str, help="mscodec aligned sortbyname read")
    parser.add_argument('outbam',  type=str,  help="output converted read bam")
    parser.add_argument('-m', "--min_mapq",  type=int, default=0,  help="min mapq to output")
    args = parser.parse_args()
    return args


def process(opt):

    #total, highconf = parse_linker_trim_log(opt.trimlog)
    bismark_alignment_bam = pysam.AlignmentFile(opt.inputbam, ignore_truncation=True)
    with pysam.AlignmentFile(opt.outbam, 'wb', header=bismark_alignment_bam.header) as outfile:
        try:
            with pysam.AlignmentFile(opt.inputbam, 'rb') as infile:
                for read in infile:
                    if read.has_tag("XR"):
                        xr = read.get_tag("XR")
                        if xr == "GA" and read.mapping_quality >= opt.min_mapq:
                            outfile.write(read)
        except Exception as error:
            print("An exception occurred:", error, file=sys.stderr)  # An exception occurred: division by zero
        finally:
            outfile.close()

if __name__ == '__main__':
    sys.exit(process(get_arguments()))