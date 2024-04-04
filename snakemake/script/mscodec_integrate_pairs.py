#!/usr/bin/env python
import argparse
import logging
import sys
import os
import pysam
from collections import defaultdict
import re
#from Bio import pairwise2

logger = logging.getLogger("{}".format(__file__))

#def check_tandem_adpt(seq):
    #linker = "AGATCGGAAGAGCTTCATCATTAGATCCATTAATGTTACACTTCAACTCTTCACCCACATCAGATTAGTACCAGCTTCGAGGATCAACACGTCAGAGTCTAGCTGGTGATAGGAAGTGTAGGTAACATAGACGAAGTTATCAACAATGTGTAACTGACTTAACGCTCTTCCGATCT"
    #res = pairwise2.align.localms(seq, linker, 1, -4, -6,-2)

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            continue
        qname = re.sub(r'/[12]$', "", read.query_name)
        if qname not in read_dict:
            if len(read_dict) > 0:
                for _, (read1, read2) in read_dict.items():
                    if read1:
                        yield read1, None
                    if read2:
                        yield read2, None
                read_dict.clear()

            if read.is_read1 or "/1" in read.query_name:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1 or "/1" in read.query_name:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            read_dict.clear()
    for qname, (read1, read2) in read_dict.items():
        if read1:
            yield read1, None
        if read2:
            yield read2, None

def is_overlapped(read1, read2):
    if read1.is_unmapped or read2.is_unmapped:
        return False
    if read1.reference_name != read2.reference_name:
        return False
    if read1.reference_start < read2.reference_end and read2.reference_start < read1.reference_end:
        return True

def is_complete_overlapped_excluding_sclips(read1, read2):
    if not is_overlapped(read1, read2):
        return False
    if read1.reference_start != read2.reference_start:
        return False
    if read1.reference_end != read2.reference_end:
        return False
    return True

def overlap_len(read1, read2):
    if not is_overlapped(read1,read2):
        return 0
    else:
        return min(read1.reference_end, read2.reference_end) - max(read1.reference_start, read2.reference_start)

def overlap_span_ratio(read1, read2):
    ol = overlap_len(read1, read2)
    if ol == 0:
        return 0;
    else:
        span = max(read1.reference_end, read2.reference_end) - min(read1.reference_start, read2.reference_start)
        return ol/span

def get_arguments():

    parser = argparse.ArgumentParser(prog="fix flag for separately aligned paired end reads", formatter_class=argparse.RawDescriptionHelpFormatter)
    args = parser.parse_args()
    return args


def reset_flag_bismark_aln(read):
    XR = read.get_tag('XR')
    XG = read.get_tag('XG')
    if XR == "CT" and XG == "CT":
        read.flag &= ~(0x10)
    if XR == "GA" and XG == "CT":
        read.flag |= 0x10
    if XR == "CT" and XG == "GA":
        read.flag |= 0x10
    if XR == "GA" and XG == "GA":
        read.flag &= ~(0x10)

def merge_single_alignments(inbam, outbam, im_dist_cutoff = 5_000, adap_v2 = False):
    unpaired_bam = pysam.AlignmentFile(inbam, "rb")
    outbam_writer = pysam.AlignmentFile(outbam, "wb", template=unpaired_bam)
    total_frag = 0
#reaplace reada as read1 and readb as read2

    for read1, read2 in read_pair_generator(unpaired_bam):
        total_frag += 1
        if not read2:
            #read1 could be read1 or read2
            read1.flag |= 1
            if read1.is_read1 or "/1" in read1.query_name:
                read1.flag |= 0x40
            else:
                read1.flag |= 0x80
            read1.query_name = re.sub(r'/[12]$', "", read1.query_name)
            outbam_writer.write(read1)
            continue

        # add assert function so that reada is always read1 and readb is always read2
        assert read1.is_read1 or "/1" in read1.query_name
        assert read2.is_read2 or "/2" in read2.query_name


        read1.flag |= 1
        read2.flag |= 1
        read1.flag |= 0x40
        read2.flag |= 0x80
        if "/1" in read1.query_name or "/2" in read1.query_name:
            read1.query_name = re.sub(r'/[12]$', "", read1.query_name)
            reset_flag_bismark_aln(read1)

        if "/1" in read2.query_name or "/2" in read2.query_name:
            read2.query_name = re.sub(r'/[12]$', "", read2.query_name)
            reset_flag_bismark_aln(read2)

        outbam_writer.write(read1)
        outbam_writer.write(read2)

    return total_frag

def process(opts):
    merge_single_alignments(sys.stdin, sys.stdout)

if __name__ == '__main__':
    sys.exit(process(get_arguments()))
