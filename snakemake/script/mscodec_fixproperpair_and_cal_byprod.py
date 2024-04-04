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
        qname = read.query_name
        if qname not in read_dict:
            if len(read_dict) > 0:
                for _, (read1, read2) in read_dict.items():
                    if read1:
                        yield read1, None
                    if read2:
                        yield read2, None
                read_dict.clear()
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
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
    parser.add_argument("rescued_paired_bam", type=str, help="methylation bam file")
    parser.add_argument("single_strand_bam", type=str, help="output bam file")
    args = parser.parse_args()
    return args



def FixProperPairFlag(inbam, single_strand_bam, rescued_paired_bam, im_dist_cutoff = 5_000, adap_v2 = False):
    unpaired_bam = pysam.AlignmentFile(inbam, "rb")
    single_strand_bam_writer = pysam.AlignmentFile(single_strand_bam, "wb", template=unpaired_bam)
    rescued_paired_bam_writer = pysam.AlignmentFile(rescued_paired_bam, "wb", template=unpaired_bam)
    total_frag = 0
    correct_frag = 0
    unmapped = 0
    intermol = 0
    double_ligation = 0
    for read1, read2 in read_pair_generator(unpaired_bam):
        total_frag += 1
        if not read2:
            if read1.infer_query_length() != read1.query_length:
                print(read1.query_name, " has different query length and infer query length", file=sys.stderr)
                continue
            if read1.is_reverse:
                read1.flag = 0x10
            else:
                read1.flag = 0
            single_strand_bam_writer.write(read1)
            if read1.query_length < 15:
                double_ligation += 1
            else:
                unmapped += 1
            continue
        if read1.infer_query_length() != read1.query_length or read2.infer_query_length() != read2.query_length:
            print(read1.query_name, " has different query length and infer query length", file=sys.stderr)
            continue
        read1.flag |= 1
        read2.flag |= 1
        if read1.query_length < 15 or read2.query_length < 15:
            double_ligation += 1
            single_strand_bam_writer.write(read1)
            single_strand_bam_writer.write(read2)
            continue
        if read1.is_unmapped or read2.is_unmapped:
            unmapped += 1
            single_strand_bam_writer.write(read1)
            single_strand_bam_writer.write(read2)
            continue

        is_correct = False
        if read1.reference_name == read2.reference_name:
            if abs(read1.template_length) < im_dist_cutoff:
                if not read1.is_reverse and read2.is_reverse:
                   if read1.reference_start <= read2.reference_start:
                       read1.flag |= 2
                       read2.flag |= 2
                       correct_frag += 1
                       is_correct = True
                   else:
                       intermol += 1
                elif not read2.is_reverse and read1.is_reverse:
                    if read2.reference_start <= read1.reference_start:
                        read1.flag |= 2
                        read2.flag |= 2
                        correct_frag += 1
                        is_correct = True
                    else:
                        intermol += 1
                else:
                    intermol+=1

            else:
                intermol += 1
        else:
            intermol += 1
        if is_correct:
            if read1.has_tag('XM') != read2.has_tag('XM'):
                rescued_paired_bam_writer.write(read1)
                rescued_paired_bam_writer.write(read2)
            else:
                print("concordant read pair ", read1.query_name, " both have or do not have XM tags", file=sys.stderr)
                single_strand_bam_writer.write(read1)
                single_strand_bam_writer.write(read2)
        else:
            single_strand_bam_writer.write(read1)
            single_strand_bam_writer.write(read2)

    return total_frag, correct_frag, intermol, unmapped,double_ligation

def process(opts):
    total_frag, correct_frag, intermol, unmapped,double_ligation = FixProperPairFlag(sys.stdin, opts.single_strand_bam, opts.rescued_paired_bam)
    print("total_frag\tcorrect_frag\tintermol\tunmapped\tdouble_ligation ")
    print(total_frag, correct_frag, intermol, unmapped, double_ligation, sep="\t")

if __name__ == '__main__':
    sys.exit(process(get_arguments()))
