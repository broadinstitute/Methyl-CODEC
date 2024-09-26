#!/usr/bin/env python
import sys
import pysam
import argparse
from collections import defaultdict


def get_arguments():
    parser = argparse.ArgumentParser(prog="foo", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sourcebam", type=str, help="bam that are ccompared to")
    parser.add_argument("targetbam", type=str, help="bam that compares")
    parser.add_argument("-p", "--pairend",default=False, action='store_true', help="if target bam is paired-end data")
    parser.add_argument('-w', "--wiggle", type=int, default=150, help="max distance between alignments")
    parser.add_argument('-m', "--min_mapq", type=int, default=20, help="min mapping quality")
    parser.add_argument('-s', "--max_softclip", type=int, default=30, help="max softclip length")
    parser.add_argument('-v', "--verobse", action="store_true", help="print discordant reads")
    parser.add_argument('-o', "--outprefix", type=str, default="NA", help="print discordant reads")
    args = parser.parse_args()
    return args

def clip5_len(read):
    return 0 if read.cigartuples[0][0] != 4 else read.cigartuples[0][1]
def clip3_len(read):
    return 0 if read.cigartuples[-1][0] != 4 else read.cigartuples[-1][1]
def _get_unclipped_read_ends(read):
    """
    Accepts read and returns unclipped ends

    Returns
    -------
    0-based start, 1-based end coordinate of read
    """
    # Check if soft clipping (softclipping = 4)
    l_adj = clip5_len(read)
    r_adj = clip3_len(read)
    unclipped_start = read.reference_start - l_adj
    unclipped_end = read.reference_end + r_adj
    return unclipped_start, unclipped_end

def get_softclip_length(cigar_string):
    """
    Parse a CIGAR string and return the total length of soft clips (S operations).

    Args:
        cigar_string (str): CIGAR string (e.g., "5S90M5S")

    Returns:
        int: Total length of soft clips in the CIGAR string.
    """
    import re

    # Regular expression to match CIGAR operations
    cigar_pattern = re.compile(r'(\d+)([MIDNSHP=X])')

    total_softclip = 0

    # Find all CIGAR operations in the string
    for length, operation in cigar_pattern.findall(cigar_string):
        if operation == 'S':
            total_softclip += int(length)

    return total_softclip

def read_signle_generator(bam, **kwargs):
    for read in bam:
        if not _keep_read(read):
            continue
        yield read

def read_pair_generator(bam, **kwargs):
    """ Generate read pairs in a BAM file

    param: bam (pysam.AlignmentFile) - bam object to iterate over

    Note: cannot use bam_iterator's read pair generator as BAM is not indexed.

    yields tuple(read1, read2)
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam:
        if not _keep_read(read):
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
def _keep_read(read):
    return (
        #read.is_paired
        #and not read.is_unmapped
        #and not read.mate_is_unmapped
        not read.is_secondary
        and not read.is_supplementary
    )
def is_proper_orientation(read1, read2):
    if read1.is_unmapped or read2.is_unmapped:
        return False
    if read1.reference_name != read2.reference_name:
        return False
    if read1.is_reverse == read2.is_reverse:
        return False
    if read1.is_reverse and read1.reference_start < read2.reference_start:
        return False
    if read2.is_reverse and read2.reference_start < read1.reference_start:
        return False
    return True
def get_unclipped_fragment_ends(read1, read2):
    """
    Accepts a read pair and returns the fragment unclipped ends
    Useful when tracking down families after Fgbio

    RUOLIN: To be consistent with FGBIO and more reasonable,  for proper oriented pairs the fragment start is
    the leftmost unclipped site of the plus strand read and fragment end is the rightmost unclipped site
    of the reverse strand read.

    Returns:
    0-based start, 1-based end coordinate of fragment
    """
    if is_proper_orientation(read1, read2):
        start1, end1 = _get_unclipped_read_ends(read1)
        start2, end2 = _get_unclipped_read_ends(read2)
        if read1.is_reverse:
            fragment_start = start2
            fragment_end = end1
        else:
            fragment_start = start1
            fragment_end = end2
        return fragment_start, fragment_end
    else:
        raise Exception("UNDEFINED")
def get_insert_size(read1, read2):
    start, end = get_unclipped_fragment_ends(read1,read2)
    return end - start

def get_forward_query_alignment_end(read):
    if read.is_reverse:
        return read.query_length - read.query_alignment_start
    else:
        return read.query_alignment_end

def process(opt):
    discordant = 0
    not_aligned = 0
    concordant = 0
    pstrand_fail_mapq = 0
    cstrand_fail_mapq = 0
    bismark_fail_mapq = 0
    both_fail_mapq = 0
    total = 0
    #total, highconf = parse_linker_trim_log(opt.trimlog)
    iter_tr = iter(read_signle_generator(pysam.AlignmentFile(opt.targetbam, ignore_truncation=True)))
    iter_sr = iter(read_signle_generator(pysam.AlignmentFile(opt.sourcebam, ignore_truncation=True)))
    with open(f'{opt.outprefix}.discordant.alignment.txt', 'w') as f1, open(f'{opt.outprefix}.concordant.alignment.txt', 'w') as f2:
        try:
            while True:
                try:
                    ## assuming the two bams are generated according to the order in the fastq. 
                    ## iterate the smaller bam first. Asuming the reads in the smaller bam is a subset of the reads in the bigger bam
                    tr = next(iter_tr)
                    if opt.pairend:
                        tr = next(iter_tr)
                except StopIteration:
                    break
                while True:
                    try:
                        sr = next(iter_sr)
                        if sr.is_duplicate:
                            continue
                        total += 1
                    except StopIteration:
                        break
                    strand = "-" if sr.is_reverse else "+"
                    mscodec_converted_strand_mapq = sr.get_tag("MQ")
                    mate_cigar = sr.get_tag("MC")
                    # Parse the CIGAR string
                    cnvt_softclip_len = get_softclip_length(mate_cigar)
                    if sr.query_name != tr.query_name:
                        # if sr.mapping_quality < opt.min_mapq:
                        #     pstrand_fail_mapq += 1
                        # elif mscodec_converted_strand_mapq < opt.min_mapq:
                        #     cstrand_fail_mapq += 1
                        # else:
                        not_aligned +=1
                        print("not-aligned", sr.reference_name, sr.reference_start, sr.reference_end, sr.query_name, sr.mapping_quality, strand, sep="\t", file=f1)
                    else:
                        # #print("same", sr.query_name, tr.query_name)
                        # if sr.mapping_quality < opt.min_mapq:
                        #     pstrand_fail_mapq += 1
                        # elif mscodec_converted_strand_mapq < opt.min_mapq:
                        #     cstrand_fail_mapq += 1
                        # else:
                        #     if abs(sr.reference_start - tr.reference_start) > opt.wiggle or sr.reference_name != tr.reference_name:
                        #         discordant += 1
                        #         print("discordant", sr.reference_name, sr.reference_start, sr.reference_end, sr.query_name, sr.mapping_quality, strand, sep="\t", file=f1)
                        #     else:
                        #         concordant += 1
                        #         print(sr.reference_name, sr.reference_start, sr.reference_end, sr.query_name, sr.mapping_quality, strand, sep="\t", file=f2)
                        if abs(sr.reference_start - tr.reference_start) > opt.wiggle or sr.reference_name != tr.reference_name:
                           if tr.mapping_quality < 40 and (sr.mapping_quality < opt.min_mapq or mscodec_converted_strand_mapq < opt.min_mapq or cnvt_softclip_len > opt.max_softclip):
                               both_fail_mapq += 1
                               print("fail_both_mapq", sr.reference_name, sr.reference_start, sr.reference_end, sr.query_name, sr.mapping_quality, strand, sep="\t", file=f1)
                           elif tr.mapping_quality < 40 and sr.mapping_quality >= opt.min_mapq and mscodec_converted_strand_mapq >= opt.min_mapq and cnvt_softclip_len <= opt.max_softclip:
                               bismark_fail_mapq += 1
                               print("fail_bismark_mapq", sr.reference_name, sr.reference_start, sr.reference_end, sr.query_name, sr.mapping_quality, strand, sep="\t", file=f1)
                           elif sr.mapping_quality < opt.min_mapq and tr.mapping_quality >= 40:
                               pstrand_fail_mapq += 1
                               print("fail_mscodec_mapq", sr.reference_name, sr.reference_start, sr.reference_end, sr.query_name, sr.mapping_quality, strand, sep="\t", file=f1)
                           elif (mscodec_converted_strand_mapq < opt.min_mapq or cnvt_softclip_len > opt.max_softclip) and tr.mapping_quality >= 40:
                               cstrand_fail_mapq += 1
                               print("fail_mscodec_mapq", sr.reference_name, sr.reference_start, sr.reference_end, sr.query_name, sr.mapping_quality, strand, sep="\t", file=f1)
                           else:
                               discordant += 1
                               print("discordant", sr.reference_name, sr.reference_start, sr.reference_end, sr.query_name, sr.mapping_quality, strand, sep="\t", file=f1)
                        else:
                            concordant += 1
                            print(sr.reference_name, sr.reference_start, sr.reference_end, sr.query_name, sr.mapping_quality, strand, sep="\t", file=f2)
                        break
        except Exception as error:
            print("An exception occurred:", error, file=sys.stderr)  # An exception occurred: division by zero
        finally:
            print("total", "bismark_not_aligned", "fail_mapq_both" ,"fail_mapq_bismark", "fail_mapq_mscodec",  "concordant", "discordant",  "concordant ratio", sep="\t")
            print(total, not_aligned, both_fail_mapq, bismark_fail_mapq,  pstrand_fail_mapq + cstrand_fail_mapq, concordant, discordant, concordant/total, sep="\t")

if __name__ == '__main__':
    sys.exit(process(get_arguments()))