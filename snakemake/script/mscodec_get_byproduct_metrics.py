#!/usr/bin/env python
import argparse
import sys
import pandas as pd
import os
from collections import defaultdict

def get_arguments():

    parser = argparse.ArgumentParser(prog="foo", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sampleid", type=str, help="sample id")
    parser.add_argument("trimlog", type=str, help="trim adapter log file")
    parser.add_argument("msalignlog", type=str, help="ms-align log file")
    args = parser.parse_args()
    return args
def process(opts):
    trimdf = pd.read_csv(opts.trimlog, sep=": ", header=None)
    tot = trimdf[trimdf[0] == "TOTAL"][1].values[0]
    adapter = trimdf[trimdf[0] == "LOST_BOTH"][1].values[0]
    double_ligation = trimdf[trimdf[0] == "LOST_READ1"][1].values[0] + trimdf[trimdf[0] == "LOST_READ2"][1].values[0]

    mslog = pd.read_csv(opts.msalignlog, sep="\t", header=None)

    correct=mslog[mslog[0] == "correct"][1].values[0]
    intermol = mslog[mslog[0] == 'intermolecular'][1].values[0]

    header = ["sample_id",
              ### CDS specific
              "pct_correct",
              "pct_double_ligation",
              "pct_adp_dimer",
              "pct_intermol",
              "pct_categorized",
              "n_correct",
              "n_double_ligation",
              "n_adp_dimer",
              "n_intermol",
              "n_categorized",
              "n_total_fastq",
              ]
    print("\t".join(header))
    print(
        opts.sampleid,
        correct/tot,
        double_ligation/tot,
        adapter/tot,
        intermol/tot,
        (correct+double_ligation+adapter+intermol)/tot,
        correct,
        double_ligation,
        adapter,
        intermol,
        (correct+double_ligation+adapter+intermol),
        tot, sep="\t")

if __name__ == '__main__':
    sys.exit(process(get_arguments()))
