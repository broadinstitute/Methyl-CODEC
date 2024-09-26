#!/usr/bin/env python
import sys
import argparse
import pandas as pd

def get_arguments():

    parser = argparse.ArgumentParser(prog="count mutations by monomer and pstrand type", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("variant_file", type=str, help="CODEC variant file")
    parser.add_argument("metrics_file", type=str, help="CODEC mutant metrics file")
    parser.add_argument("output_file", type=str, help="output monomer count file")
    args = parser.parse_args()
    return args

def absolute_path(path):
    """convert relative path to absolute path"""
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(os.getcwd(), path)



# Function to complement bases
def complement_base(base):
    complement = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
    return complement.get(base, base)

# Function to process multi-base mutations
def process_multibase(ref, alt, meth_char):
    mutations = []
    if len(ref) == len(alt):
        for r, a, m in zip(ref, alt, meth_char):
            mutations.append((r, a, m))
    else:
        mutations.append((ref, alt, meth_char))
    return mutations

# Process the mutations and count them
def process(options):

    # Filter for SNV type
    # Define all possible mutation classes
    all_mutations = ['A>T', 'A>C', 'A>5mC', 'A>G', 'T>A', 'T>C', 'T>5mC', 'T>G', 'C>A', 'C>T', 'C>G', 'G>A', 'G>T', 'G>C', 'G>5mC']
    # Initialize the mutation counts dictionary with all possible mutation classes
    mutation_counts = {mutation: 0 for mutation in all_mutations}
    with open(options.variant_file, 'r') as f:
        next(f)
        header = next(f).strip().split('\t')
        for line in f:
            row_line = line.strip().split('\t')
            row = dict(zip(header, row_line))
            if row['type'] != 'SNV':
                continue
            ref, alt, pstrand_orientation, meth_chars = row['ref'], row['alt'], int(row['pstrand_orientation']), row['meth_char']
            # Process multi-base mutations
            mutations = process_multibase(ref, alt, meth_chars)

            for mutation in mutations:
                ref_base, alt_base, meth_char = mutation[0], mutation[1], mutation[2]

                # Complement the bases if pstrand_orientation is reverse
                if pstrand_orientation == 2:
                    ref_base = complement_base(ref_base)
                    alt_base = complement_base(alt_base)

                mutation_str = f"{ref_base}>{alt_base}"
                if meth_char.isupper():
                    mutation_str = f"{ref_base}>5m{alt_base}"
                if mutation_str in mutation_counts:
                    mutation_counts[mutation_str] += 1
                else:
                    mutation_counts[mutation_str] = 1

    sorted_mutation_counts = dict(sorted(mutation_counts.items()))

    error_met_df = pd.read_csv(options.metrics_file, delimiter='\t')
    den = pd.DataFrame({
        "base": ["A", "C", "G", "T"],
        "count": [
            error_met_df["n_A_eval_f"].iloc[0] + error_met_df["n_T_eval_r"].iloc[0],
            error_met_df["n_C_eval_f"].iloc[0] + error_met_df["n_G_eval_r"].iloc[0],
            error_met_df["n_G_eval_f"].iloc[0] + error_met_df["n_C_eval_r"].iloc[0],
            error_met_df["n_T_eval_f"].iloc[0] + error_met_df["n_A_eval_r"].iloc[0]
        ]
    })

    with open(options.output_file, 'w') as fo:
    # Print the mutation counts
        print("muttype", "pstrand_o", "n_snv", "n_base_eval", "mut_rate", sep='\t', file=fo)
        for mutation, count in sorted_mutation_counts.items():
            muttype_first_base = mutation[0]
            n_base_eval = den[den["base"] == muttype_first_base]["count"].values[0]
            print(f"{mutation}\t1\t{count}\t{n_base_eval}\t{count/n_base_eval}", file=fo)

if __name__ == '__main__':
    sys.exit(process(get_arguments()))