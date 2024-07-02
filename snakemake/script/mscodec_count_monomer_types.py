#!/usr/bin/env python
import sys
import argparse
import pandas as pd

def get_arguments():

    parser = argparse.ArgumentParser(prog="count mutations by monomer and pstrand type", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_file", type=str, help="CODEC variant file")
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
def process_multibase(ref, alt):
    mutations = []
    if len(ref) == len(alt):
        for r, a in zip(ref, alt):
            mutations.append((r, a))
    else:
        mutations.append((ref, alt))
    return mutations

# Process the mutations and count them
def process(options):
    # Read the file and skip lines starting with '#'
    df = pd.read_csv(options.input_file, delimiter='\t', comment='#')

    # Filter for SNV type
    df_snv = df[df['type'] == 'SNV']
    # Define all possible mutation classes
    all_mutations = ['A>T', 'A>C', 'A>G', 'T>A', 'T>C', 'T>G', 'C>A', 'C>T', 'C>G', 'G>A', 'G>T', 'G>C']
    # Initialize the mutation counts dictionary with all possible mutation classes
    mutation_counts = {mutation: 0 for mutation in all_mutations}
    for index, row in df_snv.iterrows():
        ref, alt, pstrand_orientation = row['ref'], row['alt'], row['pstrand_orientation']

        # Process multi-base mutations
        mutations = process_multibase(ref, alt)

        for mutation in mutations:
            ref_base, alt_base = mutation[0], mutation[1]

            # Complement the bases if pstrand_orientation is reverse
            if pstrand_orientation == 2:
                ref_base = complement_base(ref_base)
                alt_base = complement_base(alt_base)

            mutation_str = f"{ref_base}>{alt_base}"
            if mutation_str in mutation_counts:
                mutation_counts[mutation_str] += 1
            else:
                mutation_counts[mutation_str] = 1

    sorted_mutation_counts = dict(sorted(mutation_counts.items()))
    with open(options.output_file, 'w') as fo:
    # Print the mutation counts
        print("muttype", "pstrand_o", "n_snv", sep='\t', file=fo)
        for mutation, count in sorted_mutation_counts.items():
            print(f"{mutation}\t1\t{count}", file=fo)

if __name__ == '__main__':
    sys.exit(process(get_arguments()))