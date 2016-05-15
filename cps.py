'''
Python implementation to calculate Codon Pair Score

J. R. Coleman, D. Papamichail, S. Skiena, B. Futcher, E. Wimmer
and S. Mueller, "Virus attenuation by genome-scale changes in
codon-pair bias: a novel method for developing viral vaccines"
Science, Vol. 320, 1784-1787, 27 June 2008

Shyam Saladi (saladi@caltech.edu)
April 2016
'''

# Only used in main()
import sys
import os.path
import warnings
import argparse

# used in codon pair calculation
import re
import collections

# used in codon pair calculation
import pandas as pd
import numpy as np
import scipy.stats.mstats

# Only used in main()
import Bio.SeqIO
import Bio.Data.CodonTable

global ref_trna
trna_table = pd.read_csv('~/ml-expression/scales/codon_abundances.csv', \
    header=0, index_col=0, comment='#')

def get_pairs(seq):
    pairs = pd.Series(re.findall('......', seq)).value_counts(sort=False)
    pairs.name = 'pair_count'
    return pairs


def calc_reference(seqlist):
    '''
    Calculate codon pairs and relative frequency for each given a genome
    '''

    all_pairs = pd.Series(name='pair_count')

    for record in seqlist:
        all_pairs = all_pairs.add(get_pairs(str(record.seq)), fill_value=0)

    return process_pairs(all_pairs)


def process_pairs(pairs_df):
    '''
    Process a set of codon pairs for a sequence

    relative frequency = f(codon1) * f(codon2) / f(codon1-codon2)
    '''
    fun_split = lambda x: pd.Series([x[:3], x[3:]], index=['codon1', 'codon2'])
    codon_names = pd.DataFrame(pairs_df.index.to_series().apply(fun_split))

    pairs_df = pd.concat([codon_names, pairs_df], axis=1)

    # Individual counts
    codon1_count = pairs_df[['codon1', 'pair_count']].groupby('codon1').sum()
    codon2_count = pairs_df[['codon2', 'pair_count']].groupby('codon2').sum()
    codon_count = codon1_count.add(codon2_count, fill_value=0)

    codon_count.columns = ['codon1_count']
    pairs_df = pairs_df.merge(codon_count, left_on='codon1',
                              right_index=True, copy=False)
    codon_count.columns = ['codon2_count']
    pairs_df = pairs_df.merge(codon_count, left_on='codon2',
                              right_index=True, copy=False)

    pairs_df['relative_freq'] = pairs_df['pair_count'] / \
        (pairs_df['codon1_count'] * pairs_df['codon2_count'])

    return pairs_df


def calc_cpb(seq, ref):
    '''
    Calculate codon pair bais for a given gene against a reference set
    '''

    pairs = process_pairs(get_pairs(seq))
    pairs['pairs_ref'] = pairs['relative_freq'].divide(ref['relative_freq'])

    # Remove those codon pairs present in one set but not the other
    pairs = pairs[pd.notnull(pairs['pairs_ref'])]

    pairs['pairs_ref'] = pairs['pairs_ref'].apply(np.log)

    return pairs['pairs_ref'].sum() / len(pairs.index)


def main():
    parser = argparse.ArgumentParser(
        description='Calculate codon pair score/bias.')

    parser.add_argument('fna_filename',
        metavar='fna_file',
        type=str,
        default=sys.stdin,
        help='FASTA-formatted nucleotide file of the Coding Sequences'
             '(stdin by default)')

    parser.add_argument('--ref',
        metavar='reference_file',
        type=str,
        help='Calculate reference pairs for provided file')

    parser.add_argument('--output',
        metavar='output_file',
        type=str,
        default=sys.stdout,
        help='Specifies name of file to write to (stdout by default)')

    parser.add_argument('--calc_reference',
        action='store_true',
        help='Calculate reference pairs for provided file')

    args = parser.parse_args()

    if isinstance(args.fna_filename, str):
        with open(args.fna_filename, 'r+') as fh:
            seq_records = list(Bio.SeqIO.parse(fh, "fasta"))
    else:
        print("Reading sequences from stdin", file=sys.stderr)
        seq_records = Bio.SeqIO.parse(args.fna_file, "fasta")

    if args.calc_reference:
        # extra blank row because of multiindex
        # https://github.com/pydata/pandas/issues/6618
        calc_reference(seq_records).to_csv(args.output)
    else:
        pair_ref = pd.read_csv(args.ref, index_col=0)
        for record in seq_records:
            print(calc_cpb(str(record.seq), pair_ref))

if __name__ == '__main__':
    main()
