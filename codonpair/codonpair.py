"""
Python implementation to calculate Codon Pair Score

J. R. Coleman, D. Papamichail, S. Skiena, B. Futcher, E. Wimmer
and S. Mueller, "Virus attenuation by genome-scale changes in
codon-pair bias: a novel method for developing viral vaccines"
Science, Vol. 320, 1784-1787, 27 June 2008

Shyam Saladi (saladi@caltech.edu)
"""

from __future__ import print_function, division

import sys
import warnings
import re
from pkg_resources import resource_filename
import argparse

import pandas as pd
import numpy as np

import Bio.SeqIO
import Bio.Data.CodonTable

def codon_pair_ref(fn):
    """Uses the codon pairs specified by fn to be used as reference
    """
    ref_pair_df = pd.read_table(fn, header=0, index_col=0, na_filter=False)
    if 'cps' not in ref_pair_df.columns:
        # reading a file created by old cps.pl
        ref_pair_df = pd.read_table(fn, header=None, usecols=[0,1,3],
                                    na_filter=False, index_col=1,
                                    names=['aa_pair', 'cps', 'cp_cnt'])

        ref_pair_df['codon1'] = ref_pair_df.index.str[:3]
        ref_pair_df['codon2'] = ref_pair_df.index.str[-3:]
        ref_pair_df['aa1'] = ref_pair_df['aa_pair'].str[0]
        ref_pair_df['aa2'] = ref_pair_df['aa_pair'].str[1]
    ref_pair_df = prep_ref_data(ref_pair_df)
    return ref_pair_df


def prep_ref_data(ref):
    """Prepares the codon pair counts for CPS calculation

    Helper for `codon_pair_ref` (above)
    """
    # each codon's individial count
    ref = summarize_type(ref, 'codon1', 'codon2', 'cp_cnt')
    # each amino acid's individual count
    ref = summarize_type(ref, 'aa1', 'aa2', 'cp_cnt')

    # aa pair count
    aa_pair = ref[['aa_pair', 'cp_cnt']].groupby('aa_pair').sum()
    aa_pair.columns = ['aa_pair_cnt']
    ref = pd.merge(ref, aa_pair, how='left', left_on='aa_pair',
                   right_index=True, sort=False, copy=False)

    # Codon Pair Score adjusted by Jeffreys-Perks Law
    # https://en.wikipedia.org/wiki/Additive_smoothing#Pseudocount
    # alpha = 0.5 gives Jefferys-Perks Law
    tot_pairs = ref['cp_cnt'].sum()
    jp_factor = tot_pairs / (tot_pairs + 2048)
    jp_summant = 1 / (2 * (tot_pairs + 2048))
    ref['cps'] = np.log((ref['cp_cnt'] * jp_factor + jp_summant) *
                            ref['aa1_cnt'] * ref['aa2_cnt'] /
                            (ref['codon1_cnt'] * ref['codon2_cnt'] *
                             ref['aa_pair_cnt']))
    return ref


def get_pairs(seq):
    """Calculates codon pair counts given a sequence as string
    """
    pairs = re.findall('......', seq)
    pairs.extend(re.findall('......', seq[3:])) # second frame
    return pd.Series(pairs, name='cp_cnt').value_counts(sort=False)


def summarize_type(df, value1, value2, total):
    """Helper function for calculating a codon pair reference
    """
    # counts are divided by 2 becuase each internal codons are counted twice,
    # but this isn't *exactly* correct becuase it doesn't properly treat
    # edge effects, i.e. the last codon isn't double counted. Done this way
    # becuase the original implementation (`cps_perl`) does this (c.a. line 224)

    value1_count = df[[value1, total]].groupby(value1).sum() # /2
    value2_count = df[[value2, total]].groupby(value2).sum() # /2
    value_count = value1_count.add(value2_count, fill_value=0)

    value_count.columns = [value1+'_cnt']

    df = df.merge(value_count, how='left', left_on=value1, right_index=True,
                  copy=False)

    value_count.columns = [value2+'_cnt']
    return df.merge(value_count, how='left', left_on=value2, right_index=True,
                    copy=False)


def calc_reference(seqlist, codon_table='Standard'):
    """Calculate codon pairs and relative frequency for a sequence set

    `seqlist` expected to be a iterable of `Bio.SeqRecord.SeqRecord`s
    `codon_table` refers to keys given of
    `Bio.Data.CodonTable.unambiguous_dna_by_name`
    """
    # Set up data frame for counting
    codons = pd.DataFrame.from_dict(
        Bio.Data.CodonTable.unambiguous_dna_by_name[codon_table].forward_table,
        orient='index')

    pairs = pd.DataFrame({'cp_cnt': [0]*len(codons.index)**2},
            index=pd.MultiIndex.from_product([codons.index, codons.index],
                                             names=['codon1', 'codon2']))
    pairs.reset_index(inplace=True)
    pairs.index = pairs['codon1'] + pairs['codon2']
    pairs.index.name = 'cp'

    # Match with amino acids
    codons.columns = ['aa1']
    pairs = pd.merge(pairs, codons, how="outer", left_on='codon1',
                     right_index=True, sort=False, copy=False)

    codons.columns = ['aa2']
    pairs = pd.merge(pairs, codons, how="outer", left_on='codon2',
                     right_index=True, sort=False, copy=False)

    for record in seqlist:
        pairs['cp_cnt'] = pairs['cp_cnt'].add(get_pairs(str(record.seq)),
                                              fill_value=0)
    return pairs


def calc_cpb(seq, ref_fn=resource_filename(__name__, 'data/ec_de3_ref.cps.tbd')):
    """Calculate codon pair bais for a given gene against a reference set

    Returns: cps_sum, pair_count, cps_sum/pair_count

    "CPB is the arithmetic mean of the individual codon pair scores of all
    pairs making up the ORF." (Caption of Fig. S1 in Coleman, 2008). Equation
    shown in part S1B seems to be written incorrectly.

    ## From the previous implementation (see `cps_perl`)
        ```perl
        $B = $powers_of_4[6]; 4**6 = 4096
        $Bl = $B/2; 4096/2 = 2048
        $jp_factor = $total_codon_pairs / ($total_codon_pairs + 2048);
        $jp_summant = 1/(2*($total_codon_pairs + $Bl));
        $my $discounted_num = &discount($4 + 0);
        $codon_pair{$2}->{num} = $4*$jp_factor + $jp_summant;
        ```
    """

    ref = codon_pair_ref(ref_fn)

    pairs = pd.DataFrame(get_pairs(seq))

    pairs['cp_cnt_cps'] = pairs['cp_cnt'] * ref['cps']

    # Codon pairs present in the sequence of interest but not the reference
    # could indicate an issue with the sequence provided
    if pairs['cp_cnt_cps'].isnull().values.any():
        warnings.warn('Codon pair found in sequence not found in reference.')
        #pairs = pairs[pd.notnull(pairs['pairs_ref'])]

    cps_sum = pairs['cp_cnt_cps'].sum()
    pair_count = np.sum(pairs['cp_cnt'] * pairs['cp_cnt_cps'].notnull())

    return (cps_sum, pair_count, cps_sum/pair_count)


def write_reference(seq_records, out_fn):
    calc_reference(seq_records).to_csv(out_fn, sep='\t')
    return


def main():
    parser = argparse.ArgumentParser(
        description='Calculate codon pair score/bias.')

    parser.add_argument('fna_filename',
        metavar='fna_file',
        type=str,
        default='-',
        help='FASTA-formatted nucleotide file of the Coding Sequences'
             '(stdin by default)')

    parser.add_argument('--ref',
        metavar='reference_file',
        type=str,
        default=None,
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

    if args.ref is not None:
        codon_pair_ref(args.ref)

    if args.fna_filename == '-':
        print("Reading sequences from stdin", file=sys.stderr)
        seq_records = Bio.SeqIO.parse(args.fna_file, "fasta")
    else:
       with open(args.fna_filename, 'r+') as fh:
            seq_records = list(Bio.SeqIO.parse(fh, "fasta"))

    if args.calc_reference:
        write_reference(seq_records, args.output)
    else:
        for record in seq_records:
            print(calc_cpb(str(record.seq)))


if __name__ == '__main__':
    main()

