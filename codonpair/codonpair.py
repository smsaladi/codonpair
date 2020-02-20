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
from collections import OrderedDict

from pkg_resources import resource_filename
import argparse

import pandas as pd
import numpy as np

import Bio.SeqIO
import Bio.Data.CodonTable

        
class CodonPair():
    def __init__(self, df_pair_count=None):
        if df_pair_count is None:
            obj = self.__class__.from_named_reference('E. coli')
            self.df_ref = obj.df_ref
        else:
            cols = df_pair_count.columns
            req_cols = ['codon1', 'codon2', 'aa1', 'aa2', 'cp_cnt']
            missing_cols = set(req_cols) - set(cols)
            if len(missing_cols) != 0:
                raise ValueError("df_pair_count is missing cols: {}. ".format(", ".join(missing_cols)))

            self.df_ref = _prep_ref_data(df_pair_count)
        return

    @classmethod
    def from_named_reference(cls, name):
        """
        """
        if name == 'E. coli':
            ref_fn = resource_filename(__name__, 'data/ec_de3_ref.cps.tbd')
        elif name == 'cps_perl':
            ref_fn = resource_filename(__name__, 'data/ec_de3_ref.cps.tbd')
        elif name == 'S. pneumoniae':
            ref_fn = resource_filename(__name__, 'data/sp_ref.cps.tbd')
        else:
            raise ValueError("Unrecognized name: " + name)
        return cls.from_reference_file(ref_fn)

    @classmethod
    def from_reference_file(cls, fn):
        """Uses the codon pairs specified by fn to be used as reference
        """
        df = pd.read_table(fn, header=0, index_col=0, na_filter=False)
        if 'cps' not in df.columns:
            # reading a file created by old cps.pl
            df = pd.read_table(fn, header=None, usecols=[0,1,3],
                               na_filter=False, index_col=1,
                               names=['aa_pair', 'cps', 'cp_cnt'])

            df['codon1'] = df.index.str[:3]
            df['codon2'] = df.index.str[-3:]
            df['aa1'] = df['aa_pair'].str[0]
            df['aa2'] = df['aa_pair'].str[1]

        return cls(df)

    @classmethod
    def from_sequences(cls, seqlist, codon_table='Standard'):
        """Calculate codon pairs and relative frequency for a sequence set

        `seqlist` expected to be a iterable of `Bio.SeqRecord.SeqRecord`s
        `codon_table` refers to keys given of
        `Bio.Data.CodonTable.unambiguous_dna_by_name`
        """
        # Set up data frame for counting
        df_codons = pd.DataFrame.from_dict(
            Bio.Data.CodonTable.unambiguous_dna_by_name[codon_table].forward_table,
            orient='index')

        df_pairs = pd.DataFrame(index=pd.MultiIndex.from_product(
            [df_codons.index, df_codons.index], names=['codon1', 'codon2']
        ))
        df_pairs['cp_cnt'] = 0
        df_pairs.reset_index(inplace=True)
        df_pairs.index = df_pairs['codon1'] + df_pairs['codon2']
        df_pairs.index.name = 'cp'

        # Match with amino acids
        df_codons.columns = ['aa1']
        df_pairs = pd.merge(df_pairs, df_codons, how="outer", left_on='codon1',
                        right_index=True, sort=False, copy=False)

        df_codons.columns = ['aa2']
        df_pairs = pd.merge(df_pairs, df_codons, how="outer", left_on='codon2',
                        right_index=True, sort=False, copy=False)

        for s in seqlist:
            df_pairs['cp_cnt'] = df_pairs['cp_cnt'].add(_get_pairs(s), fill_value=0)

        return cls(df_pairs)

    def write_reference(self, out_fn):
        self.df_ref.to_csv(out_fn, sep='\t')
        return

    def cpb(self, seq):
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
        df_pair = pd.DataFrame(_get_pairs(seq))

        df_pair['cp_cnt_cps'] = df_pair['cp_cnt'] * self.df_ref['cps']

        # Codon pairs present in the sequence of interest but not the reference
        # could indicate an issue with the sequence provided
        if df_pair['cp_cnt_cps'].isnull().values.any():
            warnings.warn('Codon pair found in sequence not found in reference.')
            #pairs = pairs[pd.notnull(pairs['pairs_ref'])]

        cps_sum = df_pair['cp_cnt_cps'].sum()
        n_pair = np.sum(df_pair['cp_cnt'] * df_pair['cp_cnt_cps'].notnull())

        return OrderedDict(total_cps=cps_sum, n_pair=n_pair, cpb=cps_sum/n_pair)


_re_pairs = re.compile('......')
def _get_pairs(seq):
    """Calculates codon pair counts given a sequence as string
    """
    pairs = _re_pairs.findall(seq)
    pairs.extend(_re_pairs.findall(seq[3:])) # second frame
    return pd.Series(pairs, name='cp_cnt').value_counts(sort=False)


def _summarize_type(df, value1, value2, total):
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

def _prep_ref_data(df_ref):
    """Prepares the codon pair counts for CPS calculation
    """
    # each codon's individial count
    df_ref = _summarize_type(df_ref, 'codon1', 'codon2', 'cp_cnt')
    # each amino acid's individual count
    df_ref = _summarize_type(df_ref, 'aa1', 'aa2', 'cp_cnt')

    if 'aa_pair' not in df_ref.columns:
        df_ref['aa_pair'] = df_ref['aa1'] + df_ref['aa2']

    # aa pair count
    aa_pair = df_ref[['aa_pair', 'cp_cnt']].groupby('aa_pair').sum()
    aa_pair.columns = ['aa_pair_cnt']
    df_ref = pd.merge(df_ref, aa_pair, how='left', left_on='aa_pair',
                      right_index=True, sort=False, copy=False)

    # Codon Pair Score adjusted by Jeffreys-Perks Law
    # https://en.wikipedia.org/wiki/Additive_smoothing#Pseudocount
    # alpha = 0.5 gives Jefferys-Perks Law
    n_pair = df_ref['cp_cnt'].sum()
    jp_factor = n_pair / (n_pair + 2048)
    jp_summant = 1 / (2 * (n_pair + 2048))
    df_ref['cps'] = np.log(
        (df_ref['cp_cnt'] * jp_factor + jp_summant) * df_ref['aa1_cnt'] * df_ref['aa2_cnt'] /
        (df_ref['codon1_cnt'] * df_ref['codon2_cnt'] * df_ref['aa_pair_cnt'])
    )
    return df_ref

def calc_cpb(seq):
    """Convenience method to calculate codon pair bias using the defaults for the provided sequence
    """
    return CodonPair().cpb(seq)

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

    if args.ref is None:
        cp = CodonPair()
    else:
        cp = CodonPair.from_reference_file(args.ref)

    if args.fna_filename == '-':
        seq_records = [str(rec.seq) for rec in Bio.SeqIO.parse(sys.stdin, "fasta")]
    else:
       with open(args.fna_filename, 'r') as fh:
            seq_records = [str(rec.seq) for rec in Bio.SeqIO.parse(fh, "fasta")]

    if args.calc_reference:
        cp = CodonPair.from_sequences(seq_records)
        cp.write_reference(args.output)
    else:
        for s in seq_records:
            print(cp.cpb(s))

    return

if __name__ == '__main__':
    main()
