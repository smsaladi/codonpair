"""
Python implementation to calculate Codon Pair Score

J. R. Coleman, D. Papamichail, S. Skiena, B. Futcher, E. Wimmer
and S. Mueller, "Virus attenuation by genome-scale changes in
codon-pair bias: a novel method for developing viral vaccines"
Science, Vol. 320, 1784-1787, 27 June 2008

Shyam Saladi (saladi@caltech.edu)
"""

import warnings
import re
from pkg_resources import resource_filename

import pandas as pd
import numpy as np
import Bio.Data.CodonTable

ref_pair_df = None
"""
"""

cur_ref_fn = None
"""
"""

def codon_pair_ref(fn=None, force_reload=False):
    """Calculate tAI for the provided nucleotide sequence

    Parameters
    ----------
    fn : Optional[str]
        Specifies the reference tRNA abundances used for calculation. `codonR`
        refers to the abundances given with the original package. If specified
        as a `dict`, abundances are taken directly where keys are codons and
        values are counts of cognate tRNAs and values are cached until
        `recalc_weights` = True.

    Returns
    -------
    pd.DataFrame
        XXXXX

    Raises
    ------
    None
    """
    global ref_pair_df, cur_ref_fn
    if ref_pair_df is None or force_reload or fn != cur_ref_fn:
        if fn is None:
            fn = resource_filename(__name__, 'data/ec_de3_ref.cps.tbd')
        else:
            cur_ref_fn = fn
        ref_pair_df = pd.read_table(fn, header=0, index_col=0, na_filter=False)

        if 'cps' not in ref_pair_df.columns:
            # reading a file created by old cps.pl
            ref_pair_df = pd.read_table(fn, header=None, usecols=[0,1,3],
                                        na_filter=False, index_col=1,
                                        names=['aa_pair', 'cp', 'cp_cnt'])

            ref_pair_df['codon1'] = ref_pair_df.index.str[:3]
            ref_pair_df['codon2'] = ref_pair_df.index.str[-3:]
            ref_pair_df['aa1'] = ref_pair_df['aa_pair'].str[0]
            ref_pair_df['aa2'] = ref_pair_df['aa_pair'].str[1]

        ref_pair_df = prep_ref_data(ref_pair_df)

    return ref_pair_df


def prep_ref_data(ref):
    """Calculate tAI for the provided nucleotide sequence

    Parameters
    ----------
    ref : pd.DataFrame
        XXXXX

    Returns
    -------
    pd.DataFrame
        tAI value for the nucseq specified with the options given

    Raises
    ------
    None
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
    tot_pairs = ref['cp_cnt'].sum()
    jp_factor = tot_pairs / (tot_pairs + 2048)
    jp_summant = 1 / (2 * (tot_pairs + 2048))
    ref['cps'] = np.log((ref['cp_cnt'] * jp_factor + jp_summant) *
                            ref['aa1_cnt'] * ref['aa2_cnt'] /
                            (ref['codon1_cnt'] * ref['codon2_cnt'] *
                             ref['aa_pair_cnt']))
    return ref


def get_pairs(seq):
    """Calculate codon pairs and relative frequency for each given a genome

    Parameters
    ----------
    seq : str
        XXXXX

    Returns
    -------
    pd.Series
        XXXXX

    Raises
    ------
    None
    """
    pairs = re.findall('......', seq)
    pairs.extend(re.findall('......', seq[3:-3])) # second frame
    return pd.Series(pairs, name='cp_cnt').value_counts(sort=False)


def summarize_type(df, value1, value2, total):
    """Calculate codon pairs and relative frequency for each given a genome

    Parameters
    ----------
    df : pd.DataFrame
        XXXXX

    value1 : str
        XXXXX

    value2 : str
        XXXXX

    total : str
        XXXXX

    Returns
    -------
    pd.DataFrame
        tAI value for the nucseq specified with the options given

    Raises
    ------
    None
    """
    # counts are divided by 2 becuase each internal codons are counted twice,
    # but this isn't *exactly* correct becuase it doesn't properly treat
    # edge effects, i.e. the last codon isn't double counted. Done this way
    # becuase the original implementation does this (c.a. line 224)

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
    """Calculate codon pairs and relative frequency for each given a genome

    Parameters
    ----------
    seqlist : iterable of `Bio.SeqRecord`
        XXXXXX

    codon_table : Optional[str]
        XXXXXX

    Returns
    -------
    pd.DataFrame
        XXXXX

    Raises
    ------
    None
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
    global ref_pair_df
    ref_pair_df = pairs
    return pairs


def calc_cpb(seq, ref_fn=None):
    """Calculate codon pairs and relative frequency for each given a genome

    Parameters
    ----------
    seq : str

    ref_fn : Optional[str]
        asdf

    Returns
    -------
    float
        tAI value for the nucseq specified with the options given

    Raises
    ------
    None

    Calculate codon pair bais for a given gene against a reference set
    "CPB is the arithmetic mean of the individual codon pair scores of all
    pairs making up the ORF." (Caption of Fig. S1 in Coleman, 2008). Equation
    shown in part S1B seems to be written incorrectly.

    #   $B = $powers_of_4[6]; 4**6 = 4096
    #   $Bl = $B/2; 4096/2 = 2048
    #   $jp_factor = $total_codon_pairs / ($total_codon_pairs + 2048);
    #   $jp_summant = 1/(2*($total_codon_pairs + $Bl));
    #   $my $discounted_num = &discount($4 + 0);
    #   $codon_pair{$2}->{num} = $4*$jp_factor + $jp_summant;
    """

    ref = codon_pair_ref(fn=ref_fn)

    pairs = pd.DataFrame(get_pairs(seq))

    pairs['cp_cnt_cps'] = pairs['cp_cnt'] * ref['cps']

    # Codon pairs present in the sequence of interest but not the reference
    # could indicate an issue with the sequence provided
    if pairs['cp_cnt_cps'].isnull().values.any():
        warnings.warn('Codon pair found in sequence not found in reference.')
        #pairs = pairs[pd.notnull(pairs['pairs_ref'])]

    cps_sum = pairs['cp_cnt_cps'].sum()
    pair_count = pairs['cp_cnt'].sum()

    return (cps_sum, pair_count, cps_sum/pair_count)


def write_reference(seq_records, out_fn):
    """Calculate codon pairs and relative frequency for each given a genome

    Parameters
    ----------
    nucseq : str

    seq_records : list of Bio.SeqRecord
        Specifies the reference tRNA abundances used for calculation. `codonR`
        refers to the abundances given with the original package. If specified
        as a `dict`, abundances are taken directly where keys are codons and
        values are counts of cognate tRNAs and values are cached until
        `recalc_weights` = True.

    out_fn : str
        Weights for each codon are calculated only once using `calc_weights`
        and then cached/stored in `trna_data`. If `True`, weights will be
        recalculated.

    Returns
    -------
    None

    Raises
    ------
    None
    """
    calc_reference(seq_records).to_csv(out_fn, sep='\t')
    return
