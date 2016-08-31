# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import Bio.SeqIO

import cps

def cpb_vs_ref(fna_file, ref_csv, **kwargs):
    # Read reference features
    calc = pd.read_table(ref_csv, header=None, names=['old_cps_sum', 'old_cpb'])
    with open(fna_file, 'r+') as fh:
        calc['seq'] = pd.Series([str(rec.seq)
                                    for rec in Bio.SeqIO.parse(fh, 'fasta')])

    calc[['new_cps_sum', 'new_pair_count', 'new_cpb']] = \
        calc['seq'].apply(cps.calc_cpb, **kwargs).apply(pd.Series)

    assert np.allclose(calc['old_cps_sum'].values, calc['new_cps_sum'].values,
                       atol=1e-3, rtol=1)
    assert np.allclose(calc['old_cpb'].values, calc['new_cpb'].values,
                       atol=1e-3, rtol=1)

    return


def cps_vs_ref(fna_file, ref_csv):
    old_cps = pd.read_table(ref_csv, header=None, usecols=[0,1,3],
                            names=['aa_pair', 'cp', 'cp_cnt'], index_col=1)

    new_cps = cps.calc_reference(Bio.SeqIO.parse(fna_file, "fasta"))

    old_cps, new_cps = old_cps.align(new_cps, axis=0)

    assert np.allclose(old_cps['cp_cnt'].values,
                       new_cps['cp_cnt'].values, atol=1e-4)

    return


def test_sp_genome_sp_cps():
    cps_vs_ref("../data/sp_genome.fna", "../data/sp_ref.cps.tbd")
    return


def test_ecolik12_ecde3():
    basename = "data/ecolik12.ffn"
    cpb_vs_ref(basename, basename + ".ec_cpb")
    return


def test_ecolik12_sp():
    basename = "data/ecolik12.ffn"
    cpb_vs_ref(basename, basename + ".sp_cpb", ref_fn="../data/sp_ref.cps.tbd")
    return


def test_daley_gfp_de3():
    basename = "data/Daley_gfp.fna"
    cpb_vs_ref(basename, basename + ".ec_cpb")
    return


def test_daley_gfp_sp():
    basename = "data/Daley_gfp.fna"
    cpb_vs_ref(basename, basename + ".sp_cpb", ref_fn="../data/sp_ref.cps.tbd")
    return
