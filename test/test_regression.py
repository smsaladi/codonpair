# -*- coding: utf-8 -*-

import os.path
import pkg_resources

import pandas as pd
import numpy as np
import Bio.SeqIO

import codonpair

TEST_DIR = os.path.dirname(os.path.abspath(__file__))

def cpb_vs_ref(fna_file, ref_csv, ref_fn=None):
    # Read reference features
    calc = pd.read_table(ref_csv, header=None, names=['old_cps_sum', 'old_cpb'])
    with open(fna_file, 'r+') as fh:
        calc['seq'] = \
            pd.Series([str(rec.seq) for rec in Bio.SeqIO.parse(fh, 'fasta')])

    if ref_fn is None:
        cp = codonpair.CodonPair()
    else:
        cp = codonpair.CodonPair.from_reference_file(ref_fn)
    
    calc[['new_cps_sum', 'new_cpb']] = \
        calc['seq'].apply(lambda x: pd.Series(cp.cpb(x))[['total_cps', 'cpb']])

    diff_vals = (
            ~np.isclose(calc['old_cps_sum'], calc['new_cps_sum']) |
            ~np.isclose(calc['old_cpb'], calc['new_cpb'])
    )
    if diff_vals.sum() > 0:
        print("Differing values:")
        print(calc[diff_vals])

    np.testing.assert_allclose(calc['old_cps_sum'], calc['new_cps_sum'])
    np.testing.assert_allclose(calc['old_cpb'], calc['new_cpb'])

    return


def cps_vs_ref(fna_file, ref_csv):
    old_cps = pd.read_table(ref_csv, header=None, usecols=[0,1,3],
                            names=['aa_pair', 'cp', 'cp_cnt'], index_col=1)

    cp = codonpair.CodonPair.from_sequences([
        str(rec.seq) for rec in Bio.SeqIO.parse(fna_file, "fasta")
    ])
    new_cps = cp.df_ref

    old_cps, new_cps = old_cps.align(new_cps, axis=0)

    np.testing.assert_allclose(old_cps['cp_cnt'], new_cps['cp_cnt'])
    return


def pkgfile(fn):
    return pkg_resources.resource_filename("codonpair", fn)

def test_sp_genome_sp_cps():
    cps_vs_ref(pkgfile('data/sp_genome.fna'),
               pkgfile('data/sp_ref.cps.tbd'))
    return

def test_ecolik12_ecde3():
    basename = TEST_DIR + "/data/ecolik12.ffn"
    cpb_vs_ref(basename, basename + ".ec_cpb")
    return

def test_ecolik12_sp():
    basename = TEST_DIR + "/data/ecolik12.ffn"
    cpb_vs_ref(basename, basename + ".sp_cpb",
               ref_fn=pkgfile("data/sp_ref.cps.tbd"))
    return

def test_daley_gfp_de3():
    basename = TEST_DIR + "/data/Daley_gfp.fna"
    cpb_vs_ref(basename, basename + ".ec_cpb")
    return

def test_daley_gfp_sp():
    basename = TEST_DIR + "/data/Daley_gfp.fna"
    cpb_vs_ref(basename, basename + ".sp_cpb",
               ref_fn=pkgfile("data/sp_ref.cps.tbd"))
    return
