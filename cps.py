'''
Python adaptation of Dimitris Papamichail's 2005-11-20 utility to calculate
Codon Pair Score (cps.pl) to calculate reference codon pairs from an input
dataset (e.g. a genome).

Papamichail's Perl scripts were written as part of the following publication:

J. R. Coleman, D. Papamichail, S. Skiena, B. Futcher, E. Wimmer
and S. Mueller, "Virus attenuation by genome-scale changes in
codon-pair bias: a novel method for developing viral vaccines"
Science, Vol. 320, 1784-1787, 27 June 2008

Shyam Saladi (saladi@caltech.edu)
April 2016
'''

# Only used in main() to read fasta files
# not for tAI calculation
import sys
import os.path
import warnings

# used in tAI calculation
import re
import collections

# used in tAI calculation
import pandas as pd
import numpy as np
import scipy.stats.mstats

# Only used in main() to read fasta files
# not for tAI calculation
import Bio.SeqIO

global ref_trna
trna_table = pd.read_csv('~/ml-expression/scales/codon_abundances.csv', \
    header=0, index_col=0, comment='#')

def test_calc_tAI():
    # do tests
    return

def calc_tAI(nucseq, ref_trna='codonR', bacteria=True, recalc_weights=False):
    '''
    This is a rewrite of Mario dos Reis's (m.reis@mail.cryst.bbk.ac.uk)
    popular codonR program in Python.

    It should be more extensible than the original program though perhaps not
    faster.

    If there are codons that have non-ATGCU letters, those codons will
    silently ignored from analysis

    There are small differences in the calculation compared to the original
    codonR program (O(1e-3)). This should be revisited.
    '''
    # Set up reference tRNA abundance source, i.e. which column
    # or a user supplied dict
    if ref_trna is not dict and ref_trna in trna_table.columns:
        ref_trna_count = trna_table[ref_trna]
    else:
        raise ValueError

    # Now we calculate the relative adaptiveness values for each codon for
    # the reference genome (if not already calculated)
    # Specify a new column to store the ln(weights) in:
    lnweights_colname = ref_trna + '_ln_weights' + bacteria*'_bact'
    if lnweights_colname not in trna_table.columns or recalc_weights:
        # Only do this once, so store after calculating
        trna_table[lnweights_colname] = \
            calc_weights(trna_count=ref_trna_count, bacteria=bacteria)

    ln_weights = trna_table[lnweights_colname]
    ln_weights.name = 'ln_weights'

    # We will ignore Methionine codons in our analysis (there is no
    # automatic way to differentiate between 'START' Met-tRNA genes and
    # normal Met-tRNAs in any genome):
    if 'ATG' in ln_weights.index:
        ln_weights.drop('ATG', inplace=True)

    # Check nucseq (just in case)
    nucseq = nucseq.upper().replace('U', 'T')
    # Now, count codons:
    # codons = collections.Counter(re.findall('...', nucseq))
    # codon_count = pd.Series(codons, name='codons')
    codon_id, codon_count = np.unique(re.findall('...', nucseq),
                                      return_counts=True)
    codon_count = pd.Series(codon_count, index=codon_id, name='codons')

    # remove non-standard codons from analysis by joining with weights series
    df = pd.concat([ln_weights, codon_count], axis = 1).loc[ln_weights.index]

    # and now we an finally calculate tAI:
    # tai is the weighted average of the weights (weighted by the codon counts)
    codon_prop = df['codons'] / df['codons'].sum()
    return np.exp(np.sum(df['ln_weights'] * codon_prop))

    # This may be faster
    # intermed = np.average(df['ln_weights'], weights=df['codons'], returned=True)
    # return np.exp(intermed[0]/intermed[1])

    # For more details, read the references!
    # [1] dos Reis et al. (2003) Nuc. Acids Res. 31:6976
    # [2] dos Reis et al. (2004) Nuc. Acids Res. 32:5036


def test_calc_weights():
    # do tests
    return

def calc_weights(trna_count, bacteria, optimized_weights=True, \
    keep_codonR_err=True):
    """
    Function to calculate relative adaptiveness values
    tRNA -- tRNA gene copy number pandas series
    """

    if optimized_weights:
        p = {'T': 0.59, 'C': 0.72, 'A': 0.0001, 'G': 0.32}

        # s is what the original script works with
        # p = 1 - s
        # s <- c(0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68)

        isoleucine_p = 1 - 0.89
    else:
        p = {'T': 0.5, 'C': 0.5, 'A': 0.25, 'G': 0.5}
        # s <- c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5)
        isoleucine_p = 1 - 0.5

    # obtain absolute adaptiveness values (Ws)

    # add zero's for stop codons since (e.g. 'TGN' may need a 'TGA' abundance)
    trna_count = pd.concat([trna_count,
        pd.Series(0, index=['TGA', 'TAA', 'TAG'])], axis=0)

    # trna file with codonR has a transposition resulting in erronous weights
    # for TGN codons
    if keep_codonR_err:
        trna_count['TGA'] = 1
        trna_count['TGC'] = 1
        trna_count['TGT'] = 0

    # for new weights to be calculated
    weights = pd.Series(0.0, index=trna_count.index)
    # number of tRNAs corresponds to anticodons (not codons)
    for codon in list(weights.index):
        base = codon[:2]
        if codon[2] == 'T':                  # INN -> NNT, NNC, NNA
            weights[codon] = trna_count[codon] + p['T']*trna_count[base+'C']
        elif codon[2] == 'C':                # GNN -> NNT, NNC
            weights[codon] = trna_count[codon] + p['C']*trna_count[base+'T']
        elif codon[2] == 'A':                # TNN -> NNA, NNG
            weights[codon] = trna_count[codon] + p['A']*trna_count[base+'T']
        elif codon[2] == 'G':                # CNN -> NNG
            weights[codon] = trna_count[codon] + p['G']*trna_count[base+'A']
        else:
            raise ValueError('Non-standard codon or notation')

    # Weight calculation in original script
    #   p[1]*tRNA[i] + p[5]*tRNA[i+1],     # INN -> NNT, NNC, NNA
    #   p[2]*tRNA[i+1] + p[6]*tRNA[i],       # GNN -> NNT, NNC
    #   p[3]*tRNA[i+2] + p[7]*tRNA[i],       # TNN -> NNA, NNG
    #   p[4]*tRNA[i+3] + p[8]*tRNA[i+2])     # CNN -> NNG

    # i correspondance based on ordering of input file (i.e. ecolik12.trna)
    # i -> T
    # i+1 -> C
    # i+2 -> A
    # i+3 -> G

    # if bacteria, modify isoleucine ATA codon
    if bacteria:
        weights['ATA'] = isoleucine_p

    # get rid of stop codons and methionine
    weights.drop(['ATG', 'TGA', 'TAA', 'TAG'], inplace=True)

    # get ws
    weights = weights/weights.max()

    # substitute 0-ws by gm
    nonzero_weights = weights[weights != 0]
    geometric_mean = scipy.stats.mstats.gmean(nonzero_weights)
    weights[weights == 0] = geometric_mean

    return np.log(weights)

def main():
    if len(sys.argv) == 1:
        # read from sys.stdin
        fh = sys.stdin
        for record in Bio.SeqIO.parse(fh, 'fasta'):
            print(record.id, calc_tAI(str(record.seq)))
    elif os.path.exists(sys.argv[1]):
        # try reading argument as file
        with open(sys.argv[1]) as fh:
            for record in Bio.SeqIO.parse(fh, 'fasta'):
                print(record.id, calc_tAI(str(record.seq)))
    elif sys.argv[1]:
#        warnings.warn('Trying argument as sequence itself')
        print(sys.argv[1], calc_tAI(sys.argv[1]))


if __name__ == '__main__':
    main()

my @powers_of_4;
$powers_of_4[0] = 1;
for (my $i = 1; $i < 16; $i++)
{
  push (@powers_of_4, $powers_of_4[$i-1]<<2);
}

my $CODON_FILE_LOCATION = '/temp/saladi/cps/codons';
my $DEFAULT_CODON_PAIR_INFO_FILE = '/temp/saladi/cps/Escherichia_coli__BL21_Gold_DE3_pLysS_AG__uid30681.cpf';
my $OUTPUT_3LINE = 1;

if ($#ARGV == -1)
{
  print "\nUsage: codon_pair_score.pl [-p codon_pair_preference_file] [-c coding_regions_file] [-w wt_seq] <input_file>\n\n";
  print "  <input_file>: input file in fasta format (Put '-' for STDIN)\n";
  print "  -c coding_regions_file: Import coding regions from file (Default: whole input)\n";
  print "  -p codon pair preference file: Import codon pair preferences\n";
  print "\n";
  exit;
}

$args{f} = $ARGV[0];

# FASTA record
$fasta_record =
{
  contig_count => 0,
  contig_name => [],
  contig_seq => [],
  contig_len => []
};

# AA record
$aa_record =
{
  aa => "na",		# Amino acid symbol
  num_of_codons => 0,	# Number of codons for this AA
  num_per_type => [],	# how many of each type of codon we have in our sequence
  location_list => [],	# locations of the codons for this AA in our sequence
  location_type => [],	# type of codon for each location
  at_least_two_different => 0,	# greater than 0 if the AA has at least two different codons in our sequence
  codon_map => [],
  perm => [],
  rev_perm => [],
  new_type => []
};

$aa_byname{$aa_record->{aa}} = $aa_record; # Create a hash of aa records

my $file_name = $args{f};
if ($file_name =~ /\/([^\/\.]+)\.fasta$/)
{
  $file_name = $1;
}
elsif ($file_name =~ /.*\/([^\/]+)$/)
{
  $file_name = $1;
}

open(SEQ_FILE, $args{f})
  or die "Could not open input_file ", $args{f}, " for reading: $!\n";
&importFASTA(*SEQ_FILE, $fasta_record);
close(SEQ_FILE);

$shuffled_sequence = $fasta_record->{contig_seq}->[0];

fasta_record # sequence record

$num_of_coding_regions = 0;
if ($args{c})
{
  &import_coding_regions;
}
else
{
  $num_of_coding_regions = 1;
  $coding_regions[0][0] = 0;
  $coding_regions[0][1] = $fasta_record->{contig_len}->[0];
}

&import_codon_info;
&create_work_space;

if ($args{p})
{
  &import_codon_pair_info;
}
elsif(not $args{c})
{
  $args{p} = $DEFAULT_CODON_PAIR_INFO_FILE;
  &import_codon_pair_info;
}

my ($cp_score, $cp_score_vs_len, $codon_length) = &calculate_codon_pair_score($shuffled_sequence);
#print "Length in codons = $codon_length\n";
if ($OUTPUT_3LINE)
{
  print "$file_name\tCPS\t", $cp_score, "\t";
  print "CPSpL\t", $cp_score_vs_len, "\n";
}
else
{
  print "\nCPS\t", $cp_score, "\t";
  print "CPSpL\t", $cp_score_vs_len, "\n";
}
# Codon pair structure
$codon_pair_record =
{
  bases => "na",	# The six bases that consist the codon pair
  AAs => "na",		# Two Amino acid symbols
  num => 0,	 	# Number of codon pairs appearing
  exp_num => 0,		# Number of codon pairs expected
  value => 0,		# value for this codon pair (log of ration observed/predicted with discounting)
  chisq => 0		# chi square value for observed and expected
};

$codon_pair{$codon_pair_record->{bases}} = $codon_pair_record; # Create a hash of codon pair records

my $B;
my $Bl;
my $jp_factor;
my $jp_summant;

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function that imports the codon pair information from the codon pair file and
# stores it in a hash.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def import_codon_pair_info():
{
  my $total_codon_pairs = 0;
  open (CODON_PAIR_FILE, $args{p})
    or die "Could not open codon pair info file ", $args{p}, " for reading: $!\n";
#  my $garbage = <CODON_PAIR_FILE>;
  while(<CODON_PAIR_FILE>)
  {
    chop($_);
    /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    $total_codon_pairs += $4;
  }
  close(CODON_PAIR_FILE);

# Calculate the Jeffreys-Perks law constants

  $B = $powers_of_4[6];
  $Bl = $B/2;
  $jp_factor = $total_codon_pairs / ($total_codon_pairs + $Bl);
  $jp_summant = 1/(2*($total_codon_pairs + $Bl));

  open (CODON_PAIR_FILE, $args{p})
    or die "Could not open codon pair info file ", $args{p}, " for reading: $!\n";
#  $garbage = <CODON_PAIR_FILE>;
  while(<CODON_PAIR_FILE>)
  {
    chop($_);
    /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    if (!$aa_pair_freq{$1})
    {
      $aa_pair_freq{$1} = $4;
    }
    else
    {
      $aa_pair_freq{$1} += $4;
    }
    my $first_codon = substr($2, 0, 3);
    my $second_codon = substr($2, 3, 3);
    if (!$codon_freq{$first_codon})
    {
      $codon_freq{$first_codon} = $4/2;
    }
    else
    {
      $codon_freq{$first_codon} += $4/2;
    }
    if (!$codon_freq{$second_codon})
    {
      $codon_freq{$second_codon} = $4/2;
    }
    else
    {
      $codon_freq{$second_codon} += $4/2;
    }
    my $first_aa = substr($1, 0, 1);
    my $second_aa = substr($1, 1, 1);
    $codon_to_aa_code{$first_codon} = $first_aa;
    $codon_to_aa_code{$second_codon} = $second_aa;
    if (!$aa_freq{$first_aa})
    {
      $aa_freq{$first_aa} = $4/2;
    }
    else
    {
      $aa_freq{$first_aa} += $4/2;
    }
    if (!$aa_freq{$second_aa})
    {
      $aa_freq{$second_aa} = $4/2;
    }
    else
    {
      $aa_freq{$second_aa} += $4/2;
    }
    $codon_pair{$2}->{AAs} = $1;
    my $discounted_num = &discount($4 + 0);
#    print STDERR "$4\t$discounted_num\n";
    $codon_pair{$2}->{num} = $discounted_num;
    my $exp_num = $3;
  }
  # Calculate the expected number of codon pairs and the codon pair values
  @codon_pair_scores = ();
  $min_codon_pair_score = 1000000;
  $max_codon_pair_score = -1000000;
  @codon_pair_chisq_values = ();
  $min_codon_pair_chisq_value = 1000000;
  $max_codon_pair_chisq_value = -1000000;
  foreach my $key (keys %codon_pair)
  {
    my $first_codon = substr($key, 0, 3);
    my $second_codon = substr($key, 3, 3);
    my $first_aa = substr($codon_pair{$key}->{AAs}, 0, 1);
    my $second_aa = substr($codon_pair{$key}->{AAs}, 1, 1);
    my $aa_pair = $codon_pair{$key}->{AAs};
    $codon_pair{$key}->{exp_num} = ($codon_freq{$first_codon}/$aa_freq{$first_aa}) *
      ($codon_freq{$second_codon}/$aa_freq{$second_aa}) * $aa_pair_freq{$codon_pair{$key}->{AAs}};
    my $value = log($codon_pair{$key}->{num}/$codon_pair{$key}->{exp_num});
    $codon_pair{$key}->{value} = $value;
    my $value_chisq = &chisq($codon_pair{$key}->{num}, $codon_pair{$key}->{exp_num});
    $codon_pair{$key}->{chisq} = $value_chisq;
    push (@codon_pair_scores, $value);
    push (@codon_pair_chisq_values, $value_chisq);
    if ($value > $max_codon_pair_score) {$max_codon_pair_score = $value;}
    if ($value < $min_codon_pair_score) {$min_codon_pair_score = $value;}
    if ($value_chisq > $max_codon_pair_chisq_value) {$max_codon_pair_chisq_value = $value_chisq;}
    if ($value_chisq < $min_codon_pair_chisq_value) {$min_codon_pair_chisq_value = $value_chisq;}
    if ($aa_pair_value_sum{$aa_pair})
    {
      $aa_pair_value_sum{$aa_pair} += $codon_pair{$key}->{value};
    }
    else
    {
      $aa_pair_value_sum{$aa_pair} = $codon_pair{$key}->{value};
    }
  }
}

def discount(arg, jp_factor, jp_summant):
  return arg*jp_factor + jp_summant

def import_codon_info(CODON_FILE):
    # Imports the codon information from the codons file and stores,
    # it in a pandas dataframe
    return pd.read_csv(CODON_FILE)

def import_coding_regions(CODING_REGION_FILE, format="fasta", **kwargs):
    # Function that imports the coding regions from the coding region file
    pairs = []
    with open(CODING_REGION_FILE, 'r+') as fh:
        for record in SeqIO.parse(fh, format, **kwargs):
            pairs.extend(re.findall('......', seq))

    return collections.counter(pairs)

def calculate_codon_pair_score(seq, ref):
    cp_initial_score = 0

    codon_pairs = re.findall('......', seq)
    for thispair in codon_pairs:
        cp_initial_score += ref[thispair]

    pair_num = len(codon_pairs)
    return (cp_initial_score, cp_initial_score/pair_num, pair_num)
