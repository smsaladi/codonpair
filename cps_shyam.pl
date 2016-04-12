#!/usr/bin/perl -w

# Utility to calculate the codon pair score of a coding sequence
#
# Author: Dimitris Papamichail
#
# Version 0.1
#
# Creation: 2005-11-20
# 
#  Copyright (c) 2008 Research Foundation of the State University of
#  New York. All rights reserved.
#
#  Redistribution and use of the source code, with or without
#  modification, are permitted provided that the following conditions are met:
#
#  1.      Redistributions must reproduce the above copyright notice, this
#  list of conditions and the following disclaimer in the  documentation
#  and/or other materials provided with the distribution.  Redistributions of
#  source code must also reproduce this information in the source code itself.
#
#  2.      If the program is modified, redistributions must include a notice
#  (in the same places as above) indicating that the redistributed program is
#  not identical to the version distributed by the original creator.
#
#  3.      The name of the original creator may not be used to endorse or
#  promote products derived from this software without specific prior written
#  permission.
#
#  We also request that use of this software be cited in publications as
#
#     J. R. Coleman, D. Papamichail, S. Skiena, B. Futcher, E. Wimmer
#     and S. Mueller, "Virus attenuation by genome-scale changes in
#     codon-pair bias: a novel method for developing viral vaccines"
#     Science, Vol. 320, 1784-1787, 27 June 2008
#
#  THIS SOFTWARE IS PROVIDED BY THE ORIGINAL AUTHOR ``AS IS'' AND  ANY
#  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  ARE
#  DISCLAIMED. IN NO EVENT SHALL THE ORIGINAL AUTHOR BE LIABLE  FOR ANY
#  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL  DAMAGES
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS  OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)  HOWEVER
#  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
#  OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
#  SUCH DAMAGE.
#

use lib "/temp/saladi/cps/lib";
use Common;
use Getopt::Std;

my @powers_of_4;
$powers_of_4[0] = 1;
for (my $i = 1; $i < 16; $i++)
{
  push (@powers_of_4, $powers_of_4[$i-1]<<2);
}

my $MAX_SEQUENCE_LENGTH = 1000000;
my $CODON_FILE_LOCATION = '/temp/saladi/cps/codons';
my $DEFAULT_CODON_PAIR_INFO_FILE = '/temp/saladi/cps/Escherichia_coli__BL21_Gold_DE3_pLysS_AG__uid30681.cpf';
my $OUTPUT_3LINE = 1;

&getopts("w:o:l:c:r:d:p:i:hs", \%args); # -v, -D, -o ARG, sets $args{v}, $args{D}, $args{o}

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

if ($fasta_record->{contig_count} > 1)
{
  die "\nCurrent version does not support more than one fasta input sequence\n\n";
}

if ($fasta_record->{contig_len}->[0] > $MAX_SEQUENCE_LENGTH)
{
  die "\nMaximum sequence length supported is $MAX_SEQUENCE_LENGTH. Limit exceeded.\n\n";
}

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
sub import_codon_pair_info
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

sub discount
{
  return $_[0]*$jp_factor + $jp_summant;
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function that imports the codon information from the codons file and stores
# it in a hash.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub import_codon_info
{
  open (CODON_FILE, "<$CODON_FILE_LOCATION")
    or die "Could not open codon info file codons for reading!\n";
  %aa = ();
  while(<CODON_FILE>)
  {
    chop($_);
    /^([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/;
    push (@{$aa{$2}}, $1);
    $codon_to_aa{$1} = $2;
  }
  close(CODON_FILE);
  $variable_aas = 0;	# AAs that have more than one codon representations
  foreach $key (keys %aa)
  {
    push (@aas, $key);
    $aa_codon_num{$key} = 0;
    foreach $codon (@{$aa{$key}})
    {
      $codon_type{$codon} = $aa_codon_num{$key};
      $aa_codon_num{$key}++;
    }
    if (($aa_codon_num{$key} > 1) && ($key ne "Ter"))
    {
      $variable_aas++;
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function that imports the coding regions from the coding region file $args{c}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub import_coding_regions
{
  open (CODING_REGION_FILE, $args{c})
    or die "Could not open coding regions file ", $args{c}, " for reading!\n";
  $all_region_seq_length = 0;
  while(<CODING_REGION_FILE>)
  {
    chop($_);
    /([^-]+)-([^-]+)$/;
    if ((($2 + 1 - $1) % 3) != 0)
    {
      die "\n ERROR: Coding region $1-$2 not a multiple of 3!\n\n";
    }
    if (($1+0) > ($2+0))
    {
      die "\nERROR: Coding region $1-$2 is reversely indicated. Please correct.\n\n";
    }
    ($coding_regions[$num_of_coding_regions][0], $coding_regions[$num_of_coding_regions][1]) = ($1-1, $2-1);
    my $reg_length = $coding_regions[$num_of_coding_regions][1] - $coding_regions[$num_of_coding_regions][0] + 1;
    $all_region_seq_length += $reg_length;
    $region_length[$num_of_coding_regions++] = $region_length;
  }
  close(CODING_REGION_FILE);
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function that isolates the codons we will work with. Basically creates
# a table with pointers to all codons we can change in the coding regions
# specified, after excluding the locked regions.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub create_work_space
{
  for (my $i = 0; $i < $num_of_coding_regions; $i++)
  {
    for (my $j = $coding_regions[$i][0]; $j < $coding_regions[$i][1]; $j +=3)
    {
      my $codon = substr($fasta_record->{contig_seq}->[0], $j, 3);
      my $codon_t = $codon_type{$codon};
      my $aa = $codon_to_aa{$codon};
      $aa_byname{$aa}->{num_of_codons}++;
      push (@{$aa_byname{$aa}->{location_list}}, $j);
      push (@{$aa_byname{$aa}->{location_type}}, $codon_t);
    }
  }
  foreach my $aa (@aas)
  {
    if (!$aa_byname{$aa}->{num_of_codons})
    {
      $aa_byname{$aa}->{num_of_codons} = 0;
    }
  }
}

sub calculate_codon_pair_score
{
  my ($seq) = @_;
  $cp_initial_score = 0;
  my $pair_num = 0;
  my $seq_print = "";
  if ($OUTPUT_3LINE)
  {
    print $file_name, "\t";
    $seq_print .= $file_name . "\t";
  }
  for (my $i= 0; $i < $num_of_coding_regions; $i++)
  {
    for (my $j = $coding_regions[$i][0]; $j < $coding_regions[$i][1] - 3; $j += 3)
    {
      my $cp = substr($seq, $j, 6);
      if (exists $codon_pair{$cp})
      {
        $cp_initial_score += $codon_pair{$cp}->{value};
	$pair_num++;

        if ($OUTPUT_3LINE)
        {
          $seq_print .= substr($cp,0,3) . "\t";
          print $codon_pair{$cp}->{value}, "\t";
        }
        else
        {
          print substr($cp,0,3), "\t", $codon_pair{$cp}->{value}, "\n";
        }
      }
      else
      {
        if ($OUTPUT_3LINE)
        {
          $seq_print .= substr($cp,0,3) . "\t";
          print "-1\t";
        }
        else
        {
          print substr($cp,0,3), "\t-1\n";
        }
        print STDERR "Codon pair $cp does not exist!\n";
      }
    }
  }
  if ($OUTPUT_3LINE)
  {
    print "\n",$seq_print,"\n";
  }
  return ($cp_initial_score, $cp_initial_score/$pair_num, $pair_num);
}
