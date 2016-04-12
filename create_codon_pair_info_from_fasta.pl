#!/usr/bin/perl -w

# Codon pair creator
#
# Copyright 2005 Dimitris Papamichail
#
# Version 0.1
#
# Creation: 2005-08-19
#
# Given a species' fasta information of the proteins, this utility
# will output the codon pair information

use lib "/temp/saladi/cps/lib";

use strict;
use Common;
use Getopt::Std;
use bignum;

my $DEBUG_CODON_IMPORT = 0;
my $CODON_FILE_LOCATION = '/temp/saladi/cps/codons';

# FASTA record
my $fasta_record =
{
  contig_count => 0,
  contig_name => [],
  contig_seq => [],
  contig_length => []
};

# Codon variables

my $variable_aas = 0;    # AAs that have more than one codon representations
my @aas = ();
my %codon_to_aa = ();
my %aa_code = ();
my %aa_codon_num = ();
my %codon_type = ();
my %type_to_codon = ();
my %codon_freq = ();
my %aa_freq = ();
my %aa_pair_freq = ();

&import_codon_info;

# Codon pair variables

my %codon_pair = ();

open(SEQ_FILE, "-")
  or die "Could not open fasta file for reading: $!\n";
&importFASTA(*SEQ_FILE, $fasta_record);
close(SEQ_FILE);

for (my $i = 0; $i < $fasta_record->{contig_count}; $i++)
{
  my $seq = $fasta_record->{contig_seq}->[$i];
  for (my $j = 0; $j <= length($seq) - 6; $j += 3)
  {
    my $cp = substr($seq, $j, 6);
    my $c1 = substr($cp, 0, 3);
    my $c2 = substr($cp, 3, 3);
    if ($c1 !~ /N/ and $c2 !~ /N/)
    {
      my $aa1 = $codon_to_aa{$c1};
      my $aa2 = $codon_to_aa{$c2};
      my $aa_pair = $aa1.$aa2;
      if (exists $codon_freq{$c1}) {$codon_freq{$c1}++;} else {$codon_freq{$c1}=1;}
      if (exists $codon_freq{$c2}) {$codon_freq{$c2}++;} else {$codon_freq{$c2}=1;}
      if (exists $aa_freq{$aa1}) {$aa_freq{$aa1}++;} else {$aa_freq{$aa1}=1;}
      if (exists $aa_freq{$aa2}) {$aa_freq{$aa2}++;} else {$aa_freq{$aa2}=1;}
      if (exists $aa_pair_freq{$aa_pair}) {$aa_pair_freq{$aa_pair}++;} else {$aa_pair_freq{$aa_pair}=1;}
#      print $codon_pair, "\n";
      if (exists $codon_pair{$cp})
      {
        $codon_pair{$cp}++;
      }
      else
      {
        $codon_pair{$cp} = 1;
      }
    }
  }
}

my @powers_of_4;
$powers_of_4[0] = 1;
for (my $i = 1; $i < 16; $i++)
{
  push (@powers_of_4, $powers_of_4[$i-1]<<2);
}

for (my $i = 0; $i < $powers_of_4[6]; $i++)
{
  my $cp = &num_to_str($i, 6);
  my $cd1 = substr($cp, 0, 3);
  my $cd2 = substr($cp, 3, 3);
  my $aa1 = $codon_to_aa{$cd1};
  my $aa2 = $codon_to_aa{$cd2};
  my $aa_pair = $aa1.$aa2;
  if ($aa1 && $aa2 && ($aa1 ne "Ter") && ($aa2 ne "Ter"))
  {
    print $aa_code{$aa1}, $aa_code{$aa2}, "\t";
    print "$cp\t";
    if (exists $codon_freq{$cd1} && exists $codon_freq{$cd2} && exists $aa_freq{$aa1} && exists $aa_freq{$aa2} && exists $aa_pair_freq{$aa_pair})
    { 
      my $exp = $codon_freq{$cd1}*$codon_freq{$cd2}/($aa_freq{$aa1}*$aa_freq{$aa2})*$aa_pair_freq{$aa_pair};
      $exp = sprintf("%.2f", $exp);
      print $exp, "\t";
    }
    else
    {
      print "0", "\t";
    }
    if (exists $codon_pair{$cp})
    {
      print $codon_pair{$cp}, "\n";
    }
    else
    {
      print "0\n";
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function that imports the codon information from the codons file and stores
# it in a hash.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub import_codon_info
{
  open (CODON_FILE, "<$CODON_FILE_LOCATION")
    or die "Could not open codon info file codons for reading: $!\n";
  my %aa = ();
  while(<CODON_FILE>)
  {
    chop($_);
    /^([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/;
    push (@{$aa{$2}}, $1);
    $codon_to_aa{$1} = $2;
    $aa_code{$2} = $3;
#    $code_to_aa{$3} = $2;
#  print "Now code $3 corresponds to three letter code $2\n";
  }
  close(CODON_FILE);
  foreach my $key (keys %aa)
  {
    push (@aas, $key);
    $aa_codon_num{$key} = 0;
    foreach my $codon (@{$aa{$key}})
    {
      $codon_type{$codon} = $aa_codon_num{$key};
      $type_to_codon{$key}->[$aa_codon_num{$key}] = $codon;
      $aa_codon_num{$key}++;
    }
    if (($aa_codon_num{$key} > 1) && ($key ne "Ter"))
    {
#      $aa_to_index{$key} = $variable_aas;      # Index number for AAs that have multiple codons
# print "AA $key has ", $aa_codon_num{$key}, " codon types\n";
#      $index_to_aa{$variable_aas} = $key;      # Translates an index to the corresponding AA
      $variable_aas++;
    }
  }
  if ($DEBUG_CODON_IMPORT)
  {
    foreach my $key (keys %aa)
    {
      print $key, "\t", $aa_code{$key}, "\n";
      foreach my $codon (@{$aa{$key}})
      {
         print "\t$codon";
         print "\t", $codon_to_aa{$codon};
         print "\t", $codon_type{$codon}, "\n";
      }
      print "\n";
      print "Codons of $key: ";
      for (my $i = 0; $i < $aa_codon_num{$key}; $i++)
      {
        print $type_to_codon{$key}->[$i], " ";
      }
      print "\n";
    }
  }
}
