# #!/usr/bin/perl -w

package Common;
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter();

@ISA = qw(Exporter);
@EXPORT = qw(	importGENBANK 
		importFASTA 
		hamming 
		corr 
		num_to_str 
		str_to_num 
		mean 
		varr 
		stddev 
		chisq
		create_fasta_record 
		powers_of_4 
		export_fasta_sequences 
		edit_dist 
		array_min 
		array_max
		fisher_exact_test
		z_value
		p_value);

# These types of records will have to be created in programs 
# using the importFASTA and importGENBANK routines

# genbank record
# my $genbank_record =
# {
#   name => "STDIN",
#   segment_count => 0,
#   segment_name => [],
#   segment_seq => [],
#   segment_len => []
# };

# FASTA record
# my $fasta_record =
# {
#   name => "STDIN",
#   contig_count => 0,
#   contig_name => [],
#   contig_seq => [],
#   contig_length => []
# };
#
# To use:
#
# open(SEQ_FILE, "<some.file")
#   or die "Could not open input_file some.file for reading: $!\n";
# &importFASTA(*SEQ_FILE, $fasta_record);
# close(SEQ_FILE);
#

my @bases = ('A', 'C', 'G', 'T');
my %base_values = ( "A" => 0, "C" => 1, "G" => 2, "T" => 3 );
my @powers_of_4;
$powers_of_4[0] = 1;
for (my $i = 1; $i < 16; $i++)
{
  push (@powers_of_4, $powers_of_4[$i-1]<<2);
}

# FASTA record to be used implicitly by calling program
my $fr =
{
  name => "",
  contig_count => 0,
  contig_name => [],
  contig_seq => [],
  contig_length => []
};

my %fs;
$fs{$fr->{name}} = $fr;
my $fs_counter = 0;

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The subrouting create_fasta_record does what its name suggests, basically
# imports a fasta file and saves it in a fasta_record structure, which it returns
# relevant information in a structure
# Arguments:
# 1. Filehandle
# Returns:
# 1. fasta_record (fr) reference
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub create_fasta_record
{
  my ($file_ptr) = @_;
  $fs{$fs_counter}->{name} = $fs_counter;
  &importFASTA($file_ptr, $fs{$fs_counter});
  return $fs{$fs_counter++};
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The subroutine export_fasta_sequences is exporting a number of sequences from
# a fasta_record, the index of which is between $min and $max. The ids for the
# fasta sequences exported are the same as they were imported.
# Arguments:
# 1. fasta_record pointer
# 2. min index for exporting
# 3. max index for exporting
# 4. file pointer where to export the data (assumed to be valid)
# Returns:
#    Nothing important
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub export_fasta_sequences
{
  my ($fasta_ptr, $min, $max, $file_ptr) = @_;

  for (my $i = $min; $i <= $max; $i++)
  {
    print $file_ptr ">", $fasta_ptr->{contig_name}->[$i], "\n";
    print $file_ptr $fasta_ptr->{contig_seq}->[$i], "\n";
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The subrouting importGENBANK imports a GENBANK file, storing the contigs and
# relevant information in a structure
# Arguments:
# 1. Filehandle
# 2. Record pointer to store the data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub importGENBANK
{
  my $filehandle = $_[0];
  my $genbank_record = $_[1];
  my $segment_count = -1;
  my $seq = "";
  my $definition = "";
  my $len = 0;
  my $mode = "header";
  while(<$filehandle>)
  {
    if (/^\/\//)
    {
      if ($mode eq "sequence")
      {
        $seq =~ y/[U]/[T]/;
        push @{$genbank_record->{segment_seq}}, $seq;
        $seq = "";
        $definition = "";
      }
      $mode = "header";
    }
    elsif ($mode eq "sequence")
    {
      s/[^ACTGU]//gi;
      $seq .= $_;
    }
    elsif (/^LOCUS\s*\S*\s*(\d*)\s.*$/)
    {
      $segment_count++;
      push @{$genbank_record->{segment_len}}, $1 + 0;
    }
    elsif (/^ORIGIN/)
    {
      $mode = "sequence";
    }
    elsif (/^DEFINITION  (.*)$/)
    {
      $mode = "definition";
      $definition .= $1;
    }
    elsif ($mode eq "definition")
    {
      if (/^ {12}(.*)/)
      {
        $definition .= " $1";
      }
      else
      {
        $genbank_record->{segment_name}->[$segment_count] = $definition;
        $mode = "header";
      }
    }
  }
  $genbank_record->{segment_count} = ++$segment_count;
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The subrouting importFASTA imports a FASTA file, storing the contigs and
# relevant information in a structure
# Arguments:
# 1. Filehandle
#
# A second argument of the title of the record maybe added if multiple file
# input will be allowed.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub importFASTA
{
  my $filehandle = $_[0];
  my $fasta_record = $_[1];
  my $contig_count = -1;
  my $seq = "";
  my $len = 0;
  while(<$filehandle>)
  {
    if (!/^>/)
    {
      s/\s//gi;
      $seq .= $_;
    }
    else
    {
      /^>(.*)$/;
      push @{$fasta_record->{contig_name}}, $1;
      if ($contig_count >= 0)
      {
        push @{$fasta_record->{contig_len}}, length($seq);
	$seq =~ tr/a-z/A-Z/;
        push @{$fasta_record->{contig_seq}}, $seq;
      }
      $seq = "";
      $contig_count++;
    }
  }
  push @{$fasta_record->{contig_len}}, length($seq);
  $seq =~ tr/a-z/A-Z/;
  push @{$fasta_record->{contig_seq}}, $seq;
  $fasta_record->{contig_count} = ++$contig_count;
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The hamming procedure calculates the hamming distance of two sequences, which
# basically is the number of "bits" (in our case nucleotides) that need to be
# changed so the one sequence resembles the other.
# Arguments:
# 1. First sequence to be compared
# 2. Second sequence     -//-
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub hamming
{
  (my $seq1, my $seq2) = @_;
  my $len = length($seq1);
  my $hamming_distance = 0;
  for (my $i = 0; $i < $len; $i++)
  {
    if (substr($seq1, $i, 1) ne substr($seq2, $i, 1))
    {
      $hamming_distance++;
    }
  }
  return $hamming_distance;
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The subroutine 'mean' will calculate the mean of the elements of a list.
# Arguments:
# 1. Pointer to list
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub mean
{
  my ($arr) = @_;
  my $num_of_elements = @$arr;
  my $sum = 0;
  for (my $i = 0; $i < $num_of_elements; $i++)
  {
    $sum += $arr->[$i];
  }
  return ($sum/$num_of_elements);
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The subroutine 'varr' will calculate the variance of the elements of a list.
# Arguments:
# 1. Pointer to list
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub varr
{
  my ($arr) = @_;
  my $num_of_elements = @$arr;
  my $sum = 0;
  my $mean = &mean($arr);
  for (my $i = 0; $i < $num_of_elements; $i++)
  {
    $sum += ($arr->[$i] - $mean)**2;
  }
  if ($num_of_elements > 1)
  {
    return ($sum/($num_of_elements-1));
  }
  else
  {
    return 0;
  }
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Given an array of numeric values, return the standard deviation of the values.
## Arguments:
## 1. Pointer to an array of numeric values (pointer)
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub stddev
{
  return (sqrt(&varr($_[0])));
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Given two values, namely the observed and expected values of some quantity
# calculates the chi square statistic.
# Arguments:
# 1. Observed value
# 2. Expected value
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub chisq
{
  return (($_[0] - $_[1])**2)/$_[1];
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The subroutine 'corr' will calculate the Pearson's correlation of two lists.
# Arguments:
# 1. List # 1
# 2. List # 2
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub corr
{
  my ($arr1, $arr2, $min) = @_;
  if (@$arr1 < @$arr2) {$min = @$arr1} else {$min = @$arr2}
  my $sx = my $sy = my $sx2 = my $sy2 = my $sxy = 0;
  for (my $i = 0; $i < $min; $i++)
  {
    $sx += $arr1->[$i];
    $sy += $arr2->[$i];
    $sx2 += ($arr1->[$i])*($arr1->[$i]);
    $sy2 += ($arr2->[$i])*($arr2->[$i]);
    $sxy += ($arr1->[$i])*($arr2->[$i]);
  }
  return ($sxy - ($sx*$sy/$min))/sqrt(($sx2 - ($sx*$sx/$min))*($sy2 - ($sy*$sy/$min)));
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The next subroutine takes as input a number and a size and returns a
# corresponding string of that size composed by DNA bases.
# Arguments:
#   1. Number to be converted to DNA string
#   2. Desired size of string
# Note: 1st_argument <= 2nd_argument**4
# Definitions:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub num_to_str
{
  my $num;

  $num = $_[0];
  my $return_value = '';
  for (my $j = 0; $j < $_[1]; $j++)
  {
    $return_value = $bases[$num % 4].$return_value;
    $num = $num>>2;
  }
  $return_value;
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The str_to_num subroutine takes as input a string returns a corresponding
# number, such that there is a one-to-one mapping from the strings based on
# the {A, C, G, T} alphabet and the natural numbers.
# Arguments:
#   1. String to be converted/indexed to a number
# Definitions:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub str_to_num
{
  my $str = $_[0];
  my $return_value = 0;
  for (my $j = 0; $j < length($str); $j++)
  {
    my $char = substr($str, $j, 1);
    if (exists($base_values{$char}))
    {
      $return_value += $base_values{substr($str, $j, 1)}*$powers_of_4[length($str) - $j - 1];
    }
    else
    {
      return -1;
    }
  }
  return $return_value;
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The subroutine edit_dist computes the edit distance between two strings.
# Arguments:
# 1. String #1
# 2. String #2
# Returned value: Edit distance
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub edit_dist
{
#  my $str1 = $_[0];
#  my $str2 = $_[1];
  my @str1 = unpack("C*", $_[0]);
  my @str2 = unpack("C*", $_[1]);
  my $len1 = length($_[0]);
  my $len2 = length($_[1]);
#  my $len1 = length($str1);
#  my $len2 = length($str2);

  # We have to be a little careful with the index of the table we are using
  # Since we do not want to reconstruct the alignment, only two rows of the
  # matrix are sufficient. Obviously we will use pointers to these rows, so
  # we can exchange them freely without copying elements.

  # Initialize old_row
  my $old_row;
  for (my $i = 0; $i <= $len1; $i++)
  {
    $old_row->[$i] = $i;
  }
  my $new_row;
  for (my $i = 1; $i <= $len2; $i++)
  {
    $new_row = [];
    $new_row->[0] = $i;
    for (my $j = 1; $j <= $len1; $j++)
    {
      $new_row->[$j] = &min3($old_row->[$j] + 1,
                             $new_row->[$j-1] + 1,
                             $old_row->[$j-1] + (($str1[$j-1] eq $str2[$i-1]) ? 0 : 1));
#                             $old_row->[$j-1] + ((substr($str1, $j-1, 1) eq substr($str2, $i-1, 1)) ? 0 : 1));
    }
    $old_row = $new_row;
  }
  return $new_row->[$len1];
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sub min3
{
  my $min_temp = ($_[0] < $_[1]) ? $_[0] : $_[1];
  return ($_[2] < $min_temp) ? $_[2] : $min_temp;
}

sub array_min
{
  my ($arr) = @_;
  my $min = $arr->[0];
  foreach my $elem (@$arr)
  {
    if ($elem < $min) { $min = $elem; }
  }
  return $min;
}

sub array_max
{
  my ($arr) = @_;
  my $max = $arr->[0];
  foreach my $elem (@$arr)
  {
    if ($elem > $max) { $max = $elem; }
  }
  return $max;
}

###############################################################################
# The function z_value calculates the standardized variable z for the 
# sampling distribution of differences in proportions. The inputs are:
# 1. P1 -> Sample proportion 1
# 2. P2 -> Sample proportion 2
# 3. N1 -> Large sample size 1
# 4. N2 -> Large sample size 2
###############################################################################
sub z_value
{
  my ($P1, $P2, $N1, $N2) = @_;
  if ($P1 < $P2)
  {
    ($P1, $P2, $N1, $N2) = ($P2, $P1, $N2, $N1);
  }
  my $p = (($N1*$P1) + ($N2*$P2))/($N1+$N2);
  my $q = 1 - $p;
#  print "P1 = $P1, P2 = $P2, N1 = $N1, N2 = $N2, p = $p, q = $q\n";
  my $stddev = sqrt($p*$q*((1/$N1)+(1/$N2)));
  return ($P1-$P2)/$stddev;
}

###############################################################################
# The function p_value returns the probability value for the one tail test
# that corresponds to the z value given as an argument. Basically it gives the
# probability that we can have such a high or higher z-value. The calculation
# is taken from the javascript code found at:
# http://www.fourmilab.ch/rpkp/experiments/analysis/zCalc.html
# The following are comments found in that code:
#     /*  The following JavaScript functions for calculating normal and
#         chi-square probabilities and critical values were adapted by
#         John Walker from C implementations
#         written by Gary Perlman of Wang Institute, Tyngsboro, MA
#         01879.  Both the original C code and this JavaScript edition
#         are in the public domain.  */
#
#     /*  POZ  --  probability of normal z value
#
#         Adapted from a polynomial approximation in:
#                 Ibbetson D, Algorithm 209
#                 Collected Algorithms of the CACM 1963 p. 616
#         Note:
#                 This routine has six digit accuracy, so it is only useful for absolute
#                 z values <= 6.  For z values > to 6.0, poz() returns 0.0.
#     */
# 
# Arguments:
# 1. z value
sub p_value
{
  my ($z) = @_;
  my ($x, $y, $w) = (0.0, 0.0, 0.0);
  my $Z_MAX = 6.0;
  if ($z == 0.0) 
  {
    $x = 0.0;
  } 
  else 
  {
    $y = 0.5 * abs($z);
    if ($y > ($Z_MAX * 0.5)) 
    {
      $x = 1.0;
    } 
    elsif ($y < 1.0) 
    {
      $w = $y * $y;
      $x = ((((((((0.000124818987 * $w
        - 0.001075204047) * $w + 0.005198775019) * $w
        - 0.019198292004) * $w + 0.059054035642) * $w
        - 0.151968751364) * $w + 0.319152932694) * $w
        - 0.531923007300) * $w + 0.797884560593) * $y * 2.0;
    } 
    else 
    {
      $y -= 2.0;
      $x = (((((((((((((-0.000045255659 * $y
        + 0.000152529290) * $y - 0.000019538132) * $y
        - 0.000676904986) * $y + 0.001390604284) * $y
        - 0.000794620820) * $y - 0.002034254874) * $y
        + 0.006549791214) * $y - 0.010557625006) * $y
        + 0.011630447319) * $y - 0.009279453341) * $y
        + 0.005353579108) * $y - 0.002141268741) * $y
        + 0.000535310849) * $y + 0.999936657524;
    }
  }
  return 1 - ($z > 0.0 ? (($x + 1.0) * 0.5) : ((1.0 - $x) * 0.5));
}

sub fisher_exact_test
{
  my ($n1, $n2, $n3, $n4) = @_;
  my $p_value = 0.0;

  open(INPUT_FILE, ">temp.R")
    or die "Could not open file temp.R for writing: $!\n";
  print INPUT_FILE "x <- cbind(c($n1, $n2), c($n3, $n4))\n";
  if ($n1 > $n2)
  {
    print INPUT_FILE "y <- fisher.test(x, alt=\"g\")\n";
  }
  else
  {
    print INPUT_FILE "y <- fisher.test(x, alt=\"l\")\n";
  }
  print INPUT_FILE "mylist <- (\"y\")\n";
  print INPUT_FILE "dump(mylist, file=\"temp.Rout\")\n";
  close(INPUT_FILE);
  `R CMD BATCH < temp.R`;
  open(OUTPUT_FILE, "<temp.Rout")
    or die "Could not open file temp.Rout for reading: $!\n";
  while (<OUTPUT_FILE>)
  {
    if (/p.value = ([^,]*),/)
    {
      $p_value = $1 + 0.0;
    }
  }
  close(OUTPUT_FILE);
  return $p_value;
}


1;
