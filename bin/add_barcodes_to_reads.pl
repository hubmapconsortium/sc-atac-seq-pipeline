#!/usr/bin/perl -w

use strict;
use warnings;
use 5.010;
use Getopt::Long qw(GetOptions);


my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';
$rcTable{'N'}='N';

sub revComp{
  my $seq = shift;
  my $rcSeq='';
  for(my $i=0; $i<length($seq); $i++){
    $rcSeq = $rcTable{uc(substr($seq,$i,1))} . $rcSeq;
  }
  return $rcSeq;
}

my $read1_file;
my $read2_file;
my $read3_file;
my $prefix;
GetOptions(
	'input-fastq1=s' => \$read1_file,
	'input-fastq2=s' => \$read3_file,
	'input-barcode-fastq=s' => \$read2_file,
	'output-fastq-prefix=s' => \$prefix)
	or die "Usage: $0 --input-fastq1 NAME\n";


say $read1_file;
say $read2_file;
say $read3_file;
say $prefix;


#my $read1_file = $ARGV[0];
#my $read2_file = $ARGV[1];
#my $read3_file = $ARGV[2];
#my $prefix = $ARGV[3];

say 'opening input files';

open(R1, "zcat -f $read1_file |") or die("error reading $read1_file\n");
open(R2, "zcat -f $read2_file |") or die("error reading $read2_file\n");
open(R3, "zcat -f $read3_file |") or die("error reading $read3_file\n");

say 'opening output files';

open(R1_OUT, ">>$prefix.R1.fastq") or die("error writing to $prefix.R1.fastq");
open(R3_OUT, ">>$prefix.R3.fastq") or die("error writing to $prefix.R3.fastq");

while(my $r1_1 = <R1>){

  my $r1_2 = <R1>;
  my $r1_3 = <R1>;
  my $r1_4 = <R1>;
  my $r2_1 = <R2>;
  my $r2_2 = <R2>;
  my $r2_3 = <R2>;
  my $r2_4 = <R2>;
  my $r3_1 = <R3>;
  my $r3_2 = <R3>;
  my $r3_3 = <R3>;
  my $r3_4 = <R3>;

  my @fields = split / /, $r1_1;

  my $barcode = substr($r2_2, 0, 16);

  $r1_1=~ s/@//g;
  $r3_1=~ s/@//g;

  print R1_OUT "@", $barcode, "::", $r1_1, $r1_2, $r1_3, $r1_4;
  print R3_OUT "@", $barcode, "::", $r3_1, $r3_2, $r3_3, $r3_4;
}

