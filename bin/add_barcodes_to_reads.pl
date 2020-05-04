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
  
  #my $barcode4_and_1 = $fields[1];
  #$barcode4_and_1 =~ s/1:N:0://g;
  #chomp($barcode4_and_1);
  #my $barcode4 = substr($barcode4_and_1, 0, 8);
  #if($illumina_workflow eq "B"){
  #  $barcode4 = revComp($barcode4);
  #}
  #my $barcode1 = substr($barcode4_and_1, 8, 8);
  
  my $barcode1 = substr($r2_2, 0, 8);
  my $barcode2 = substr($r2_2, 38, 8);
  my $barcode3 = substr($r2_2, 76, 8);
  my $umi = substr($r2_2, 84, 10);
  
  $r1_1=~ s/@//g;
  $r3_1=~ s/@//g;
  
  #my $cell_barcode = $barcode1 . $barcode2 . $barcode3 . $barcode4;
  my $cell_barcode = $barcode1 . $barcode2 . $barcode3;
  my $rc_cell_barcode = revComp($cell_barcode);

  # Don't add the prefix to the read id string because SnapTools
  # snap-pre will fail to construct a correct snap file.
  # This probably should be debugged in snap_pre.py 
  #print R1_OUT "@", $prefix, "_", $rc_cell_barcode, ":$umi:", $r1_1, $r1_2, $r1_3, $r1_4;
  #print R3_OUT "@", $prefix, "_", $rc_cell_barcode, ":$umi:", $r3_1, $r3_2, $r3_3, $r3_4;

  print R1_OUT "@", $rc_cell_barcode, ":$umi:", $r1_1, $r1_2, $r1_3, $r1_4;
  print R3_OUT "@", $rc_cell_barcode, ":$umi:", $r3_1, $r3_2, $r3_3, $r3_4;


}

