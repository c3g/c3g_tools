#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
addSampleNamesToVcf.pl

PURPOSE:
Takes x into input and generates y as output.

INPUT:
--infile <string>      : Vcf file
--sampleNames <string> : text file with one sample name per line.
	
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Your Name - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $infile, $sampleNames);
my $verbose = 0;

GetOptions(
  'infile=s'      => \$infile,
  'sampleNames=s' => \$sampleNames,
  'verbose'       => \$verbose,
  'help'          => \$help
);
if ($help) { print $usage; exit; }

die "--infile <string> missing\n" unless($infile);
die "--sampleNames <string> missing\n" unless($sampleNames);

## MAIN

# Store sample name in array;
#
my @names;
open(SN, "<".$sampleNames) or die "Can't open $sampleNames\n";
while(<SN>){
  chomp;
  push(@names, $_);
}

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
  chomp;
  if($_ =~ m/#/){
    print $_."\n";
    next;
  }
  my @row = split(/\t/, $_);

  if($row[0] =~ m/Chrom/i){ #Header desc. line.
    print STDOUT "##HEADER For each sampleName columns, values are => Cons:Cov:Reads1:Reads2:Freq:P-value\n"; 
    print STDOUT join("\t", @row[0..9]);
    print STDOUT "\t";
    print STDOUT join("\t", @names);
    print STDOUT "\n";
    

  }else{
    my $sampleRow = $row[10];
    my @sampleRow = split(/\s/, $sampleRow);
    
    die "--sampleNames file contains @names entries while  vcf file contains @sampleRow samples...\n" if(@names != @sampleRow);
    
    print STDOUT join("\t", @row[0..9]);
    print "\t";
    print STDOUT join("\t", @sampleRow);
    print STDOUT "\n";
  
  }
}
exit;
