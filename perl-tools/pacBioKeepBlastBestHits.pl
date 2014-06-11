#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
pacBioKeepBlastBestHits.pl

PURPOSE:
Takes keeps best n hits from blast table.

INPUT:
--infile <string> : Sequence file
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $infile, $n);
my $verbose = 0;

GetOptions(
  'infile=s'  => \$infile,
  'n=i'       => \$n,
  'verbose'   => \$verbose,
  'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $header = "";
my %hash;
my $i=1;
my $j=0;

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
  chomp;
  if($_ =~ m/^#/){
    $header .= $_."\n";
    next;
  }

  my $id = (split /\t/, $_)[0];
  

  if(exists $hash{$i."===".$id}){
    $hash{$i."===".$id}{$j} = $_;
    next if($n == $j+1);
    $j++;

  }else{
    $hash{$i."===".$id}{$j} = $_;
    $j = 0;
    $i = 1;

  }  
}

print STDOUT $header;

for my $k1 (sort keys %hash) { 
  for my $k2 (sort keys %{ $hash{$k1} }){
    print STDOUT $hash{$k1}{$k2}."\n";

  }
}
exit;
