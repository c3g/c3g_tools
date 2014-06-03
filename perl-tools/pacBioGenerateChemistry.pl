#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Find;

my $usage=<<'ENDHERE';
NAME:
pacBioGenerateChemistry.pl

PURPOSE:
Loop in a directory, find every .h5 file and
write the path in a file.

INPUT:
--fofn <string>      : fofn
--chemistry <string> : P4-C2 or P5-C3				

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $fofn, $chemistry);
my $verbose = 0;

GetOptions(
  'fofn=s'      => \$fofn,
  'chemistry=s' => \$chemistry,
  'verbose'     => \$verbose,
  'help'        => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my %hash;
open(IN, "<".$fofn) or die "Can't open $fofn\n";
while(<IN>){
  chomp;
  $_ =~ m/(\S\d+_\d+_\d+_\S\d+_\S\d_\S\d)/;
  $hash{$1} = $_;  
}
close(IN);

print STDOUT "<Map>\n";
for my $key ( keys %hash ) {
  print STDOUT "   <Mapping><Movie>$key</Movie><SequencingChemistry>$chemistry</SequencingChemistry></Mapping>\n";
}
print STDOUT "</Map>\n";

exit;

