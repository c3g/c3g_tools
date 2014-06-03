#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
scriptName.pl

PURPOSE:
Takes x into input and generates y as output.

INPUT:
--infile <string> : Sequence file
				
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Your Name - your.name@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
  'infile=s'  => \$infile,
  'verbose'   => \$verbose,
  'help'      => \$help
);
if ($help) { print $usage; exit; }

## MAIN
my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
	$counter++;
}

