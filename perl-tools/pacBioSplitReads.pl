#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use List::Util qw(sum);
use File::Slurp;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
pacBioAssemblyStats.pl

PURPOSE:
Generate relevant plots and table(s) of current PacBio assembly.

INPUT:
--infile <string>       : infile (subreads)
--cutoff <int>          : read length cutoff

OUTPUT:
--outfileShort <string> : Short reads
--outfileLong <string>  : Long reads

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $infile, $cutoff, $outfileShort, $outfileLong);
my $verbose = 0;

GetOptions(
    'infile=s' 	     		=> \$infile,
    'cutoff=i' 	     		=> \$cutoff,
	'outfileShort=s'        => \$outfileShort,
	'outfileLong=s'         => \$outfileLong,
    'verbose' 	     		=> \$verbose,
    'help' 		     		=> \$help
);
if ($help) { print $usage; exit; }

## Validate
die "--infile missing\n" unless($infile);
die "--cutoff missing\n" unless($cutoff);
die "--outfileShort missing\n" unless($outfileShort);
die "--outfileLong missing\n" unless($outfileLong);

## MAIN
open(OUT_S, ">".$outfileShort) or die "Can't open $outfileShort\n";
open(OUT_L, ">".$outfileLong) or die "Can't open $outfileLong\n";

my $readCount = 0;
my $shortReads = 0;
my $longReads = 0;

if($infile =~ m/\.fastq|\.fq/){
	my $ref_fastq_db = Iterator::FastqDb->new($infile) or die("Unable to open Fasta file, $infile\n");
	while( my $curr = $ref_fastq_db->next_seq() ) {
		if(length($curr->seq()) >= $cutoff){
			print OUT_L $curr->header()."\n".$curr->seq()."\n";	
			$longReads++;
		}else{
			print OUT_S $curr->header()."\n".$curr->seq()."\n";	
			$shortReads++;
		}
		$readCount++;	
	}
}elsif($infile =~ m/\.fasta|\.fa|\.fsa|\.fna/){
	my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
	while( my $curr = $ref_fasta_db->next_seq() ) {
		if(length($curr->seq()) >= $cutoff){
			print OUT_L $curr->header()."\n".$curr->seq()."\n";	
			$longReads++;
		}else{
			print OUT_S $curr->header()."\n".$curr->seq()."\n";	
			$shortReads++;
		}
		$readCount++;	
	}
}
close(OUT_S);
close(OUT_L);

print STDERR "Processed ".$readCount." reads...\nShort reads: ".$shortReads."\tLong reads: ".$longReads."\n";

exit;
