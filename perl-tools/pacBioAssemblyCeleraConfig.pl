#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use List::Util qw(sum);
use File::Slurp;

my $usage=<<'ENDHERE';
NAME:
pacBioAssemblyCeleraConfig.pl

PURPOSE:

INPUT:
--infile <string>  : Sequence file
--coverage <int>   : Coverage
--genomeSize <int> : Estimated genome size	

--coverageCutoff   : To get cutoff based on cummulative length
                     of sequences up until a (Coverage * 
                     Genome size).
--xml <string>     : XML file to parse.
--xmlOut <string>  : XML outfile with modified value in the minReadLength field.
--merylMemory      : Set the merylMemory parameter...

TODO
--hgapCutoff       : To get cutoff based upon a more 
                     sophisticated algorithm.
			
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $infile, $num_threads, $minReadSize, $overlapper, $merCompression, $merSize, $merylMemory);
my $verbose = 0;

GetOptions(
    'infile=s' 	   		=> \$infile,
	'num_threads=i'		=> \$num_threads,
	'minReadSize=i'		=> \$minReadSize,
	'overlapper=s' 		=> \$overlapper,
	'merCompression=i' 	=> \$merCompression,
	'merSize=i'			=> \$merSize,
	'merylMemory=i'		=> \$merylMemory,
    'verbose' 	   		=> \$verbose,
    'help' 		   		=> \$help
);
if ($help) { print $usage; exit; }

## Validate
die "--infile missing\n" 			unless($infile);
warn "--merylMemory missing\n"		unless($merylMemory);
warn "--num_threads missing\n" 		unless($num_threads);
warn "--minReadSize missing\n" 		unless($minReadSize);
warn "--overlapper missing\n" 		unless($overlapper);
warn "--merCompression missing\n" 	if(!defined($merCompression));

## MAIN
#my $text = read_file( $infile );
open(IN, '<'.$infile) or die "Can't open $infile\n";
while(<IN>){
	chomp;
	if($_ =~ m/overlapper=/ && defined($overlapper)){
		print STDOUT "overlapper=$overlapper\n";

	}elsif($_ =~ m/frgMinLen=/ && defined($minReadSize)){
		print STDOUT "frgMinLen=$minReadSize\n"; 

	}elsif($_ =~ m/merSize=/ && defined($merSize)){
		print STDOUT "merSize=$merSize\n";

	}elsif($_ =~ m/ovlThreads=/ && defined($num_threads)){
		print STDOUT "ovlThreads=$num_threads\n";

	}elsif($_ =~ m/merCompression=/ && defined($merCompression)){
		print STDOUT "merCompression=$merCompression\n";

	}elsif($_ =~ m/frgCorrThreads=/ && defined($num_threads)){
		print STDOUT "frgCorrThreads=$num_threads\n";

	}elsif($_ =~ m/merylThreads=/ && defined($num_threads)){
		print STDOUT "merylThreads=$num_threads\n";

	}elsif($_ =~ m/merOverlapperThreads=/ && defined($num_threads)){
		print STDOUT "merOverlapperThreads=$num_threads\n";

	}elsif($_ =~ m/cnsConcurrency=/ && defined($num_threads)){
		print STDOUT "cnsConcurrency=$num_threads\n";

	}elsif($_ =~ m/merylMemory=/ && defined($num_threads)){
		print STDOUT "merylMemory=$merylMemory\n";

	}else{ 
		print STDOUT $_."\n";

	}
}
close(IN);

exit;
