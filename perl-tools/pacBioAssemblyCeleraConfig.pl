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
--infile <string>    : Sequence file
--merylThreads       : 
--frgMinLen          : 
--overlapper         : 
--merCompression     : 
--merSize            : 
--merylMemory        : Set the merylMemory parameter...
--ovlErrorRate       : 
--ovlMinLen          :
--ovlStoreMemory     :
--ovlConcurrency     :
--ovlCorrConcurrency :
--cnsConcurrency     :
--frgCorrThreads     :
--ovlErrorRate       :
--ovlThreads         : 

TODO
			
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $infile, $merylThreads, $frgMinLen, $frgCorrThreads, $cnsConcurrency, $ovlCorrConcurrency, $ovlConcurrency, $ovlStoreMemory, $ovlMinLen, $ovlErrorRate, $ovlThreads, $overlapper, $merCompression, $merSize, $merylMemory);
my $verbose = 0;

GetOptions(
    'infile=s' 	   		    => \$infile,
	'merylThreads=i'		=> \$merylThreads,
	'ovlThreads=i'		    => \$ovlThreads,
	'frgCorrThreads=i'		=> \$frgCorrThreads,
	'cnsConcurrency=i'		=> \$cnsConcurrency,
	'ovlCorrConcurrency=i'	=> \$ovlCorrConcurrency,
	'ovlConcurrency=i'		=> \$ovlConcurrency,
	'ovlStoreMemory=i'		=> \$ovlStoreMemory,
	'ovlErrorRate=f'		=> \$ovlErrorRate,
	'ovlMinLen=i'           => \$ovlMinLen,
	'frgMinLen=i'	        => \$frgMinLen,
	'overlapper=s' 		    => \$overlapper,
	'merCompression=i'  	=> \$merCompression,
	'merSize=i'			    => \$merSize,
	'merylMemory=i'		    => \$merylMemory,
    'verbose' 	   		    => \$verbose,
    'help' 		   		    => \$help
);
if ($help) { print $usage; exit; }

## Validate
die "--infile missing\n" 			unless($infile);
warn "--merylMemory missing\n"		unless($merylMemory);
warn "--merylThreads missing\n"     unless($merylThreads);
warn "--frgMinLen missing\n" 		unless($frgMinLen);
warn "--overlapper missing\n" 		unless($overlapper);
warn "--merCompression missing\n" 	if(!defined($merCompression));

## MAIN
#my $text = read_file( $infile );
open(IN, '<'.$infile) or die "Can't open $infile\n";
while(<IN>){
	chomp;
	if($_ =~ m/overlapper=/ && defined($overlapper)){
		print STDOUT "overlapper=$overlapper\n";

	}elsif($_ =~ m/frgMinLen=/ && defined($frgMinLen)){
		print STDOUT "frgMinLen=$frgMinLen\n"; 

	}elsif($_ =~ m/merSize=/ && defined($merSize)){
		print STDOUT "merSize=$merSize\n";

	}elsif($_ =~ m/ovlThreads=/ && defined($ovlThreads)){
		print STDOUT "ovlThreads=$ovlThreads\n";

	}elsif($_ =~ m/merCompression=/ && defined($merCompression)){
		print STDOUT "merCompression=$merCompression\n";

	}elsif($_ =~ m/frgCorrThreads=/ && defined($frgCorrThreads)){
		print STDOUT "frgCorrThreads=$frgCorrThreads\n";

	}elsif($_ =~ m/merylThreads=/ && defined($merylThreads)){
		print STDOUT "merylThreads=$merylThreads\n";

	#}elsif($_ =~ m/merOverlapperThreads=/ && defined($)){
	#	print STDOUT "merOverlapperThreads=$num_threads\n";
    #
	}elsif($_ =~ m/cnsConcurrency=/ && defined($cnsConcurrency)){
		print STDOUT "cnsConcurrency=$cnsConcurrency\n";

	}elsif($_ =~ m/merylMemory=/ && defined($merylMemory)){
		print STDOUT "merylMemory=$merylMemory\n";

	}elsif($_ =~ m/ovlCorrConcurrency=/ && defined($ovlCorrConcurrency)){
		print STDOUT "ovlCorrConcurrency=$ovlCorrConcurrency\n";

	}elsif($_ =~ m/ovlConcurrency=/ && defined($ovlConcurrency)){
		print STDOUT "ovlConcurrency=$ovlConcurrency\n";

	}elsif($_ =~ m/ovlStoreMemory=/ && defined($ovlStoreMemory)){
		print STDOUT "ovlStoreMemory=$ovlStoreMemory\n";

	}elsif($_ =~ m/ovlErrorRate=/ && defined($ovlErrorRate)){
		print STDOUT "ovlErrorRate=$ovlErrorRate\n";

	}elsif($_ =~ m/ovlMinLen=/ && defined($ovlMinLen)){
		print STDOUT "ovlMinLen=$ovlMinLen\n";

	}else{ 
		print STDOUT $_."\n";

	}
}
close(IN);

exit;
