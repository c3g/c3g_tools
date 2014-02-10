#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
getBarcodes.pl

PURPOSE:

INPUT:
--runId <string> : MiSeq run ID. ex: M00833_0173
		
OUTPUT:
STDOUT           : Fasta file of barcode sequences.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $infile, $runId);
my $verbose = 0;

GetOptions(
    'infile=s' 	=> \$infile,
	'runId=s'   => \$runId,
    'verbose' 	=> \$verbose,
    'help' 		=> \$help
);
if ($help) { print $usage; exit; }

## MAIN

# Hack to add _\d+ at the end of $root
my @dirname;
push(@dirname, "/lb/robot/miSeqSequencer/miSeqRuns/");
push(@dirname, "/lb/robot/miSeqSequencer/miSeqRuns/2013");
push(@dirname, "/lb/robot/miSeqSequencer/miSeqRuns/2014");
push(@dirname, "/lb/robot/miSeqSequencer/miSeqRuns/2015");

my $root;
foreach my $dirname (@dirname){

	next if(!-d $dirname);

	#print STDERR "dirname: ".$dirname."\n";

	opendir(DIR, $dirname);
	my @files = readdir(DIR);
	closedir DIR;

	#140129_M00833_0173_000000000-A7RED
	foreach my $file (@files){
		my $searchString = "(".$runId.")";
		if($file =~ m/$searchString/){
			$root = "$dirname/$file/SampleSheet.nanuq.csv";
			last;
		}
	}
}

print STDERR "Found " . $root . " file.\n";

open(IN, "<".$root) or die "Can't open file ".$root."\n";
while(<IN>){
	chomp;
	next if $. == 1;
	my @row = split(/,/, $_);
	print STDOUT ">".$row[2]."\n".$row[4]."\n";
}
close(IN);
exit;
