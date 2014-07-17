#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

## MUGQIC LIBS:
use SampleSheet;

my $usage=<<'ENDHERE';
NAME:
getMiSeqBarcodes.pl

PURPOSE:

INPUT:
--runId <string>       : MiSeq run ID. ex: M00833_0173
--sampleSheet <string> : Sample sheet
		
OUTPUT:
STDOUT                 : Fasta file of barcode sequences.

NOTES:
Please load the following modules in order to run:
mugqic/perl/5.18.2
mugqic/tools/1.9
mugqic/pipeline/1.3

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $infile, $runId, $sampleSheet);
my $verbose = 0;

GetOptions(
  'infile=s'      => \$infile,
  'runId=s'       => \$runId,
  'sampleSheet=s' => \$sampleSheet,
  'verbose'       => \$verbose,
  'help'          => \$help
);
if ($help) { print $usage; exit; }

## MAIN
#system("module load mugqic/perl/5.18.2");
#system("module load mugqic/tools/1.9");
#system("module load mugqic/pipeline/1.3");

# First get sample names from sample sheet

my $rHoAoH_sampleInfo = SampleSheet::parseSampleSheetAsHash($sampleSheet);
my %sampleNamesHash;
for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
  $sampleName =~ s/\.//g;
  $sampleNamesHash{$sampleName} = $sampleName;
}

# Then get path to SampleSheet.nanuq.csv
# Hack to add _\d+ at the end of $root
my @dirname;
push(@dirname, "/lb/robot/miSeqSequencer/miSeqRuns/");
push(@dirname, "/lb/robot/miSeqSequencer/miSeqRuns/2013");
push(@dirname, "/lb/robot/miSeqSequencer/miSeqRuns/2014");
push(@dirname, "/lb/robot/miSeqSequencer/miSeqRuns/2015");
push(@dirname, "/lb/robot/miSeqSequencer/miSeqRuns/2016");

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
  my $headerName = $row[2];
  my $sequence = $row[4];
  print STDERR $headerName."\n";
  for my $key (keys %sampleNamesHash){
    print STDERR $key."\n";
    
    if($headerName =~ m/($key)/){
	    print STDOUT ">".$1."\n".$row[4]."\n";
      delete $sampleNamesHash{$key};
    }
  }
}
close(IN);
exit;
