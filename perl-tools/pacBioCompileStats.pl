#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use List::Util qw( sum min max);
use File::Slurp;
use Iterator::FastaDb;
use File::Basename;
use File::Find;

my $usage=<<'ENDHERE';
NAME:
pacBioAssemblyStats.pl

PURPOSE:
Generate relevant plots and table(s) of current PacBio assembly.

INPUT:
--indir <string>    : Directory containing sample's assemblies.

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $indir);
my $verbose = 0;

GetOptions(
	'indir=s'                 => \$indir,
    'verbose' 	     		  => \$verbose,
    'help' 		     		  => \$help
);
if ($help) { print $usage; exit; }

## Validate
die "--indir missing\n" unless($indir);

## Global scope variable
my %hQc;
my %hFasta;
my $sampleName = fileparse($indir);
print STDERR "[DEBUG] sampleName: $sampleName\n";

## SUBS
sub eachFile{
	my $filename = $_;
	my $fullpath = $File::Find::name;
	#remember that File::Find changes your CWD, 
	#so you can call open with just $_

	if (-e $filename) { 
		
		if($fullpath =~ m/9-terminator\/$sampleName\S+\.qc$/){
			if($fullpath =~ m/(\d+percent)_(merSize\d+)/){
				print STDERR "[DEBUG] $fullpath\n";
				$hQc{"$1_$2"} = $fullpath;
			}else{
				die "Can't find proper nomenclature on .qc files.\n";	
			}
		}
		#if($fullpath =~ m/9-terminator\/$sampleName\S+\.ctg.fasta$/){
		#	if($fullpath =~ m/(\d+percent)_(merSize\d+)/){
		#		print STDERR "[DEBUG] $fullpath\n";
		#		$hFasta{"$1_$2"} = $fullpath;
		#	}else{
		#		die "Can't find proper nomenclature on .ctg.fasta files.\n";	
		#	}
		#}
	}
}

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text
}

## MAIN

# Find .qc and ctg.fasta files.
find (\&eachFile, $indir);

my %hResults;

foreach my $key (sort{$a cmp $b} keys %hQc){
	print STDERR "[DEBUG] ".$key."\t".$hQc{$key}."\n";
	my $file = $hQc{$key};

	# Get assembly stats
	my $totalContigs;    #		TotalContigsInScaffolds
	my $totalBases;      #  	TotalBasesInScaffolds
	my $minContigLength; #		MinContigLength
	my $maxContigLength; #		MaxContigLength
	my $N50Bases;        #		N50ContigBases
	my $N25Bases;        #		N25ContigBases
	my $N75Bases;        #		N75ContigBases
	my $contigCoverage;  #		$contigCoverage = ContigsOnly
	my $gcContent;       #		Content
	
	open(IN, '<'.$file) or die "Can't open $file\n";
	while(<IN>){
		chomp;
		if($_ =~ m/TotalContigsInScaffolds=(.*)$/){
			$totalContigs = $1;	
		}
		if($_ =~ m/TotalBasesInScaffolds=(.*)$/){
			 $totalBases = commify($1);
		}
		if($_ =~ m/MinContigLength=(.*)$/){
			$minContigLength = commify($1);
		}
		if($_ =~ m/MaxContigLength=(.*)$/){
			$maxContigLength = commify($1);
		}
		if($_ =~ m/N50ContigBases=(.*)$/){
			$N50Bases = commify($1);
		}
		if($_ =~ m/N25ContigBases=(.*)$/){
			$N25Bases = commify($1);
		}
		if($_ =~ m/N75ContigBases=(.*)$/){
			$N75Bases = commify($1);
		}
		if($_ =~ m/ContigsOnly=(.*)$/){
			$contigCoverage = $1."X";
		}
		if($_ =~ m/Content=(.*)$/){
			$gcContent = $1;
		}
	}
	close(IN);

	$hResults{$key}{'totalContigsInScaffolds'} = $totalContigs;
	$hResults{$key}{'totalBasesInScaffolds'}   = $totalBases;
	$hResults{$key}{'minContigLength'}         = $minContigLength;
	$hResults{$key}{'maxContigLength'}         = $maxContigLength;
	$hResults{$key}{'N25'}                     = $N25Bases;
	$hResults{$key}{'N50'}                     = $N50Bases;
	$hResults{$key}{'N75'}                     = $N75Bases;
	$hResults{$key}{'contigs coverage'}             = $contigCoverage;
	$hResults{$key}{'gc%'}                      = $gcContent;	
}	

my @string;
$string[0] = "#Assembly name\t";
# Then loop through result hash and print results in a table
my $i = 1;
foreach my $key (sort{$a cmp $b} keys %hResults) {
	$string[$i] = $key;	

    foreach my $key2 (sort{$a cmp $b} keys %{ $hResults{$key} }) {
		if($i == 1){
			$string[0] .= "\t".$key2;
		}
		$string[$i] .= "\t".$hResults{$key}{$key2}
	}
	$i++;
}

foreach my $line (@string){
	print STDOUT "$line\n";
}

exit;

