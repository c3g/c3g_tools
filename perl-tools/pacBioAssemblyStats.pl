#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use List::Util qw( sum min max);
use File::Slurp;

my $usage=<<'ENDHERE';
NAME:
pacBioAssemblyStats.pl

PURPOSE:
Generate relevant plots and table(s) of current PacBio assembly.

INPUT:
--filteredSummary <string>       : File path (usually /bla/bla/filtered_summary.csv)
--assemblyQc <string>            : Celera assembly QC stats file
--contigs <string>               : Celera contigs.
--sampleName <string>            : sample name
--suffix <string>                : suffix
--estimatedGenomeSize <int>      : Estimated genome size

OUTPUT:
--outdir <string>  : ouput directory where table files will be written.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $outdir, $sampleName, $suffix, $estimatedGenomeSize, $smrtCells, $filteredSummary, $assemblyQc, $contigs);
my $verbose = 0;

GetOptions(
	'filteredSummary=s'       => \$filteredSummary,
	'assemblyQc=s'            => \$assemblyQc,
	'contigs=s'               => \$contigs,
    'outdir=s' 	     		  => \$outdir,
	'sampleName=s'   		  => \$sampleName,
	'suffix=s'			 	  => \$suffix,
	'estimatedGenomeSize=i'   => \$estimatedGenomeSize,
	'smrtCells=i'			  => \$smrtCells,
    'verbose' 	     		  => \$verbose,
    'help' 		     		  => \$help
);
if ($help) { print $usage; exit; }

## Validate
die "--filteredSummary missing\n" unless($filteredSummary);
die "--assemblyQc missing\n" unless($assemblyQc);
die "--outdir missing\n" unless($outdir);

## MAIN

my $cmd;

# Generate plots of raw data.
#$cmd .= " Rscript \$R_TOOLS/pacBioAssemblyPlots.R -i ".$indir."/filtering/data/filtered_summary.csv -o ".$outdir;
$cmd .= " Rscript \$R_TOOLS/pacBioAssemblyPlots.R -i ".$filteredSummary." -o ".$outdir;
print STDERR "[DEBUG] ".$cmd."\n";
system($cmd);
die "Command failed: $!\n" if($? != 0);

# Get assembly stats
my $totalContigs;    #		TotalContigsInScaffolds
my $totalBases;      #  	TotalBasesInScaffolds
my $minContigLength; #		MinContigLength
my $maxContigLength; #		MaxContigLength
my $N50Bases;        #		N50ContigBases
my $contigCoverage;  #		$contigCoverage = ContigsOnly
my $gcContent;       #		Content

#my $contigStats = "$indir/$suffix/assembly/9-terminator/$sampleName.qc";
my $contigStats = $assemblyQc;
open(IN, '<'.$contigStats) or die "Can't open $contigStats\n";
while(<IN>){
	chomp;
	if($_ =~ m/TotalContigsInScaffolds=(.*)$/){
		$totalContigs = $1;	
	}
	if($_ =~ m/TotalBasesInScaffolds=(.*)$/){
		 $totalBases = $1;
	}
	if($_ =~ m/MinContigLength=(.*)$/){
		$minContigLength = $1;
	}
	if($_ =~ m/MaxContigLength=(.*)$/){
		$maxContigLength = $1;
	}
	if($_ =~ m/N50ContigBases=(.*)$/){
		$N50Bases = $1;
	}
	if($_ =~ m/ContigsOnly=(.*)$/){
		$contigCoverage = $1;
	}
	if($_ =~ m/Content=(.*)$/){
		$gcContent = $1;
	}
}
close(IN);

# Print to file relevant values of the assembly process.
open(OUT, '>'.$outdir."/summaryTableAssembly.tsv") or die "Can't open ".$outdir."/summaryTableAssembly.tsv";
print OUT "\"Description\"\t\"Value\"\n";
print OUT "\"Sample name\"\t\"$sampleName\"\n";
print OUT "\"Assembly name\"\t\"$sampleName-$suffix\"\n";
print OUT "\"Estimated genome size (bp)\"\t\"$estimatedGenomeSize\"\n";
print OUT "\"Total contigs\"\t\"$totalContigs\"\n";
print OUT "\"Total bases in contigs (bp)\"\t\"$totalBases\"\n";
print OUT "\"Minimum contig length (bp)\"\t\"$minContigLength\"\n";
print OUT "\"Maximum contig length (bp)\"\t\"$maxContigLength\"\n";
print OUT "\"N50 contigs bases\"\t\"$N50Bases\"\n";
print OUT "\"Contigs coverage (X)\"\t\"$contigCoverage\"\n";
print OUT "\"GC content (%)\"\t\"$gcContent\"\n";

# Print to file relevant values from raw reads.
# Compute N25, N50, N75 and display length of each contig.

$totalBases = 0;
my $counter = 1;

my %hash;

my $ref_fasta_db = Iterator::FastaDb->new($contigs) or die("Unable to open Fasta file, $contigs\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
	my $length = length($curr->seq);
	my $header = $curr->header;
	$header =~ s/>//;

	#print STDOUT "Sequence: ".$header."\t".$length." bp\n";
	$totalBases += $length;
	$hash{$header} = $length;
	$counter++;
}
$counter = ($counter - 1);
#print OUT "Total of ".$counter." sequences\n";

foreach my $key (sort {$hash{$b} <=> $hash{$a}} (keys %hash)){
	print OUT "\"".$key . " length (bp)\"\t\"" .$hash{$key}."\"\n";	
}

my ($N25, $N50, $N75);
my $cummSum = 0;
$counter = 1;
foreach my $key (sort {$hash{$b} <=> $hash{$a}} (keys %hash)){
	$cummSum += $hash{$key};
	my $ratio = int( ($cummSum / $totalBases) * 100);

	if($ratio >= 25 && !$N25){
		print OUT "\"N25 - 25% of total sequence length is contained in the ".$counter." sequence(s) having a length (bp) >= \"\t\"".$hash{$key}."\"\n";
		$N25=1;
	}
	if($ratio >= 50 && !$N50){
		print OUT "\"N50 - 50% of total sequence length is contained in the ".$counter." sequence(s) having a length (bp) >= \"\t\"".$hash{$key}."\"\n";
		$N50=1;
	}
	if($ratio >= 75 && !$N75){
		print OUT "\"N75 - 75% of total sequence length is contained in the ".$counter." sequence(s) having a length (bp) >= \"\t\"".$hash{$key}."\"\n";
		$N75=1;
	}
	$counter++;
}
close(OUT);

# Raw reads/filtered reads stats.
if( (-e $filteredSummary) and (-s $filteredSummary) ){
	open(OUT, '>'.$outdir."/summaryTableReads.tsv") or die "Can't open ".$outdir."/summaryTableReads.tsv";
	print OUT "\"Description\"\t\"Value\"\n";
	#my $readsStats = "$indir/filtering/results/filterReports_filterStats.xml";
	my $readsStats = $filteredSummary;

	my $totalReads = 0;
	my $passedReads = 0;
	my $failedReads = 0;
	my @readLength;
	my @passedReads;

	open(IN, '<'.$readsStats) or die "Can't open ".$readsStats."\n";
	while(<IN>){
		chomp;
		next if($. == 1);
		my @row = split(/,/, $_);
		my $readLength = $row[3];					
		push(@row, $readLength);
		if($row[8] == 1){
			$passedReads++;
			push(@passedReads, $readLength);
		}
		$failedReads++ if($row[8] == 0);
		$totalReads++;
		push(@readLength, $readLength);
	}
		
	my $shortestRead = min(@passedReads);
	my $longestRead = max(@passedReads);
	my $averageReadLength = sum(@readLength)/@readLength;
	my $averageReadLengthQCpassed = sum(@passedReads)/@passedReads;
	
	$averageReadLength = int $averageReadLength;
	$averageReadLengthQCpassed = int $averageReadLengthQCpassed;

	print OUT "\"Number of SMRT cells\"\t\"".$smrtCells."\"\n";
	print OUT "\"Total subreads\"\t\"".$totalReads."\"\n";
	print OUT "\"Subreads that passed QC\"\t\"".$passedReads."\"\n";
	print OUT "\"Subreads that failed QC\"\t\"".$failedReads."\"\n";
	print OUT "\"Average subread length\"\t\"".$averageReadLength."\"\n";
	print OUT "\"Average subread length that passed QC\"\t\"".$averageReadLengthQCpassed."\"\n";
	print OUT "\"Shortest subread length that passed QC\"\t\"".$shortestRead."\"\n";
	print OUT "\"Longest subread length that passed QC\"\t\"".$longestRead."\"\n";
	close(IN);
	close(OUT);
}

exit;


