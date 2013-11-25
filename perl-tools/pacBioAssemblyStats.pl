#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use List::Util qw(sum);
use File::Slurp;

my $usage=<<'ENDHERE';
NAME:
pacBioAssemblyStats.pl

PURPOSE:
Generate relevant plots and table(s) of current PacBio assembly.

INPUT:
--filteredSubreadsTable <string> : File path (usually /bla/bla/filtered_summary.csv)
--filteredSummary <string>       : File path (usually /bla/bla/filterReports_filterStats.xml)
--assemblyQc <string>            : Celera assembly QC stats file
--sampleName <string>            : sample name
--suffix <string>                : suffix
--estimatedGenomeSize <int>      : Estimated genome size

OUTPUT:
--outdir <string>  : ouput directory

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $indir, $outdir, $sampleName, $suffix, $estimatedGenomeSize, $smrtCells, $filteredSubreadsTable, $filteredSummary, $assemblyQc);
my $verbose = 0;

GetOptions(
    'indir=s' 	     		  => \$indir,
	'filteredSubreadsTable=s' => \$filteredSubreadsTable,
	'filteredSummary=s'       => \$filteredSummary,
	'assemblyQc=s'            => \$assemblyQc,
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
die "--indir missing\n" unless($indir);
die "--outdir missing\n" unless($outdir);

## MAIN

my $cmd;

# Generate plots of raw data.
#$cmd .= " Rscript \$R_TOOLS/pacBioAssemblyPlots.R -i ".$indir."/filtering/data/filtered_summary.csv -o ".$outdir;
$cmd .= " Rscript \$R_TOOLS/pacBioAssemblyPlots.R -i ".$filteredSubreadsTable." -o ".$outdir;
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
close(OUT);

# Print to file relevant values from raw reads.

open(OUT, '>'.$outdir."/summaryTableReads.tsv") or die "Can't open ".$outdir."/summaryTableReads.tsv";
print OUT "\"Description\"\t\"Value\"\n";
#my $readsStats = "$indir/filtering/results/filterReports_filterStats.xml";
my $readsStats = $filteredSummary;
open(IN, '<'.$readsStats) or die "Can't open ".$readsStats."\n";
while(<IN>){
	chomp;
	
	if($_ =~ m/name=\"(Pre-Filter Polymerase Read Bases)\" value=\"(\d+)\"/){
		print OUT "\"$1 (bp)\"\t\"$2\"\n";
	}
	if($_ =~ m/name=\"(Post-Filter Polymerase Read Bases)\" value=\"(\d+)\"/){
		print OUT "\"$1 (bp)\"\t\"$2\"\n";
	}
	if($_ =~ m/name=\"(Pre-Filter Polymerase Reads)\" value=\"(\d+)\"/){
		print OUT "\"$1\"\t\"$2\"\n";
	}
	if($_ =~ m/name=\"(Post-Filter Polymerase Reads)\" value=\"(\d+)\"/){
		print OUT "\"$1\"\t\"$2\"\n";
	}
	if($_ =~ m/name=\"(Pre-Filter Polymerase Read Length)\" value=\"(\d+)\"/){
		print OUT "\"$1 (bp)\"\t\"$2\"\n";
	}
	if($_ =~ m/name=\"(Post-Filter Polymerase Read Length)\" value=\"(\d+)\"/){
		print OUT "\"$1 (bp)\"\t\"$2\"\n";
	}
	if($_ =~ m/name=\"(Pre-Filter Polymerase Read Quality)\" value=\"(\d+)\"/){
		print OUT "\"$1\"\t\"$2\"\n";
	}
	if($_ =~ m/name=\"(Post-Filter Polymerase Read Quality)\" value=\"(\d+)\"/){
		print OUT "\"$1\"\t\"$2\"\n";
	}
}
print OUT "\"Number of SMRT cells\"\t\"".$smrtCells."\"\n";
close(IN);
close(OUT);


exit;


