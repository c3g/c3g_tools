#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use List::Util qw( sum min max);
use File::Slurp;
use Iterator::FastaDb;
use POSIX;

my $usage=<<'ENDHERE';
NAME:
pacBioAssemblyStats.pl

PURPOSE:
Generate relevant plots and table(s) of current PacBio assembly.

INPUT:
--filteredSummary <string>       : File path (usually /bla/bla/filtered_summary.csv)
--contigs <string>               : Celera contigs.
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
my ($help, $outdir, $sampleName, $suffix, $estimatedGenomeSize, $smrtCells, $filteredSummary, $contigs);
my $verbose = 0;

GetOptions(
	'filteredSummary=s'       => \$filteredSummary,
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
die "--outdir missing\n" unless($outdir);
die("--filteredSummary file is empty or does not exists! (Typed wrong filename?)\n") if((!-e $filteredSummary) and (!-s $filteredSummary));	

## MAIN

my $cmd;

# Generate plots of raw data.
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
my $totalSequencedBases = 0;

#my $contigStats = "$indir/$suffix/assembly/9-terminator/$sampleName.qc";
#my $contigStats = $assemblyQc;
#open(IN, '<'.$contigStats) or die "Can't open $contigStats\n";
#while(<IN>){
#	chomp;
#	if($_ =~ m/TotalContigsInScaffolds=(.*)$/){
#		$totalContigs = $1;	
#	}
#	if($_ =~ m/TotalBasesInScaffolds=(.*)$/){
#		 $totalBases = $1;
#	}
#	if($_ =~ m/MinContigLength=(.*)$/){
#		$minContigLength = $1;
#	}
#	if($_ =~ m/MaxContigLength=(.*)$/){
#		$maxContigLength = $1;
#	}
#	if($_ =~ m/N50ContigBases=(.*)$/){
#		$N50Bases = $1;
#	}
#	if($_ =~ m/ContigsOnly=(.*)$/){
#		$contigCoverage = $1;
#	}
#	#if($_ =~ m/Content=(.*)$/){
#	#	$gcContent = $1;
#	#}
#}
#close(IN);


# Print to file relevant values from raw reads.

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
		if($row[7] == 1){
			$passedReads++;
            $totalSequencedBases += $readLength;					
			push(@passedReads, $readLength);
		}
		$failedReads++ if($row[7] == 0);
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

# Compute N25, N50, N75 and display length of each contig.
$totalBases = 0;
my $counter = 1;
my %hash;	
my $gcCount=0;
my @seqLengths;
my $opt_i = 100;
my %len= ();
my $n = 0;
my $int;
my $totalLength = 0;
my @length;

my $ref_fasta_db = Iterator::FastaDb->new($contigs) or die("Unable to open Fasta file, $contigs\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
	my $length = length($curr->seq);
	my $header = $curr->header;
	$header =~ s/>//;
    push(@length, $length);

	###################
	push @seqLengths, $length; # record length for N50 calc's 
	$n++; 
	$int = floor( $length/$opt_i );  
	$totalLength += $length; 
	if( !defined($len{$int}) ) { 
		$len{$int} = 1;  
	} else { 
		$len{$int}++; 
	}   
	$gcCount += ($curr->seq()  =~ tr/gGcC/gGcC/);
	###################

	#print STDOUT "Sequence: ".$header."\t".$length." bp\n";
	$totalBases += $length;
	$hash{$header} = $length;
	$counter++;
}
$maxContigLength = max @length;
$minContigLength = min @length;
$gcContent = sprintf "%.2f", ($gcCount/$totalBases * 100);
$counter = ($counter - 1);
print STDOUT "Total of ".$counter." sequences\n";
$totalContigs = $counter;

# Calculate N25, N50, N75, and N90 and counts 
my $N25; my $N50; my $N75; my $N90; 
my $N25count=0; my $N50count=0; my $N75count=0; my $N90count=0; 
my $frac_covered = $totalLength; 
@seqLengths = reverse sort { $a <=> $b } @seqLengths; 
$N25 = $seqLengths[0]; 
while ($frac_covered > $totalLength*3/4) { 
	$N25 = shift(@seqLengths); 
	$N25count++; $N50count++; $N75count++; $N90count++; 
	$frac_covered -= $N25; 
} 
$N50 = $N25; 
while ($frac_covered > $totalLength/2) { 
	$N50 = shift(@seqLengths); 
	$N50count++; $N75count++; $N90count++; 
	$frac_covered -= $N50; 
} 
$N75 = $N50; 
while ($frac_covered > $totalLength/4) { 
	$N75 = shift(@seqLengths); 
	$N75count++; $N90count++; 
	$frac_covered -= $N75; 
} 
$N90 = $N75; 
while ($frac_covered > $totalLength/10) { 
	$N90 = shift(@seqLengths); 
	$N90count++; 
	$frac_covered -= $N90; 
}

print STDOUT "N25 - 25% of total sequence length is contained in the ".$N25count." sequence(s) having a length >= ".$N25." bp\n";
print STDOUT "N50 - 50% of total sequence length is contained in the ".$N50count." sequence(s) having a length >= ".$N50." bp\n";
print STDOUT "N75 - 75% of total sequence length is contained in the ".$N75count." sequence(s) having a length >= ".$N75." bp\n";
print STDOUT "N90 - 90% of total sequence length is contained in the ".$N90count." sequence(s) having a length >= ".$N90." bp\n";

# Compute X coverage value used by hgap.
my $percent;
if($suffix =~ m/(\d+)percent/){
	$percent = $1;
}else{
	die "Something wrong with assembly name...\n";
}
my $estimatedCov = ($totalSequencedBases / $estimatedGenomeSize);
my $covCutoff = sprintf "%.2f", ($estimatedCov * $percent / 100);

# Print to file relevant values of the assembly process.
open(OUT, '>'.$outdir."/summaryTableAssembly.tsv") or die "Can't open ".$outdir."/summaryTableAssembly.tsv";
print OUT "\"Description\"\t\"Value\"\n";
print OUT "\"Sample name\"\t\"$sampleName\"\n";
print OUT "\"Assembly name\"\t\"$sampleName-$suffix\"\n";
print OUT "\"HGAP cutoff (X)\"\t\"$covCutoff\"\n";
print OUT "\"Estimated genome size (bp)\"\t\"$estimatedGenomeSize\"\n";
print OUT "\"Total contigs\"\t\"$totalContigs\"\n";
print OUT "\"Total bases in contigs (bp)\"\t\"$totalBases\"\n";
print OUT "\"Minimum contig length (bp)\"\t\"$minContigLength\"\n";
print OUT "\"Maximum contig length (bp)\"\t\"$maxContigLength\"\n";
print OUT "\"N25 contigs bases\"\t\"$N25\"\n";
print OUT "\"N50 contigs bases\"\t\"$N50\"\n";
print OUT "\"N75 contigs bases\"\t\"$N75\"\n";
print OUT "\"N90 contigs bases\"\t\"$N90\"\n";
#print OUT "\"Contigs coverage (X)\"\t\"$contigCoverage\"\n";
print OUT "\"GC content (%)\"\t\"$gcContent\"\n";
close(OUT);


#foreach my $key (sort {$hash{$b} <=> $hash{$a}} (keys %hash)){
#	#print STDOUT $key . "\t" .$hash{$key}."\n";	
#}
#my ($N25, $N50, $N75);
#my $cummSum = 0;
#$counter = 1;
#foreach my $key (sort {$hash{$b} <=> $hash{$a}} (keys %hash)){
#	$cummSum += $hash{$key};
#	my $ratio = int( ($cummSum / $totalBases) * 100);
#
#	if($ratio >= 25 && !$N25){
#		print STDOUT "N25 - 25% of total sequence length is contained in the ".$counter." sequence(s) having a length >= ".$hash{$key}." bp\n";
#		$N25=1;
#	}
#	if($ratio >= 50 && !$N50){
#		print STDOUT "N50 - 50% of total sequence length is contained in the ".$counter." sequence(s) having a length >= ".$hash{$key}." bp\n";
#		$N50=1;
#	}
#	if($ratio >= 75 && !$N75){
#		print STDOUT "N75 - 75% of total sequence length is contained in the ".$counter." sequence(s) having a length >= ".$hash{$key}." bp\n";
#		$N75=1;
#	}
#	$counter++;
#}

exit;


