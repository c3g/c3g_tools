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
--shortReads <string>            : File path to short reads
--longReads <string>             : File path to long reads.
--correctedReads <string>        : File path to corrected reads.
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
my ($help, $outdir, $sampleName, $suffix, $estimatedGenomeSize, $smrtCells, $filteredSummary, $contigs, $shortReads, $longReads, $correctedReads);
my $verbose = 0;

GetOptions(
	'filteredSummary=s'     => \$filteredSummary,
  'shortReads=s'          => \$shortReads,
  'longReads=s'           => \$longReads,
  'correctedReads=s'      => \$correctedReads,
	'contigs=s'             => \$contigs,
  'outdir=s' 	     		    => \$outdir,
	'sampleName=s'   		    => \$sampleName,
	'suffix=s'			 	      => \$suffix,
	'estimatedGenomeSize=i' => \$estimatedGenomeSize,
	'smrtCells=i'			      => \$smrtCells,
  'verbose' 	     		    => \$verbose,
  'help' 		     	        => \$help
);
if ($help) { print $usage; exit; }

## Validate
die "--filteredSummary missing\n" unless($filteredSummary);
die "--outdir missing\n" unless($outdir);
die("--filteredSummary file is empty or does not exists! (Typed wrong filename?)\n") if((!-e $filteredSummary) and (!-s $filteredSummary));	
die "--shortReads missing\n" unless($shortReads);
die "--longReads missing\n" unless($longReads);
die "--correctedReads missing\n" unless($correctedReads);

## MAIN

my $cmd;

# Generate plots of raw data.
$cmd .= " Rscript \$R_TOOLS/pacBioAssemblyPlots.R -i ".$filteredSummary." -o ".$outdir;
print STDERR "[DEBUG] ".$cmd."\n";
system($cmd);
die "Command failed: $!\n" if($? != 0);

my %hashSummary;

# Get assembly stats
my $totalContigs;    #		TotalContigsInScaffolds
my $totalBases;      #  	TotalBasesInScaffolds
my $minContigLength; #		MinContigLength
my $maxContigLength; #		MaxContigLength
my $N50Bases;        #		N50ContigBases
my $contigCoverage;  #		$contigCoverage = ContigsOnly
my $gcContent;       #		Content
my $totalSequencedBases = 0;
my $totalSequencedBasesRaw = 0;

open(OUT_TABLE_2, '>'.$outdir."/summaryTableReads2.tsv") or die "Can't open ".$outdir."/summaryTableReads2.tsv";

# Print to file relevant values from raw reads.
if( (-e $filteredSummary) and (-s $filteredSummary) ){
	open(OUT, '>'.$outdir."/summaryTableReads.tsv") or die "Can't open ".$outdir."/summaryTableReads.tsv";
	print OUT "\"Description\"\t\"Value\"\n";
	#my $readsStats = "$indir/filtering/results/filterReports_filterStats.xml";
	my $readsStats = $filteredSummary;

	my $totalReads = 0;
	my $passedReads = 0;
	my $failedReads = 0;
  my $ge3Kb = 0;
  my $ge7Kb = 0;
  my $ge9Kb = 0;
  my $ge12Kb = 0;
  my $ge15Kb = 0;
  my $ge3KbRaw = 0;
  my $ge7KbRaw = 0;
  my $ge9KbRaw = 0;
  my $ge12KbRaw = 0;
  my $ge15KbRaw = 0;
	my @readLength;
	my @passedReads;
  #my $shortestReadRaw = 0;
  #my $longestReadRaw = 0;

	open(IN, '<'.$readsStats) or die "Can't open ".$readsStats."\n";
	while(<IN>){
		chomp;
		next if($. == 1);
		my @row = split(/,/, $_);
		my $readLength = $row[3];
    next if($readLength == 0);
		#push(@row, $readLength);
		if($row[7] == 1){
      $passedReads++;
      $totalSequencedBases += $readLength;
      $ge3Kb++ if($readLength >= 3000);
      $ge7Kb++ if($readLength >= 7000);
      $ge9Kb++ if($readLength >= 9000);
      $ge12Kb++ if($readLength >= 12000);
      $ge15Kb++ if($readLength >= 15000);
			push(@passedReads, $readLength);
		}else{

    }
    $ge3KbRaw++ if($readLength >= 3000);
    $ge7KbRaw++ if($readLength >= 7000);
    $ge9KbRaw++ if($readLength >= 9000);
    $ge12KbRaw++ if($readLength >= 12000);
    $ge15KbRaw++ if($readLength >= 15000);

		$failedReads++ if($row[7] == 0);
		$totalReads++;
    $totalSequencedBasesRaw += $readLength;
		push(@readLength, $readLength);
	}
		
	my $shortestRead = min(@passedReads);
	my $longestRead = max(@passedReads);
	my $shortestReadRaw = min(@readLength);
	my $longestReadRaw = max(@readLength);
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
	print OUT "\"Number of reads greater than 3 Kb that passed QC\"\t\"".$ge3Kb."\"\n";
	print OUT "\"Number of reads greater than 7 Kb that passed QC\"\t\"".$ge7Kb."\"\n";
	print OUT "\"Number of reads greater than 9 Kb that passed QC\"\t\"".$ge9Kb."\"\n";
	print OUT "\"Number of reads greater than 12 Kb that passed QC\"\t\"".$ge12Kb."\"\n";
	print OUT "\"Number of reads greater than 15 Kb that passed QC\"\t\"".$ge15Kb."\"\n";
	close(IN);
	close(OUT);

  # Fill hash for alternate table
  $hashSummary{'1- raw polymerase reads'}{'1- totalReads'}         = $totalReads;
  $hashSummary{'1- raw polymerase reads'}{'2- sequencedBases'}     = $totalSequencedBasesRaw;
  $hashSummary{'1- raw polymerase reads'}{'3- averageReadLength'}  = $averageReadLength;
  #$hashSummary{'1- raw polymerase reads'}{'failedReads'}          = $failedReads;
  #$hashSummary{'1- raw polymerase reads'}{'passedReads'}          = $passedReads;
  $hashSummary{'1- raw polymerase reads'}{'4- shortestRead'}       = $shortestReadRaw;
  $hashSummary{'1- raw polymerase reads'}{'5- longestRead'}        = $longestReadRaw;
  $hashSummary{'1- raw polymerase reads'}{'6- >=3Kb'}              = $ge3KbRaw; 
  $hashSummary{'1- raw polymerase reads'}{'7- >=7Kb'}              = $ge7KbRaw; 
  $hashSummary{'1- raw polymerase reads'}{'8- >=9Kb'}              = $ge9KbRaw; 
  $hashSummary{'1- raw polymerase reads'}{'9- >=12Kb'}             = $ge12KbRaw;
  $hashSummary{'1- raw polymerase reads'}{'91- >=15Kb'}            = $ge15KbRaw;
 
  #$hashSummary{'2- filtered'}{'1- totalReads'}        = $passedReads;
  #$hashSummary{'2- filtered'}{'2- sequencedBases'}    = $totalSequencedBases;
  #$hashSummary{'2- filtered'}{'3- averageReadLength'} = $averageReadLengthQCpassed;
  #$hashSummary{'2- filtered'}{'4- shortestRead'}      = $shortestRead;
  #$hashSummary{'2- filtered'}{'5- longestRead'}       = $longestRead;
  #$hashSummary{'2- filtered'}{'6- >=3Kb'}             = $ge3Kb;
  #$hashSummary{'2- filtered'}{'7- >=7Kb'}             = $ge7Kb;
  #$hashSummary{'2- filtered'}{'8- >=9Kb'}             = $ge9Kb;
  #$hashSummary{'2- filtered'}{'9- >=12Kb'}            = $ge12Kb;
  #$hashSummary{'2- filtered'}{'91- >=15Kb'}           = $ge15Kb;
}

getReadsStats($shortReads, "short subreads", "3-");
getReadsStats($longReads, "long subreads", "4-");
getReadsStats($correctedReads, "corrected long subreads", "5-");

# Compute N25, N50, N75 and display length of each contig for final contigs.
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

# Compute X coverage value used by hgap.
#my $percent;
my $covCutoff;
if($suffix =~ m/(\d+)X/){
	$covCutoff = $1;
}else{
	die "Something wrong with assembly name...\n";
}
my $estimatedCov = ($totalSequencedBases / $estimatedGenomeSize);
my $estimatedCovRaw = ($totalSequencedBasesRaw / $estimatedGenomeSize);
#my $covCutoff = sprintf "%.2f", ($estimatedCov * $percent / 100);
$estimatedCov = sprintf "%.0f", $estimatedCov;

# Print to file relevant values of the assembly process.
open(OUT, '>'.$outdir."/summaryTableAssembly.tsv") or die "Can't open ".$outdir."/summaryTableAssembly.tsv";
print OUT "\"Description\"\t\"Value\"\n";
print OUT "\"Sample name\"\t\"$sampleName\"\n";
print OUT "\"Assembly name\"\t\"$sampleName-$suffix\"\n";
print OUT "\"HGAP cutoff (X)\"\t\"$covCutoff\"\n";
print OUT "\"Estimated genome size (bp)\"\t\"$estimatedGenomeSize\"\n";
print OUT "\"Estimated coverage (X)\"\t\"$estimatedCov\"\n";
print OUT "\"Total contigs\"\t\"$totalContigs\"\n";
print OUT "\"Total bases in contigs (bp)\"\t\"$totalBases\"\n";
print OUT "\"Minimum contig length (bp)\"\t\"$minContigLength\"\n";
print OUT "\"Maximum contig length (bp)\"\t\"$maxContigLength\"\n";
#print OUT "\"N25 contigs bases\"\t\"$N25\"\n";
#print OUT "\"N50 contigs bases\"\t\"$N50\"\n";
#print OUT "\"N75 contigs bases\"\t\"$N75\"\n";
#print OUT "\"N90 contigs bases\"\t\"$N90\"\n";
#print OUT "\"Contigs coverage (X)\"\t\"$contigCoverage\"\n";
print OUT "\"GC content (%)\"\t\"$gcContent\"\n";
print OUT "\"N25 - 25% of total sequence length is contained in the ".$N25count." sequence(s) having a length >= \"\t".$N25."\" bp\"\n";
print OUT "\"N50 - 50% of total sequence length is contained in the ".$N50count." sequence(s) having a length >= \"\t".$N50."\" bp\"\n";
print OUT "\"N75 - 75% of total sequence length is contained in the ".$N75count." sequence(s) having a length >= \"\t".$N75."\" bp\"\n";
print OUT "\"N90 - 90% of total sequence length is contained in the ".$N90count." sequence(s) having a length >= \"\t".$N90."\" bp\"\n";
close(OUT);

# Compute N25, N50, N75 and display length of each contig for short reads (before correction).
sub getReadsStats{
  
  my $infile = shift;
  my $prefix = shift;
  my $indice = shift;	

  open(OUT, '>>'.$outdir."/summaryTableReads.tsv") or die "Can't open ".$outdir."/summaryTableReads.tsv";

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
	my $totalReads = 0;
	my $passedReads = 0;
	my $failedReads = 0;
	my $ge3Kb = 0;
	my $ge7Kb = 0;
	my $ge9Kb = 0;
	my $ge12Kb = 0;
	my $ge15Kb = 0;
	my @readLength;
	my @passedReads;
	
	my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
	while( my $curr = $ref_fasta_db->next_seq() ) {
		my $length = length($curr->seq);
		my $header = $curr->header;
		$header =~ s/>//;
	  push(@length, $length);
	  $ge3Kb++ if($length >= 3000);
	  $ge7Kb++ if($length >= 6000);
	  $ge9Kb++ if($length >= 9000);
	  $ge12Kb++ if($length >= 12000);
	  $ge15Kb++ if($length >= 15000);
	
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
	
	my $shortestRead = min(@length);
	my $longestRead = max(@length);
	my $averageReadLength = sprintf "%.2f", sum(@length)/@length;
  #$averageReadLength = sprintf "%.2f", ($averageReadLength);
		
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
	
	print OUT "\"Total number of $prefix\"\t\"".$counter."\"\n";
	print OUT "\"Total bases in $prefix\"\t\"".$totalBases."\"\n";
	print OUT "\"Average $prefix lenth\"\t\"".$averageReadLength."\"\n";
	print OUT "\"Number of $prefix greater than 3 Kb\"\t\"".$ge3Kb."\"\n";
	print OUT "\"Number of $prefix greater than 7 Kb\"\t\"".$ge7Kb."\"\n";
	print OUT "\"Number of $prefix greater than 9 Kb\"\t\"".$ge9Kb."\"\n";
	print OUT "\"Number of $prefix greater than 12 Kb\"\t\"".$ge12Kb."\"\n";
	print OUT "\"Number of $prefix greater than 15 Kb\"\t\"".$ge15Kb."\"\n";
	print OUT "\"N25 - 25% of total $prefix sequence length is contained in the ".$N25count." sequence(s) having a length >= \"\t\"".$N25." bp\"\n";
	print OUT "\"N50 - 50% of total $prefix sequence length is contained in the ".$N50count." sequence(s) having a length >= \"\t\"".$N50." bp\"\n";
	print OUT "\"N75 - 75% of total $prefix sequence length is contained in the ".$N75count." sequence(s) having a length >= \"\t\"".$N75." bp\"\n";
	print OUT "\"N90 - 90% of total $prefix sequence length is contained in the ".$N90count." sequence(s) having a length >= \"\t\"".$N90." bp\"\n";
  close(OUT);
  
  $hashSummary{"$indice $prefix"}{'1- totalReads'}        = $counter;
  $hashSummary{"$indice $prefix"}{'2- sequencedBases'}    = $totalBases;
  $hashSummary{"$indice $prefix"}{'3- averageReadLength'} = $averageReadLength;
  $hashSummary{"$indice $prefix"}{'4- shortestRead'}      = $shortestRead;
  $hashSummary{"$indice $prefix"}{'5- longestRead'}       = $longestRead;
  $hashSummary{"$indice $prefix"}{'6- >=3Kb'}             = $ge3Kb;
  $hashSummary{"$indice $prefix"}{'7- >=7Kb'}             = $ge7Kb;
  $hashSummary{"$indice $prefix"}{'8- >=9Kb'}             = $ge9Kb;
  $hashSummary{"$indice $prefix"}{'9- >=12Kb'}            = $ge12Kb;
  $hashSummary{"$indice $prefix"}{'91- >=15Kb'}           = $ge15Kb;
}

my @string = 1;
$string[0] = "#PacBio assembly summary count table";
# Then loop through hashSummary hash and print results in a table
my $i = 1;
foreach my $prefix (sort{$a cmp $b} keys %hashSummary) {
  $string[$i] = $prefix;	

  foreach my $field (sort{$a cmp $b} keys %{ $hashSummary{$prefix} }) {
    if($i == 1){
      $string[0] .= "\t".$field;
    }
    $string[$i] .= "\t".$hashSummary{$prefix}{$field};
	}
	$i++;
}

foreach(@string){
  my $line = $_;
  $line =~ s/\d{1,2}\- //g;
  print OUT_TABLE_2 $line."\n";
}
close(OUT_TABLE_2);

exit;

### OLD CODE TO DELETE...
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


#Filtered sequence (total) total length (bp), 306232954
#Filtered sequence (total) Coverage on 2100000bp (X), 145.825
#Filtered sequence (total) average length, 3841.45
#Filtered sequence (total) N50, 4879
#Filtered sequence (total) number, 79718
#Filtered sequence (total) number greater than 3Kb, 40780
#Filtered sequence (total) number greater than 7Kb, 13893
#Long sequences (into corr) total length (bp), 63003917
#Long sequences (into corr) Coverage on 2100000bp (X), 30.001
#Long sequences (into corr) Minimum Long Readlength (bp), 8961
#Long sequences (into corr) average length, 11285
#Long sequences (into corr) N50, 11178
#Long sequences (into corr) number, 5583
#Long sequences (into corr) number greater than 3Kb, 5583
#Long sequences (into corr) number greater than 7Kb, 5583
#Short sequences (into corr) total length (bp), 243229037
#Short sequences (into corr) Coverage on 2100000bp (X), 115.823
#Short sequences (into corr) average length, 3280.89
#Short sequences (into corr) N50, 4070
#Short sequences (into corr) number, 74135
#Short sequences (into corr) number greater than 3Kb, 35197
#Short sequences (into corr) number greater than 7Kb, 8310
#Corrected sequences (total) total length (bp), 46392459
#Corrected sequences (total) Coverage on 2100000bp (X), 22.091
#Corrected sequences (total) average length, 7919.5
#Corrected sequences (total) N50, 9844
#Corrected sequences (total) number, 5858
#Corrected sequences (total) number greater than 3Kb, 4993
#Corrected sequences (total) number greater than 7Kb, 4118
#Corrected sequence (intoasm) total length (bp), 46392459
#Corrected sequence (intoasm) Coverage on 2100000bp (X), 22.091
#Corrected sequence (intoasm) Minimum Readlength (bp), 500
#Corrected sequence (intoasm) average length, 7919.5
#Corrected sequence (intoasm) N50, 9844
#Corrected sequence (intoasm)s, 5858
#Corrected sequence (intoasm)s greater than 3Kb, 4993
#Corrected sequence (intoasm)s greater than 7Kb, 4118

