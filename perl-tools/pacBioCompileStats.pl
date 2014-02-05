#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use List::Util qw( sum min max);
use File::Slurp;
use Iterator::FastaDb;
use File::Basename;
use File::Find;
use POSIX;

my $usage=<<'ENDHERE';
NAME:
pacBioAssemblyStats.pl

PURPOSE:
Generate relevant plots and table(s) of current PacBio assembly.

INPUT:
--indir <string>            : Directory containing sample's assemblies.
--estimatedGenomeSize <int> : Estimated genome size.
--sampleName <string>       : Sample name.

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $indir, $sampleName, $estimatedGenomeSize);
my $verbose = 0;

GetOptions(
	'indir=s'                 => \$indir,
	'sampleName=s'            => \$sampleName,
	'estimatedGenomeSize=i'   => \$estimatedGenomeSize,
    'verbose' 	     		  => \$verbose,
    'help' 		     		  => \$help
);
if ($help) { print $usage; exit; }

## Validate
die "--indir missing\n" unless($indir);

## Global scope variable
my %hQc;
my $readsStats;
my %hQcPolished;
my %hFasta;
$sampleName = "unknown" unless($sampleName);
print STDERR "[DEBUG] sampleName: $sampleName\n" if($verbose);

## SUBS
sub eachFile{
	my $filename = $_;
	my $fullpath = $File::Find::name;
	#remember that File::Find changes your CWD, 
	#so you can call open with just $_

	if (-e $filename) { 
		if($fullpath =~ m/\d+percent\/merSize\d+\/polishing\d+\/data\/consensus\.fasta$/){
			if($fullpath =~ m/(\d+percent).*(merSize\d+).*(polishing\d+)/){
				print STDERR "[DEBUG] $fullpath\t$1\t$2\t$3\n" if($verbose);
				$hQcPolished{"$1_$2_$3"} = $fullpath;
			}else{
				die "Can't find proper nomenclature on .fasta files.\n";	
			}
		
		}elsif($fullpath =~ m/9-terminator\/$sampleName\S+\.qc$/){
			if($fullpath =~ m/(\d+percent)_(merSize\d+)/){
				print STDERR "[DEBUG] $fullpath\n" if($verbose);
				$hQc{"$1_$2"} = $fullpath;
			}else{
				die "Can't find proper nomenclature on .qc files.\n";	
			}
		# Only supposed to be one file in data structure.
		}elsif($fullpath =~ m/filtering\/data\/filtered_summary.csv/){
			print STDERR "[DEBUG] $fullpath\n" if($verbose);
			$readsStats = $fullpath;
		}
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

my $totalSequencedBases = 0;
open(IN, '<'.$readsStats) or die "Can't open ".$readsStats."\n";
while(<IN>){
	chomp;
	next if($. == 1);
	my @row = split(/,/, $_);
	my $readLength = $row[3];
	push(@row, $readLength);
	if($row[7] == 1){
        $totalSequencedBases += $readLength;					
	}
}
close(IN);
########################################################
# Compute stats directly from fasta
#
my %hash;
my %hResultsPolished;
foreach my $key (sort{$a cmp $b} keys %hQcPolished){
	my $contigs = $hQcPolished{$key};
	my $totalBases = 0;
	my $counter = 1;
	my @length;
	my $gcCount = 0;

	my @seqLengths;
	my $opt_i = 100;
	my %len= ();
	my $n = 0;
	my $int;
	my $totalLength = 0;
	my $covCutoff;

	if($estimatedGenomeSize){	
		# Compute X coverage value used by hgap.
		my $percent;
		if($key =~ m/(\d+)percent/){
			$percent = $1;
		}else{
			die "Something wrong with assembly name...\n";
		}
		my $estimatedCov = ($totalSequencedBases / $estimatedGenomeSize);
		$covCutoff = sprintf "%.2f", ($estimatedCov * $percent / 100);
	}

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
		###################
	
		$totalBases += $length;
		$hash{$header} = $length;
		$gcCount += ($curr->seq()  =~ tr/gGcC/gGcC/);
		$counter++;
	}
	
	my $gcContent = sprintf "%.2f", ($gcCount/$totalBases * 100);
	
	my $maxContigLength = max @length;
	my $minContigLength = min @length;
	
	$counter = ($counter - 1);
	$hResultsPolished{$key}{'totalContigsInScaffolds'} = $counter;
	$hResultsPolished{$key}{'totalBasesInScaffolds'}   = commify($totalBases);
	$hResultsPolished{$key}{'minContigLength'}         = commify($minContigLength);
	$hResultsPolished{$key}{'maxContigLength'}         = commify($maxContigLength);
	$hResultsPolished{$key}{'HGAP_cutoff(X)'}          = commify($covCutoff) if($estimatedGenomeSize);
	
	print STDOUT "Total of ".$counter." sequences\n" if($verbose);
		
	#foreach my $key4 (sort {$hash{$b} <=> $hash{$a}} (keys %hash)){
	#	print STDOUT $key4 . "\t" .$hash{$key4}."\n";	
	#}
	
	my $cummSum = 0;
	$counter = 1;
	
	print STDERR "[DEBUG] ".$totalLength."\n" if($verbose);
	
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
 
	$hResultsPolished{$key}{'N25'}                     = commify($N25);
	$hResultsPolished{$key}{'N50'}                     = commify($N50);
	$hResultsPolished{$key}{'N75'}                     = commify($N75);
	#$hResultsPolished{$key}{'contigs coverage'}        = $contigCoverage;
	$hResultsPolished{$key}{'gc%'}                     = $gcContent;	
}

print STDERR "[DEBUG] Done with fasta files.\n" if($verbose);

###########################################################
# Compute N25, N50, N75 and display length of each contig.
# from celera contigs (unpolished)
my %hResults;

foreach my $key (sort{$a cmp $b} keys %hQc){
	print STDERR "[DEBUG] ".$key."\t".$hQc{$key}."\n" if($verbose);
	my $file = $hQc{$key};

	# Get assembly stats
	my $totalContigsCelera;    #		TotalContigsInScaffolds
	my $totalBasesCelera;      #     	TotalBasesInScaffolds
	my $minContigLengthCelera; #		MinContigLength
	my $maxContigLengthCelera; #		MaxContigLength
	my $N50BasesCelera;        #		N50ContigBases
	my $N25BasesCelera;        #		N25ContigBases
	my $N75BasesCelera;        #		N75ContigBases
	my $contigCoverageCelera;  #		$contigCoverageCelera = ContigsOnly
	my $gcContentCelera;       #		Content
	
	open(IN, '<'.$file) or die "Can't open $file\n";
	while(<IN>){
		chomp;
		if($_ =~ m/TotalContigsInScaffolds=(.*)$/){
			$totalContigsCelera = $1;	
		}
		if($_ =~ m/TotalBasesInScaffolds=(.*)$/){
			 $totalBasesCelera = commify($1);
		}
		if($_ =~ m/MinContigLength=(.*)$/){
			$minContigLengthCelera = commify($1);
		}
		if($_ =~ m/MaxContigLength=(.*)$/){
			$maxContigLengthCelera = commify($1);
		}
		if($_ =~ m/N50ContigBases=(.*)$/){
			$N50BasesCelera = commify($1);
		}
		if($_ =~ m/N25ContigBases=(.*)$/){
			$N25BasesCelera = commify($1);
		}
		if($_ =~ m/N75ContigBases=(.*)$/){
			$N75BasesCelera = commify($1);
		}
		if($_ =~ m/ContigsOnly=(.*)$/){
			$contigCoverageCelera = $1."X";
		}
		if($_ =~ m/Content=(.*)$/){
			$gcContentCelera = $1;
		}
	}
	close(IN);

	$hResults{$key}{'totalContigsInScaffolds'} = $totalContigsCelera;
	$hResults{$key}{'totalBasesInScaffolds'}   = $totalBasesCelera;
	$hResults{$key}{'minContigLength'}         = $minContigLengthCelera;
	$hResults{$key}{'maxContigLength'}         = $maxContigLengthCelera;
	$hResults{$key}{'N25'}                     = $N25BasesCelera;
	$hResults{$key}{'N50'}                     = $N50BasesCelera;
	$hResults{$key}{'N75'}                     = $N75BasesCelera;
	$hResults{$key}{'contigs coverage'}        = $contigCoverageCelera;
	$hResults{$key}{'gc%'}                     = $gcContentCelera;	
}	

print STDOUT "Stats for sample: ".$sampleName."\n";

my @string;
$string[0] = "#Unpolished assembly name";
# Then loop through result hash and print results in a table
my $i = 1;
foreach my $keyResults (sort{$a cmp $b} keys %hResults) {
	$string[$i] = $keyResults;	

    foreach my $keyResults2 (sort{$a cmp $b} keys %{ $hResults{$keyResults} }) {
		if($i == 1){
			$string[0] .= "\t".$keyResults2;
		}
		$string[$i] .= "\t".$hResults{$keyResults}{$keyResults2}
	}
	$i++;
}

foreach my $line (@string){
	print STDOUT "$line\n";
}
print STDOUT "\n";

@string = 1;
$string[0] = "#Polished Assembly name";
# Then loop through result hash and print results in a table
$i = 1;
foreach my $keyResults1 (sort{$a cmp $b} keys %hResultsPolished) {
	$string[$i] = $keyResults1;	

    foreach my $keyResults2 (sort{$a cmp $b} keys %{ $hResultsPolished{$keyResults1} }) {
		if($i == 1){
			$string[0] .= "\t".$keyResults2;
		}
		$string[$i] .= "\t".$hResultsPolished{$keyResults1}{$keyResults2}
	}
	$i++;
}

foreach my $line (@string){
	print STDOUT "$line\n";
}

exit;

