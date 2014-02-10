#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use List::Util qw(sum);
use File::Slurp;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
pacBioGetCutoff.pl

PURPOSE:
Loop through fasta sequences. put the length of each
sequence in an array, sort it by descinding order, 
loop through it again and compute the cummulative length 
covered by each sequence as we loop though the array. 
Once that length is > (coverage * genome size) * covFraction, 
we have our threshold. The idea is to consider all reads 
above that threshold to be seeding reads to which will be 
align lower shorter subreads.

INPUT:
--infile <string>          : Sequence file (either fasta or fastq)
--coverage <int>           : Coverage
--genomeSize <int>         : Estimated genome size
--coverageFraction <float> : At what cutoff of the estimated coverage you 
                             want the cutoff value. Default: 0.30

--coverageCutoff           : To get cutoff based on cummulative length
                             of sequences up until a (Coverage * 
                             Genome size) * coverageFraction.
--xml <string>             : XML file to parse.
--xmlOut <string>          : XML outfile with modified value in the 
                             minReadLength field.

TODO
--hgapCutoff               : To get cutoff based on a more 
                             sophisticated algorithm.
			
OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $infile, $coverage, $genomeSize, $coverageFraction, $coverageCutoff, $hgapCutoff, $xml, $xmlOut);
my $verbose = 0;

GetOptions(
    'infile=s' 	   		 => \$infile,
	'xml=s'				 => \$xml,
	'xmlOut=s'			 => \$xmlOut,
	'genomeSize=i' 		 => \$genomeSize,
	'coverageCutoff' 	 => \$coverageCutoff,
	'coverageFraction=f' => \$coverageFraction,
	'hgapCutoff'		 => \$hgapCutoff,
	'coverage=i'   		 => \$coverage,
    'verbose' 	   		 => \$verbose,
    'help' 		   		 => \$help
);
if ($help) { print $usage; exit; }

## Validate
die "--infile missing\n" unless($infile);
die "--genomeSize missing\n" unless($genomeSize);
die "--coverageCutoff or --hgapCutoff missing\n" unless($coverageCutoff || $hgapCutoff);
if($xml){
	#die "--xmlOut missing\n" unless($xmlOut);
}
if($xmlOut){
	#die "--xml missing\n" unless($xml);
}

$coverageFraction = 0.30 unless($coverageFraction);

## MAIN

# You should estimate the overall coverage and length distribution for putting in
# the correct options in the configuration file.  You will need to decide a
# length cutoff for the seeding reads. The optimum cutoff length will depend on
# the distribution of the sequencing read lengths, the genome size and the
# overall yield. The general guideline is the coverage of the seeding sequences
# should be above 20x of the genome and the overall coverage should be at least
# 3x of the coverage of the seeding sequences.
# ===============================================================================
# Here loop through fasta sequences. put the length of each
# sequence in an array, sort it, loop through it again and 
# compute the cummulative length covered
# by each sequence as we loop though the array.
# Once that length is > (coverage * genome size), we have
# our threshold. The idea is to consider all reads above
# that threshold to be seeding reads to which will be 
# align lower shorter subreads.

my $readCount = 0;
my %hash;
my @array;
my $totalBases=0;

if($infile =~ m/\.fastq|\.fq/){
	my $ref_fastq_db = Iterator::FastqDb->new($infile) or die("Unable to open Fasta file, $infile\n");
	while( my $curr = $ref_fastq_db->next_seq() ) {
		$hash{length($curr->seq())}++;
		push(@array, length($curr->seq()));
	}
}elsif($infile =~ m/\.fasta|\.fa|\.fsa|\.fna/){
	my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
	while( my $curr = $ref_fasta_db->next_seq() ) {
		$hash{length($curr->seq())}++;
		$totalBases += length($curr->seq());
		push(@array, length($curr->seq()));
	}
}else{
	die "Couldn't determine if file was fastq or fasta.\n";
}

$coverage = sprintf "%0.2f", $totalBases / $genomeSize;

my $cutoff=0;
if($coverageCutoff){
	@array = sort {$b <=> $a} @array; # Descending order. First reads in the array (before to be found cutoff) will be seeding reads.
	my $sum=0;
	for (my $i=0;$i<@array;$i++){
		$sum = $sum + $array[$i];
		if($sum > ( ($coverage * $genomeSize) * $coverageFraction )){ # Coult be refined... this cutoff depends on the read length distribution.
			$cutoff = $array[$i]; 
			last;
		}
	}

	my $cov = sprintf "%0.2f", $totalBases / $genomeSize;
	my $X = sprintf "%0.2f", $sum / $genomeSize;

	print STDERR "Estimated coverage with filtered reads only: ".$cov."X\n";
	print STDERR "Estmated coverage cutoff value: ".$X."X\n";
	print STDERR "Read length cutoff value: ".$cutoff." bp\n";	

	print STDOUT $cutoff."\n";

}elsif($hgapCutoff){

	@array = sort {$b <=> $a} @array; # Descending order. First reads in the array (before to be found cutoff) will be seeding reads.
	my $total = sum(@array);
	my @coverageArray;

	my @x = (500,100,6000);
	foreach my $x (@x){
		my @y;
		foreach my $seqLength (@array){
			push(@y, $seqLength) if($seqLength > $x);
		}
		my $psum = sum(@y);
		
		my $coverage = 0.5 * $psum / $genomeSize / $x * exp( $coverage *-1 );
		my $contigCount = $coverage * $genomeSize / $x * exp ($coverage * -1);
		my $contigLength = (exp($coverage) -1) * $x / $coverage;

		my @el = ($x, $psum, $total / $psum, $coverage, $contigCount, $contigLength / $genomeSize);
		push @coverageArray, [ @el ];
	}
		
	#seq_lengths = sort(seq_lengths,'descend');
	#% coverage = sum(seq_lengths(seq_lengths>=lencut))/estimatedgenomesize;
	#
	#total = sum(seq_lengths)
	#coverage_array = []
	#
	#for x=500:100:6000,
	#	psum = sum(seq_lengths(seq_lengths>x));
	#    coverage = 0.5 * psum / estimatedgenomesize; % we loss 50% bases after the pre-assembly step
	#    contig_count = coverage * estimatedgenomesize / x * exp( -coverage );
	#    contig_length = (exp(coverage) - 1) * x /coverage;
	#	el=[ x, psum, 1.0*total/psum, coverage,  contig_count,  contig_length/estimatedgenomesize];
	#    disp(sprintf('%i, %i, %f, %0.2f, %f, %f', x, psum, 1.0*total/psum, coverage,  contig_count,  contig_length/estimatedgenomesize))
	#	coverage_array = [coverage_array; el];
	#end
	#
	#% 1. The ratio of the total number bases to the long read bases is larger than 3.
	#% 2. Estimated Lander-Waterman contig number less than 0.25. 
	#% 3. The estimated Lander-Waterman contig size is larger than 0.25x of the genome size.
	#
	#answer = coverage_array(coverage_array(:,3)>3 & coverage_array(:,5)<0.25 & coverage_array(:,6)>0.25,:);
	#if(~isempty(answer))
	#	cutoff = answer(1);
	#	disp([ 'Longest start at,' num2str(cutoff) ]);
	#else
	#	% [v,i] = min(sum([abs(coverage_array(:,3)-3)/mean(coverage_array(:,3)) abs(coverage_array(:,5)-0.25)/mean(coverage_array(:,5)) abs(coverage_array(:,6)-0.25)/mean(coverage_array(:,6))],2));
	#	[v,i] = min(abs(coverage_array(:,3)-3));
	#	%
	#	% coverage_array(:,3) - 3 >0
	#	% 0.25 - coverage_array(:,5) >0
	#	% coverage_array(:,6) - 0.25x >0
	#	% 
	#	cutoff = coverage_array(i,1);
	#	disp([ 'Longest start at,' num2str(cutoff) ]);
	#end
}

if($xml && $xmlOut){
	my $text = read_file( $xml );
	$text =~ s/(<param name=\"minLongReadLength\">\n.*<value>)(.*)(<\/value>\n.*<\/param>)/$1$cutoff$3/;
	$text =~ s/(<title>Compute \"Minimum Seed Read Length\"<\/title>\n.*<value>)True(<\/value>\n.*<input type=\"checkbox\" \/>)/$1False$2/;
	open(XMLOUT, '>'.$xmlOut) or die "Can't open $xmlOut\n";
	print XMLOUT $text;
	close(XMLOUT);
}

exit;
