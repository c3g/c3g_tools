#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Find;

my $usage=<<'ENDHERE';
NAME:
sampleSetupPacBio.pl

PURPOSE:
The script will take in input a run name and will loop through
that run directory and create symlinks in raw_reads/
Usually in a PAcBio folder, there are ~8 wells.

INPUT:
--runName <string> : Example: Run044_379 
--infile <string>  : run Id csv file from nanuq.
				
OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $runName, $nanuqAuthFile);
my $verbose = 0;

GetOptions(
    'runName=s'	      => \$runName,
	#'infile=s'       => \$infile,
	'nanuqAuthFile=s' => \$nanuqAuthFile,
    'verbose' 	      => \$verbose,
    'help' 		      => \$help
);
if ($help) { print $usage; exit; }

# Get sample sheet from Run name.
sub getSheet {
	#my $tech = shift;
	my $projectFile = shift;
	my $runName = shift;
	my $nanuqAuthFile = shift;

	# Delete previous project.nanuq.csv
	if(-e "./$projectFile") {
		system("rm ./$projectFile");
	}
 
	# RUN NAME 
	#my $command = 'wget --no-cookies --post-file '.$nanuqAuthFile.' https://genomequebec.mcgill.ca/nanuqMPS/csv/technology/Pacbio/run/'.$runName.'/filename/'.$projectFile."\n";
	
	# OR BY PROJECT NAME
	my $command = 'wget --no-cookies --post-file '.$nanuqAuthFile.' https://genomequebec.mcgill.ca/nanuqMPS/csv/technology/Pacbio/project/'.$runName.'/filename/'.$projectFile."\n";
	print STDERR '#'.$command;
	system($command);
	if ($? == -1) {
		print STDERR "failed to execute: $!\n";
		exit(1);
	}
	elsif ($? & 127) {
		printf STDERR "child died with signal %d, %s coredump\n", ($? & 127), ($? & 128) ? 'with' : 'without';
		exit(1);
	}
	else {
		my $childValue = $? >> 8;
		if($childValue != 0) {
			printf STDERR "child exited with value %d\n", $childValue;
			exit(1);
		}
	}
}

## MAIN
die "--runName arg missing.\n" unless($runName);
#die "--infile arg missing.\n" unless($infile);
die "--nanuqAuthFile arg missing.\n" unless($nanuqAuthFile);

#my $projectFile = "run.nanuq.csv";
my $projectFile = "project.nanuq.csv";
# Get sample Sheet from RunId (and not projectId)
getSheet($projectFile, $runName, $nanuqAuthFile);
#exit;
#my $root = "/lb/robot/pacbioSequencer/pacbioRuns/$runName";
my $root;
my $infile = "./".$projectFile;

#/lb/robot/pacbioSequencer/pacbioRuns/Run044_379/A04_1/Analysis_Results#

# First parse csv file to find well name associated with what sample.
my %hash;
my @prefixes;
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
	chomp;
	next if($. == 1);
	my @row = split(/","/, $_);
	my $sampleName = $row[0];
	my $protocol = $row[4];
	my $well = $row[11];
	my $run = $row[10];
	my $bp = $row[22];

	$sampleName =~ s/_MPS.*//;
	$well =~ s/"//g;
	$sampleName =~ s/"//g;
	$protocol =~ s/"//g;
	$run =~ s/"//g;
	$bp =~ s/\(|\)//g;
	my @bp = split(/\//, $bp);
	$bp = $bp[2];
	$bp =~ s/,//g;

	$hash{$run}{$well}{run} = $run;
	$hash{$run}{$well}{sampleName} = $sampleName;
	$hash{$run}{$well}{protocol} = $protocol;
	$hash{$run}{$well}{well} = $well;
	$hash{$run}{$well}{bp} = $bp;
}

close(IN);

foreach my $run (%hash){
	foreach my $well (keys %{ $hash{$run} }) {

		# Hack to add _\d+ at the end of $root
		my $dirname = "/lb/robot/pacbioSequencer/pacbioRuns";
		opendir(DIR, $dirname);
		my @files = readdir(DIR);
		closedir DIR;

		foreach my $key (@files){
			if(-d "$dirname/$key"){
				my $searchString = $run."(_\\d+)";
				if($key =~ m/$searchString/){
					$root = "$dirname/$run$1/$well/Analysis_Results/";

					opendir(DIR, $root);
					@files = readdir(DIR);
					closedir DIR;
					foreach $key (@files){
						my $searchString = "(.*.bas.h5)";
						if($key =~ m/$searchString/){
							$root .= "$1";
							$hash{$run}{$well}{path} = $root;
							
							print STDOUT $hash{$run}{$well}{sampleName}."\t".$well."\t".$hash{$run}{$well}{path}."\t".$hash{$run}{$well}{protocol}."\t".$hash{$run}{$well}{bp}."\n"; 	
						}
					}
				}
			}
		}
	}
}
	
sub eachFile{
	my $filename = $_;

	#remember that File::Find changes your CWD, 
	#so you can call open with just $_
	my $fullpath = $File::Find::name;
	
	my @path = split(/\//, $fullpath);
	my $run = pop(@path);
	$run =~ m/(Run\d+)_\d+/;
	$run = $1;
	print "RUN:".$run."\n";
	
	if($filename =~ m/qc/ && ( -d $filename)){
		# Just want to avoid .../qc/... dir.	
	}elsif( substr($filename, -7) eq ".bas.h5" ){
		
		if($fullpath =~ m/\/((\w\d{2,3})_\d)\//){
			print STDERR "run:$run\n\$1:$1\n";
			print STDOUT $hash{$run}{$1}{sampleName}."\t".$1."\t".$fullpath."\t".$hash{$run}{$1}{protocol}."\n"; 	
			#$hash{$1} = $fullpath; 
			exit;
		}
	}
}

exit;



