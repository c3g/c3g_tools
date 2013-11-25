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
--runName <string>       : Example: Run044_379 
Or 
--projectId <int>        : ex: 9999
Or
--sampleSheet <string>   : Path to your sample sheet already downloaded from nanuq.

--nanutAuthFile <string> : File having username and pass on one line.
	
OUTPUT:
Sample sheet needed to launch the PacBioi assembly pipeline.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $runName, $projectId, $sampleSheet, $nanuqAuthFile);
my $verbose = 0;

GetOptions(
    'runName=s'	      => \$runName,
	'projectId=i'     => \$projectId,
	'sampleSheet=s'   => \$sampleSheet,
	'nanuqAuthFile=s' => \$nanuqAuthFile,
    'verbose' 	      => \$verbose,
    'help' 		      => \$help
);
if ($help) { print $usage; exit; }

# Get sample sheet from Run name.
sub getSheet {
	my $projectFile   = shift;
	my $nanuqAuthFile = shift;
	my $runName       = shift;
	my $projectId     = shift;

	# Delete previous project.nanuq.csv
	if(-e "./$projectFile") {
		system("rm ./$projectFile");
	}
 
	# Run name or project ID.
	my $command;
	if(defined($runName)){
		$command = 'wget --no-cookies --post-file '.$nanuqAuthFile.' https://genomequebec.mcgill.ca/nanuqMPS/csv/technology/Pacbio/run/'.$runName.'/filename/'.$projectFile."\n";
	}
	if(defined $projectId){
		$command = 'wget --no-cookies --post-file '.$nanuqAuthFile.' https://genomequebec.mcgill.ca/nanuqMPS/csv/technology/Pacbio/project/'.$projectId.'/filename/'.$projectFile."\n";
	}

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
#unless( $runName || $projectId ){
#	die "--runName or --projectId arg missing.\n";
#}
die "--nanuqAuthFile arg missing.\n" unless($nanuqAuthFile);

my $projectFile = "project.nanuq.csv";
my $root;
my $infile;

# If no sampleSheet infile was provided.
if(!$sampleSheet){
	# Get sample Sheet from RunId (and not projectId)
	if($runName){
		getSheet($projectFile, $nanuqAuthFile, $runName, $projectId);
	}elsif($projectId){
		getSheet($projectFile, $nanuqAuthFile, $runName, $projectId);
	}
	$infile = "./".$projectFile;

}else{ # Else take provided sampleSheet.
	$infile = $sampleSheet;
}

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
	my $well = $row[17];
	my $run = $row[16];
	my $bp = $row[28];

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



