#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Find;

my $usage=<<'ENDHERE';
NAME:
pacBioGenerateFofns.pl

PURPOSE:
Loop in a directory, find every .h5 file and
write the path in a file.

INPUT:
--indir <string> : Indir containing (.h5)
				
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
    'indir=s' 	=> \$indir,
    'verbose' 	=> \$verbose,
    'help' 		=> \$help
);
if ($help) { print $usage; exit; }

## MAIN
# Loop through directory
find (\&eachFile, $indir);

sub eachFile{
	my $filename = $_;
	my $fullpath = $File::Find::name;
	#remember that File::Find changes your CWD, 
	#so you can call open with just $_
	
	if($filename =~ m/qc/ && ( -d $filename)){
		# Just want to avoid .../qc/... dir.	
	}elsif( substr($filename, -3) eq ".h5" ){
		
		#if($fullpath =~ m/\/((\w\d{2,3})_\d)\//){
			print STDOUT $fullpath."\n"; 	
			#$hash{$1} = $fullpath; 
		#}
	}
}

exit;

