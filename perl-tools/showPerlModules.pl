#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use ExtUtils::Installed;

my $usage=<<'ENDHERE';
NAME:
showPerlModules.pl

PURPOSE:

INPUT:
				
OUTPUT:
STDOUT : print a list of perl modules and their version.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $infile);
my $verbose = 0;

GetOptions(
    'infile=s' 	=> \$infile,
    'verbose' 	=> \$verbose,
    'help' 		=> \$help
);
if ($help) { print $usage; exit; }

my $instmod = ExtUtils::Installed->new();
foreach my $module($instmod->modules()){
	my $version = $instmod->version($module) || "???";
	print STDOUT "$module -- $version\n";
}
exit;
