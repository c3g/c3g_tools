#!/usr/bin/env perl

use strict;
use warnings;

use Env qw(SCRATCH TMPDIR PICARD_HOME);
use Getopt::Long;
use Data::Dumper;
use File::Find;

my $usage=<<'ENDHERE';
NAME:
bamToFastq.pl

PURPOSE:
To loop recursively through a directory and convert each .bam
to fastq.gz file.

INPUT:
--indir <string>  : Directory where itags runs are to be backed-up.

OUTPUT:
No output options.

NOTES:
module load mugqic/java/oracle-jdk1.7.0_15
module load mugqic/picard/1.108 

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $indir);
my $verbose = 0;

GetOptions(
    'indir=s'   => \$indir,
    'verbose'   => \$verbose,
    'help'      => \$help
);
if ($help) { print $usage; exit; }

## Validate
die "--indir arg missing\n" unless($indir);

## It is intended to recursively through a directory tree
## and compress all *.fastq in *.fastq.gz


## MAIN
sub eachFile{
	my $filename = $_;
	my $fullpath = $File::Find::name;
	#remember that File::Find changes your CWD, 
	#so you can call open with just $_

	if (-e $filename) { 
		
		if(substr($filename, -4) eq ".bam"){
			print STDERR "Decompressing ".$fullpath." into .gz archive...\n";
      print STDERR "Don't forget to : module load mugqic/java/oracle-jdk1.7.0_15 mugqic/picard/1.108\n";
     
      my $pair1 = $fullpath;
      my $pair2 = $fullpath;
      $pair1 =~ s/\.bam/\.pair1\.fastq/;
      $pair2 =~ s/\.bam/\.pair2\.fastq/;
  
      my $cmd = "";
      $cmd .= " java -Djava.io.tmpdir=$SCRATCH -XX:ParallelGCThreads=3 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar";
      $cmd .= " INPUT=".$fullpath;
      $cmd .= " FASTQ=".$pair1;
      $cmd .= " SECOND_END_FASTQ=".$pair2;
      $cmd .= " && gzip ".$pair1." && gzip ".$pair2;
      print STDERR "Converting $fullpath bam to fastq\n";
      print STDERR "============\n".$cmd."\n=============\n";
			system($cmd);
			$? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly decompressed ".$fullpath." into .fastq archive...\n";
		}
	}
}

## MAIN

# Compress .fastq into .gz
find (\&eachFile, $indir);

exit;
