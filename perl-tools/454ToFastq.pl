#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use File::Find;
use File::Basename;
use File::Spec;
binmode STDOUT, ":utf8";

my $usage=<<'ENDHERE';
NAME:
scriptName.pl

PURPOSE:

INPUT:
--cutAt <int>       : Final length of output file. Usually, reads quality
                      for 454 data drops after position 200 to 350. Bases
                      after these positions should be discarder because of
                      their low quality.
Either both:
--fasta <string>    : Sequence file
--qual <string>     : Qual file
OR 
--indir <string>    : absolute path to fasta and qual files.


OUTPUT:
--outIndex <string> : Indexes (barcodes)        
--outFasta <string> : Indexes in fasta format (barcodes)        
STDOUT (Fastq reads).

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, @fasta, @qual, $outIndex, $outFasta, $indir, $cutAt);
my $verbose = 0;

GetOptions(
  'indir=s'    => \$indir,
  'fasta=s'    => \@fasta,
  'qual=s'     => \@qual,
  'outIndex=s' => \$outIndex,
  'outFasta=s' => \$outFasta,
  'cutAt=i'    => \$cutAt,
  'verbose'    => \$verbose,
  'help'       => \$help
);
if ($help) { print $usage; exit; }

my %hashFiles;

# SUB
sub eachFile{
  my $filename = $_;
  my $fullpath = $File::Find::name;
  #remember that File::Find changes your CWD, 
  #so you can call open with just $_

  if (-e $filename) { 
    
    if(substr($filename, -5) eq ".qual"){
      my $base = basename($fullpath);
      $base =~ s{.*/}{};
      $base =~ s{\.[^.]+$}{};
      $hashFiles{$base}{qual} = $fullpath;
    }
    if(substr($filename, -4) eq ".fna"){
      my $base = basename($fullpath);
      $base =~ s{.*/}{};
      $base =~ s{\.[^.]+$}{};
      $hashFiles{$base}{fasta} = $fullpath;
    }
  }
}

## MAIN
die("--outIndex missing\n") unless($outIndex);
die("--outFasta missing\n") unless($outFasta);
open(OUTINDEX, ">".$outIndex) or die "Can't open $outIndex\n";
open(OUTFASTA, ">".$outFasta) or die "Can't open $outFasta\n";
#die("--fasta and --qual must be of equal length\n") if(@fasta != @qual);

# Get barcodes
my $i=0;
my $j=0;
my @barcodes;
for(glob '{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}'){
  if($i % 30 == 0){ 
    $barcodes[$j] = $_; 
    $j++;
  }else{

  }   
  $i++;
} 

my $counter=0;

if($indir){
  
  # Make sure that indir is absolute path.
  $indir = File::Spec->rel2abs($indir);

  find (\&eachFile, $indir);
  for my $key (keys %hashFiles){
    my $fasta = $hashFiles{$key}{fasta};
    my $qual = $hashFiles{$key}{qual};

    my %hash; 
  
    print STDERR "[DEBUG] processing $fasta\n";
    print STDERR "[DEBUG] processing $qual\n";

    #open(FASTA, "<".$fasta) or die "Can't open $fasta\n";
    my $ref_fasta_db = Iterator::FastaDb->new($fasta) or die("Unable to open Fasta file, $fasta\n");
    while( my $curr = $ref_fasta_db->next_seq() ) {
      my $header = $curr->header;
      $header =~ s/>//;
      $hash{$header}{fasta} = $curr->seq;
    }

    open(QUAL, "<".$qual) or die "Can't open $qual\n";
    my $tmp_title;
    while(<QUAL>){
      chomp;
      if ($_ =~ /^>/) {
        $tmp_title = $_;
        $tmp_title =~ s/>//;
        #print STDERR $tmp_title."\n";
      } else {
        $hash{$tmp_title}{qual} .= $_." ";
        #print STDERR "$_\n";
      }          
    }
    close(QUAL);
  
    # THEN print as fastq
    for my $key (sort{$a cmp $b} keys %hash){
      my $asciiQual = "";
      for( split / /, $hash{$key}{qual} ){
        $asciiQual .= chr($_ + 33);
      }
      my $header = $key;
      $header =~ s/ /_/g;
      $header .= "#".$barcodes[$counter]."/1";
    
      if($cutAt){
        if(length($hash{$key}{fasta}) >= $cutAt){
          print STDOUT "@".$header."\n".substr($hash{$key}{fasta}, 0, $cutAt)."\n+\n".substr($asciiQual, 0, $cutAt)."\n";
        }else{
          # Do nothin', discard read.
        }
      }else{
        print STDOUT "@".$header."\n".$hash{$key}{fasta}."\n+\n".$asciiQual."\n";
      }
  
    }
  
    print OUTINDEX $fasta."\t".$barcodes[$counter]."\n";
    print OUTFASTA $fasta."\n".$barcodes[$counter]."\n";
  
    $counter++;
  
  }
  exit;
}

foreach my $fasta (@fasta){
  my %hash;
  my $qual = shift(@qual);

  print STDERR "[DEBUG] processing $fasta\n";
  print STDERR "[DEBUG] processing $qual\n";

  #open(FASTA, "<".$fasta) or die "Can't open $fasta\n";
  my $ref_fasta_db = Iterator::FastaDb->new($fasta) or die("Unable to open Fasta file, $fasta\n");
  while( my $curr = $ref_fasta_db->next_seq() ) {
    my $header = $curr->header;
    $header =~ s/>//;
    $hash{$header}{fasta} = $curr->seq;
  }

  open(QUAL, "<".$qual) or die "Can't open $qual\n";
  my $tmp_title;
  while(<QUAL>){
    chomp;
    die "Length = 0 \n" if(length($_) <= 2);
    print STDERR length($_)."\n";
    if ($_ =~ /^>/) {
      $tmp_title = $_;
      $tmp_title =~ s/>//;
      #print STDERR $tmp_title."\n";
    } else {
      $hash{$tmp_title}{qual} .= $_." ";
      #print STDERR "$_\n";
    }          
  }
  close(QUAL);

  # THEN print as fastq
  for my $key (sort{$a cmp $b} keys %hash){
    my $asciiQual = "";
    for( split / /, $hash{$key}{qual} ){
      $asciiQual .= chr($_ + 33);
      print STDERR $key."\n";
    }
    my $header = $key;
    $header =~ s/ /_/g;
    $header .= "#".$barcodes[$counter]."/1";
   
    print STDERR $cutAt."\n"; 
    if($cutAt){
      if(length($hash{$key}{fasta}) >= $cutAt){
        print STDOUT "@".$header."\n".substr($hash{$key}{fasta}, 0, $cutAt)."\n+\n".substr($asciiQual, 0, $cutAt)."\n";
      }else{
        # Do nothin', discard read.
      }
    }else{
      print STDOUT "@".$header."\n".$hash{$key}{fasta}."\n+\n".$asciiQual."\n";
    }
  }

  print OUTINDEX $fasta."\t".$barcodes[$counter]."\n";

  $counter++;
}

close(OUTINDEX);
close(OUTFASTA);
exit;     
 
