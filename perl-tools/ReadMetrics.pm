#!/usr/bin/env perl

=head1 NAME

I<ReadMetrics>

=head1 SYNOPSIS

ReadMetrics-> parseFlagstats()
ReadMetrics-> parseTrimOutput()
ReadMetrics-> mergeStats()

=head1 DESCRIPTION

B<ReadMetrics> is a library to generate QC, stats and metrics


=head1 AUTHOR
B<Johanna Sandoval> - I<johanna.sandoval@mail.mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package ReadMetrics;

# Strict Pragmas
#--------------------------
use strict;
use warnings;


#--------------------------

# Dependencies
#-----------------------
#use LoadConfig;

# SUB
#-----------------------

sub parseFlagstats{
	my $sampleName    = shift;
	my $flagStatsFile = shift;
  my %stats         ;
  my $nbQCPassedReads = 0;
  my $nbDuplicateReads= 0;
  my $nbAlignedReads  = 0;
  open(FLAGSTATS, "$flagStatsFile") or return;
  $stats{"_HEADER_"}=["Number of QC Passed Reads", "Number of Aligned Reads", "Number of Duplicate Reads"];
  while(my $line = <FLAGSTATS>) {
    chomp($line);
    if($line =~ /^([0-9]+) \+ [0-9]+ in total/) {
       $nbQCPassedReads = $1;
    }
    elsif($line =~ /^([0-9]+) \+ [0-9]+ duplicates/) {
      $nbDuplicateReads = $1;
    }
    elsif($line =~ /^([0-9]+) \+ [0-9]+ mapped/) {
      $nbAlignedReads = $1;
    }
  }
  $stats{$sampleName}=[$nbQCPassedReads, $nbAlignedReads, $nbDuplicateReads];
  close(FLAGSTATS);
  return \%stats;
}

sub parseTrimOutput{
  my $sampleName      = shift;
  my $trimOutputFile  = shift;
  my %stats           ;
  my $nbRawReads      = -1;
  my $nbFilteredReads = -1;
  my $nbSingleFiltered= -1;
  
  $stats{"_HEADER_"}=["Raw Fragments", "Fragment Surviving", "Single Surviving"];
  open(TRIMSTATS, "$trimOutputFile") or return;

  while(my $line = <TRIMSTATS>) {
    chomp($line);
    if($line =~ /^Raw Fragments\,([0-9]+)/){
      $nbRawReads = $1;
    }elsif($line =~ /^Fragment Surviving\,([0-9]+)/ ){
      $nbFilteredReads = $1;
    }elsif($line =~ /^Single Surviving\,([0-9]+)/ ){
      $nbSingleFiltered= $1;
    }    
  }
  $stats{$sampleName}=[$nbRawReads, $nbFilteredReads, $nbSingleFiltered];
  close(TRIMSTATS);  
  return \%stats;
}

sub parseHomerAnnotations {
  my $annotationFile  = shift;
  my $outputFile      = shift;
  my $proximalDistance= shift;
  my $distalDistance = shift;
  my $distance5d1    = shift;
  my $distance5d2    = shift;
  my $geneDesertSize = shift;
  
  my $nbExon=0;
  my $nbIntron=0;
  my $nbProximal=0;
  my $nbDistal=0;
  my $nb5d=0;
  my $nbgeneDesert=0;
  my $nbOther=0;
  my $nbAnnotated=0;
  my @nbExonsperTranscript=();
  my @nbIntronsperTranscript=();
  my @AoTssDistances=();
    
  # Parse annotations.csv
  
  open(ANNOTATIONS, "$annotationFile") or return;
  while(my $line = <ANNOTATIONS>) {
    # skip header
    next if $. < 2;
    chomp($line);
    my @words         = split( /\t/, $line);
    # skip line without distance to transcription factor
    next if $words[9] =~ /NA/ ;
    # Compute distance to TSS
    my $distanceToTSS = $words[9];
    $nbAnnotated +=1;
    push(@AoTssDistances, $distanceToTSS); 
    if( $words[7] =~ /^exon / ){
      $nbExon += 1;
      my @wsplit=split(" ", $words[7]);
      $wsplit[$#wsplit] =~ s/\)//g;
      push(@nbExonsperTranscript, $wsplit[$#wsplit] ); 
      
    }elsif( $words[7] =~ /^intron / ){
      $nbIntron += 1;
      my @wsplit=split(" ", $words[7]);
      $wsplit[$#wsplit] =~ s/\)//g;
      push(@nbIntronsperTranscript, $wsplit[$#wsplit] ); 
      
    }else{
      if( $distanceToTSS < 0 && $distanceToTSS > $proximalDistance  ){
        #print 'proximal '. $distanceToTSS. "\n";
        $nbProximal += 1;
      }elsif( $distanceToTSS <= $proximalDistance && $distanceToTSS > $distalDistance ){
        #print 'distal '. $distanceToTSS. "\n";
        $nbDistal += 1;
      }elsif( $distanceToTSS <= $distance5d1 && $distanceToTSS > $distance5d2){
        #print '5d '. $distanceToTSS. "\n";
        $nb5d += 1;
      }elsif( $distanceToTSS < (-1*$geneDesertSize) || $distanceToTSS > $geneDesertSize ){
        #print 'gene_desert '. $distanceToTSS. "\n";
        $nbgeneDesert += 1;
      }else{
        #print 'other '. $distanceToTSS. "\n";
        # do nothing
      }      
    }    
  }
  $nbOther = $nbAnnotated - $nbExon -   $nbIntron- $nbProximal- $nbDistal - $nb5d - $nbgeneDesert;
  close(ANNOTATIONS);
  # Print general stats
  open(OUTFILE, ">".$outputFile.".tss.stats.csv");
  print OUTFILE "exon,intron,proximal,distal,5d,gene_desert,other"."\n";
  print OUTFILE  $nbExon .",". $nbIntron .",". $nbProximal .",". $nbDistal .",". $nb5d .",". $nbgeneDesert .",". $nbOther  . "\n";
  close(OUTFILE);
  # Print exon stats
  open(OUTFILE, ">".$outputFile.".exon.stats.csv");
  print OUTFILE join ("\n", @nbExonsperTranscript). "\n";
  close(OUTFILE);
  # Print intron stats
  open(OUTFILE, ">".$outputFile.".intron.stats.csv");
  print OUTFILE join("\n", @nbIntronsperTranscript). "\n";
  close(OUTFILE);
  # Print intron stats
  open(OUTFILE, ">".$outputFile.".tss.distance.csv");
  print OUTFILE join("\n", @AoTssDistances). "\n";
  close(OUTFILE);
}

sub mergeStats{
  my $sampleName           = shift;
  my $outputFile           = shift;
  my $rHoA_trimStats       = shift;
  my $rHoA_alignmentStats  = shift;
  my $printHeader          = shift;
  
  # Print a header by default
  if (!defined($printHeader) ){
    $printHeader=1;
  }
  # Open an output file
  open(OUTFILE, ">".$outputFile) or die "ERROR : Unable to open". $outputFile;
  
  if (defined($rHoA_trimStats) && defined($rHoA_alignmentStats)){
    if (defined($printHeader) ){
     print OUTFILE "SampleName".",".join(',', @{$rHoA_trimStats->{"_HEADER_"}}). ','. join(',', @{$rHoA_alignmentStats->{"_HEADER_"}})."\n";
    }
    print OUTFILE $sampleName.",".join(',', @{$rHoA_trimStats->{$sampleName}}). ','. join(',', @{$rHoA_alignmentStats->{$sampleName}})."\n";
  }elsif(defined($rHoA_trimStats) && !defined($rHoA_alignmentStats)){
    if (defined($printHeader) ){
     print OUTFILE "SampleName".",".join(',', @{$rHoA_trimStats->{"_HEADER_"}})."\n";
    }
    print OUTFILE $sampleName.",".join(',', @{$rHoA_trimStats->{$sampleName}})."\n";

  }elsif(!defined($rHoA_trimStats) && defined($rHoA_alignmentStats)){
    if (defined($printHeader) ){
      print OUTFILE "SampleName".",".join(',', @{$rHoA_alignmentStats->{"_HEADER_"}})."\n";
    }
    print OUTFILE $sampleName.",". join(',', @{$rHoA_alignmentStats->{$sampleName}})."\n";
  }else{
    print OUTFILE "\n";
  }
	close(OUTFILE);
}

1;
