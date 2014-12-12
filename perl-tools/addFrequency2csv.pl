#!/usr/bin/env perl

use Time::localtime;
use strict;

my $EXPERIMENT_TYPE = $ENV{EXPERIMENT_TYPE};

&main();


sub initVcfFiles {
  my $rA_vcfFiles = shift;  
  my $rH_normals = shift;
  my $rH_tumors = shift;

  my %fileDetails;
  my %samples;
  for my $vcfFile (@{$rA_vcfFiles}) {
    #my ($sampleName) = $vcfFile =~ /^([^.\/]+)/;
    my ($sampleName) = $vcfFile =~ /.*alignment\/(.*)\/.*/;
    if(!defined($sampleName) || length($sampleName) == 9) {
      die "Couldn't extract sample name from : $vcfFile\n";
    }

    my $fh;
    open($fh, "<", $vcfFile) or die "Can't open $vcfFile\n";
    my $line = <$fh>;
    $fileDetails{$sampleName} = {};
    $fileDetails{$sampleName}->{'fh'} = $fh;
    $line = <$fh>;
    chomp($line);
    my @values = split(",", $line);
    $fileDetails{$sampleName}->{'currentLine'} = $line;
    $fileDetails{$sampleName}->{'split'} = \@values;
  }
  return \%fileDetails;
}

sub main {
  my ($variantFile, $normalTumorFile, @vcfFiles) = @ARGV;

  open(VARFILE, $variantFile) or die "Can't open $variantFile\n";
  my ($rH_normals, $rH_tumors, $rHoA_alias) = parseSampleNormalTumor($normalTumorFile);

  my $rHoH_samplesFiles = initVcfFiles(\@vcfFiles, $rH_normals, $rH_tumors);
  my @samplesGiven = keys(%{$rHoH_samplesFiles});

  my $line = <VARFILE>;
  my @headers = split("\t", $line);
  my %columnIdxToSample;
  my $shiftOffset;
  print $headers[0];
  for(my $i=1; $i < @headers; $i++) {
    print "\t";
    if(defined($rHoA_alias->{$headers[$i]})) {
      $columnIdxToSample{$i} = $headers[$i];
      print $headers[$i]." Allele Frequency\t";
      print $headers[$i]." Variant normal ratio\t";
      print $headers[$i]." Variant tumor ratio\t";
    }
    print $headers[$i];
  }
  
  while($line = <VARFILE>){
    my @values = split("\t", $line);
    my $chr = $values[0];
    my $pos = $values[1];
    my $ref = $values[2];
    my $alt = $values[3];
#print STDERR $chr.",".$pos."\n";
    print $values[0];
    for(my $i=1; $i < @values; $i++) {
      print "\t";
      if(defined($columnIdxToSample{$i})){
        my $sample = $columnIdxToSample{$i};
        my $rA_normalVarCounts = printSampleFreq($chr, $pos, $ref, $alt, $rHoA_alias->{$sample}->[0], $rHoH_samplesFiles);
        print ',';
        my $rA_tumorVarCounts = printSampleFreq($chr, $pos, $ref, $alt, $rHoA_alias->{$sample}->[1], $rHoH_samplesFiles);
        print "\t";
        if(defined($rA_normalVarCounts) && defined($rA_tumorVarCounts) && ($rA_normalVarCounts->[0]+$rA_normalVarCounts->[1]) > 0 && ($rA_tumorVarCounts->[0]+$rA_tumorVarCounts->[1]) > 0) {
          #print log(($rA_normalVarCounts->[0]/$rA_normalVarCounts->[1])/($rA_tumorVarCounts->[0]/$rA_tumorVarCounts->[1]))/log(2);
          print ($rA_normalVarCounts->[1]/($rA_normalVarCounts->[0]+$rA_normalVarCounts->[1]));
          print "\t";
          print ($rA_tumorVarCounts->[1]/($rA_tumorVarCounts->[0]+$rA_tumorVarCounts->[1]));
        }
        else {
          print "NA\tNA";
        }
        print "\t";

      }
      print $values[$i];
    }
  }
  close(VARFILE);
}

sub printSampleFreq {
  my $chr = shift;
  my $pos = shift;
  my $ref = shift;
  my $alt = shift;
  my $sample = shift;
  my $rHoH_samplesFiles = shift;

  my $rH_fileDetails = $rHoH_samplesFiles->{$sample};
  while($rH_fileDetails->{'split'}->[0] ne $chr || $rH_fileDetails->{'split'}->[1] != $pos) {
    my $fh = $rH_fileDetails->{'fh'};
    my $newLine = <$fh>;
    if(!defined($newLine)){
      close($rH_fileDetails->{'fh'});
    }
    else {
      chomp($newLine);
      my @newValues = split(",", $newLine);
      $rH_fileDetails->{'currentLine'} = $newLine;
      $rH_fileDetails->{'split'} = \@newValues;
    }
  }

  if($rH_fileDetails->{'split'}->[0] ne $chr || $rH_fileDetails->{'split'}->[1] != $pos) {
    die "\n\nWrong $sample $chr $pos ".$rH_fileDetails->{'split'}->[0]." ".$rH_fileDetails->{'split'}->[1]."\n\n";
  }

  print $rH_fileDetails->{'split'}->[2];

  if(length($ref) > 1 || length($alt) > 1){
    return undef;
  }

  my @bases = split(/ /, $rH_fileDetails->{'split'}->[2]);
  my @retVal;
  for my $baseCnt (@bases) {
    my $base = substr($baseCnt, 0, 1);
    my $count = substr($baseCnt, 2);

    if($base eq $ref) {
      $retVal[0] = $count;
    }
    elsif($base eq $alt) {
      $retVal[1] = $count;
    }
  }

  if(@retVal < 1) {
    die "Not found: $chr:$pos\n";
  }
  return \@retVal;
}

sub parseSampleNormalTumor {
  my $normalTumorFile = shift;
  my %normals;
  my %tumors;
  my %alias;

  open(SAMPLE, "$normalTumorFile") or die "Can't open sample normal tumor file: $normalTumorFile\n";
  while(my $line =<SAMPLE>) {
    chomp $line;
    if($line =~ /^#/) {
      next;
    }
    my @values = split(/,/, $line);
    $normals{$values[1]} = $values[0];
    $tumors{$values[2]} = $values[0];
    $alias{$values[0]} = [$values[1],$values[2]];
  }
  close(SAMPLE);

  return (\%normals, \%tumors, \%alias);
}
