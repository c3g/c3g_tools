#!/usr/bin/env perl

use SampleSheet;

&main();

sub main {

  my $samplesDescriptionFile = $ARGV[0];
  my %allSamples;

  my $rAoH_SampleLaneInfos = SampleSheet::parseSampleSheet($samplesDescriptionFile);
  my %sampleInfo;
  for my $rH_Sample (@$rAoH_SampleLaneInfos) {
    if(!defined $sampleInfo{ $rH_Sample->{'name'} }) {
      $sampleInfo{ $rH_Sample->{'name'} } = [];
      $allSamples{$rH_Sample->{'name'}} = 1;
    }

    push(@{$sampleInfo{ $rH_Sample->{'name'} }}, $rH_Sample);
  }

  for(my $idx=1; $idx < @ARGV; $idx++) {
    print STDERR "Testing: " . $ARGV[$idx] . "\n";
    my ($sampleName) = $ARGV[$idx] =~ /([^\/]+)\.sorted.*.bam$/;
    my $rAoH_sample = $sampleInfo{$sampleName};
    if(!defined($rAoH_sample)) {
      warn "Missing sample in sample sheet: ".$sampleName."\n";;
    }
    delete($allSamples{$sampleName});
    

    my @header = `samtools view -H $ARGV[$idx]`;
    my %bamIDs;
    my %rgIDs;
    foreach my $line (@header){
      if($line =~ /^\@RG.*ID:([^\t]+)\t.*PU:run([^_]+)_(.)\t.*LB:([^\t]+)\t.*SM:([^\t]+)(\t.*)?/) {
        my $rgId = $1;
        my $flowcell = $2;
        my $lane = $3;
        my $library = $4;
        my $bamSampleName = $5;
        chomp($bamSampleName);
        my $key = $bamSampleName.'_'.$flowcell.'_'.$lane;
        if(defined($rgIDs{$rgId})) {
          warn "ID present multiple times. ID:".$rgId." File: ".$ARGV[$idx]."\n";
        }
        if(defined($bamIDs{$key})) {
          warn "Key present multiple times. ID:".$rgId." File: ".$ARGV[$idx]."\n";
        }

        if($bamSampleName ne $sampleName) {
          warn "Sample doesn't match. ".$bamSampleName." vs ".$sampleName." ID:".$rgId." File: ".$ARGV[$idx]."\n";
        }
        $bamIDs{$key} = $library;
        $rgIDs{$rgId} = 1;
      }
      
    }

    foreach my $rH_sample (@{$rAoH_sample}) {
      #my $key = $rH_sample->{'libraryBarcode'}.'_'.$rH_sample->{'runId'}.'_'.$rH_sample->{'lane'};
      my $key = $rH_sample->{'name'}.'_'.$rH_sample->{'runId'}.'_'.$rH_sample->{'lane'};
      if(defined($bamIDs{$key})) {
        if($rH_sample->{'libraryBarcode'} ne $bamIDs{$key}) {
          warn "Library not the same for entry: ".$key." in file ".$ARGV[$idx]." ".$rH_sample->{'libraryBarcode'}." vs ".$bamIDs{$key}."\n";
        }
        delete($bamIDs{$key});
      }
      else {

#        $key = 'A'.$rH_sample->{'name'}.'_'.$rH_sample->{'runId'}.'_'.$rH_sample->{'lane'};
#        if(defined($bamIDs{$key})) {
#          if($rH_sample->{'libraryBarcode'} ne $bamIDs{$key}) {
#            warn "Library not the same for entry: ".$key." in file ".$ARGV[$idx]."\n";
#          }
#          delete($bamIDs{$key});
#        }
#        else {
          warn "Missing entry: ".$key." in file ".$ARGV[$idx]."\n";
#        }
      }
    }
    my @remainingRGIDs = keys(%bamIDs);
    if(@remainingRGIDs > 0) {
      for my $rgId (@remainingRGIDs) {
        warn "Remaining hit: ".$rgId." in file ".$ARGV[$idx]."\n";
      }
    }
  }

  for my $remainingSample (keys(%allSamples)) {
    warn "Sample without BAM: $remainingSample\n";
  }
}


