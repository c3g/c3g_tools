#!/usr/bin/perl

&main();

sub main {

  my $samplesDescriptionFile = $ARGV[0];
  my %allSamples;

  my $rAoH_SampleLaneInfos = parseSheet($samplesDescriptionFile);
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
      if($line =~ /^@RG.*ID:([^\t]+)\t.*PU:run([^_]+)_(.)\t.*LB:([^\t]+)\t.*SM:([^\t]+)(\t.*)?/) {
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
#      my $key = $rH_sample->{'libraryBarcode'}.'_'.$rH_sample->{'runId'}.'_'.$rH_sample->{'lane'};
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
    warn "Sample wuthout BAM: $remainingSample\n";
  }
}

sub parseSheet {
  my $fileName = shift;

  my @retVal;
  open(SAMPLE_SHEET, "$fileName") or die "Can't open $fileName\n";
  my $line = <SAMPLE_SHEET>;
  my @headers = split(/",/, $line);
  my $nameIdx=-1;
  my $libraryBarcodeIdx=-1;
  my $runIdIdx=-1;
  my $laneIdx=-1;
  my $runTypeIdx=-1;
  my $statusIdx=-1;
  my $qualOffsetIdx=-1;
  for(my $idx=0; $idx < @headers; $idx++) {
    $headers[$idx] =~ s/"//g;
    if($headers[$idx] eq "Name") {
      $nameIdx=$idx;
    }
    elsif($headers[$idx] eq "Library Barcode") {
      $libraryBarcodeIdx=$idx;
    }
    elsif($headers[$idx] eq "Run") {
      $runIdIdx=$idx;
    }
    elsif($headers[$idx] eq "Region") {
      $laneIdx=$idx;
    }
    elsif($headers[$idx] eq "Run Type") {
      $runTypeIdx=$idx;
    }
    elsif($headers[$idx] eq "Status") {
      $statusIdx=$idx;
    }
    elsif($headers[$idx] eq "Quality Offset") {
      $qualOffsetIdx=$idx;
    }
  }

  my $sampleSheetErrors="";
  if($nameIdx==-1) {
    $sampleSheetErrors.="Missing Sample Name\n";
  }
  if($libraryBarcodeIdx==-1) {
    $sampleSheetErrors.="Missing Library Barcode\n";
  }
  if($runIdIdx==-1) {
    $sampleSheetErrors.="Missing Run ID\n";
  }
  if($laneIdx==-1) {
    $sampleSheetErrors.="Missing Lane\n";
  }
  if($runTypeIdx==-1) {
    $sampleSheetErrors.="Missing Run Type\n";
  }
  if($statusIdx==-1) {
    $sampleSheetErrors.="Missing Status\n";
  }
  if($qualOffsetIdx==-1) {
    $sampleSheetErrors.="Missing Quality Offset\n";
  }
  if(length($sampleSheetErrors) > 0) {
    die $sampleSheetErrors;
  }

  while($line = <SAMPLE_SHEET>) {
    $line =~ s/"//g;
    my @values = split(/,/, $line);
    if($values[$statusIdx] ne 'Data is valid') {
      warn "Invalid: $values[$nameIdx] $values[$runIdIdx] $values[$laneIdx]\n";
      next;
    }

    my %sampleInfo;
    $sampleInfo{'name'} = $values[$nameIdx];
    $sampleInfo{'libraryBarcode'} = $values[$libraryBarcodeIdx];
    $sampleInfo{'runId'} = $values[$runIdIdx];
    $sampleInfo{'lane'} = $values[$laneIdx];
    $sampleInfo{'runType'} = $values[$runTypeIdx];
    $sampleInfo{'qualOffset'} = $values[$qualOffsetIdx];

    if($values[$runTypeIdx] eq "PAIRED_END") {
      $sampleInfo{'read1File'} = $sampleInfo{'name'}.'.'.$sampleInfo{'qualOffset'}.".pair1.fastq.gz";
      $sampleInfo{'read2File'} = $sampleInfo{'name'}.'.'.$sampleInfo{'qualOffset'}.".pair2.fastq.gz";
    }
    elsif($values[4] eq "SINGLE_END") {
      $sampleInfo{'read1File'} = $sampleInfo{'name'}.'.'.$sampleInfo{'qualOffset'}.".single.fastq.gz";
    }
    else {
      print "Unrecognized run type $values[$runTypeIdx] \n";
      exit 1;
    }

    push(@retVal, \%sampleInfo);
  }

  return \@retVal;
}
