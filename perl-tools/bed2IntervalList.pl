#!/usr/bin/env perl


use strict;

use Getopt::Long;

my $version = "1.0";

&main();

sub main {
  my $dictFile;
  my $bedFile;

  my $result = GetOptions(
    "dict=s"  => \$dictFile,
    "bed=s"   => \$bedFile
  );

  my $errMsg = "";
  if(!defined($dictFile) || !-e $dictFile) {
    $errMsg .= "Missing dict\n";
  }
  if(!defined($bedFile) || !-e $bedFile) {
    $errMsg .= "Missing bed\n";
  }
  if(length($errMsg)) {
    die $errMsg;
  }

  my $line;
  my $startsWithChr=0;
  open(DICT, $dictFile) or die "Can't open dictionary: $dictFile\n";
  while($line = <DICT>) {
    if($line =~ /^\@SQ.+SN:([^\t]+)\t.*/) {
      my $chr = $1;
      if(substr($chr, 0, 3) eq 'chr') {
        $startsWithChr = 1;
      }
    }
    print $line;
  }
  close(DICT);


  open(BED, $bedFile) or die "Can't open BED: $bedFile\n";
  while($line = <BED>) {
    chomp($line);
    my @values = split(/\t/, $line);
 
    if(substr($values[0], 0, 3) eq 'chr' && $startsWithChr == 0) {
      if($values[0] eq 'chrM') {
        print 'MT';
      }
      else {
        print substr($values[0], 3);
      }
    }
    elsif(substr($values[0], 0, 3) ne 'chr' && $startsWithChr == 1) {
      if($values[0] =~ /^\d+$/ || $values[0] =~ /^[XYM]T?$/) {
        if($values[0] eq 'MT') {
          print 'chrM';
        }
        else {
          print 'chr'.$values[0]
        }
      }
      else {
        # GL
        print $values[0];
      }
    }
    else {
      # Unknown, will probably be ignored by picard
      print $values[0];
    }

    print "\t".($values[1]+1)."\t".$values[2];
    if(@values >= 6) {
      print "\t".$values[5]."\t".$values[3];
    }
    else {
      print "\t+\t";
      if(@values >= 4) {
        print $values[3];
      }
      else {
        print 'NA';
      }
    }
    print "\n";
  }
  close(BED);
}
