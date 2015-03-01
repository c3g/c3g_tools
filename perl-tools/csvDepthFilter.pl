#!/usr/bin/env perl
use strict;
use String::Util 'trim';

my $inputFile = $ARGV[0];
my $outputFile = $ARGV[1];
my $minDepth = $ARGV[2];

open(FILE, $inputFile) or die "Can't open $inputFile\n";
open(OFILE, '>'.$outputFile) or die "Can't write to $outputFile\n";
my $line = <FILE>;

my @values = split("\t",$line);
my @depthColumns;
for(my $i=0; $i < @values; $i++) {
  if($values[$i] =~ /Depth$/) {
    push(@depthColumns, $i);
  }
}
print OFILE $line;
while($line = <FILE>) {
  my $doPrint = 0;
  @values = split("\t",$line);
  for my $i (@depthColumns) {
    my @depths = split(',', $values[$i]);
    my $no = trim($depths[0]);
    my $tu = trim($depths[1]);
    if($no > $minDepth || $tu > $minDepth) {
      $doPrint = 1;
      last;
    }
  }

  if($doPrint != 0) {
    print OFILE $line;
  }
}
close(OFILE);
close(FILE);
