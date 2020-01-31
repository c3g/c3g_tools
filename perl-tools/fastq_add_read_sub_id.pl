#!/usr/bin/env perl

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Script adding a read indicator (/1 or /2) at the end of the first and third line of each read entry

for (my $i = 1; $i <= 2; $i++) {
    # get arguments by pairs incrementing the read counter (max 2). The order of input must fit the read order.
    $infile = shift @ARGV;
    $outfile = shift @ARGV;
    print "Outputfile: $outfile\n";

    unlink($outfile);
    open(OUTPUT, '>>', $outfile) or die "Error opening $outfile : $!\n";
    open(INPUT, $infile) or die "Error opening $infile : $!\n";
    while  ( my $line = <INPUT>) {
        chomp($line);
        my @tmp1 = split(' ', $line);
        my $read1 = join('/', $tmp1[0], $i); # Add the read indicator at the end of the first line of the read entry

        $line = <INPUT>;
        chomp($line);
        my $seq1 = $line;

        $line = <INPUT>;
        chomp($line);
        my @tmp2 = split(' ', $line);
        my $read2 = join('/', $tmp2[0], $i); # Add the read indicator at the end of the third line of the read entry

        $line = <INPUT>;
        chomp($line);
        my $qual1 = $line;

        print OUTPUT "$read1\n$seq1\n$read2\n$qual1\n"; # Format of the read entry. Indicators are in read1 and read2.
    }
    close INPUT;
    close OUTPUT;
}

exit 0;
