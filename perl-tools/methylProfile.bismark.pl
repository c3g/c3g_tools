#!/usr/bin/env perl
use warnings;
################################################################################
## file name: methylProfile.bismark.pl
## by Xiaojian  @Dec 19, 2016
## description: get all CpGs profile by combining both forward and reverse strands
##
################################################################################

use strict;
use FileHandle;
use Env;
use Getopt::Long;
use PerlIO::gzip;

my $help = "";
my $input = "";
my $output = "";
my $total_reads = 0;

GetOptions (
    'help'       => \$help,
    'i|input=s'  => \$input,
    'o|output=s' => \$output,
);

if ($help) {
    print_helpfile();
   exit();
}

die "--input missing\n" unless($input);
die "--output missing\n" unless($output);

die "Input file $input does not exist...\n" unless(-e $input);

my $temp = "";
my @temparray = ();
my %CpGprofileHash = ();
my %CpGfwHash = ();
my %CpGrvHash = ();
my %CpGfwReadsHash = ();
my %CpGrvReadsHash = ();
my %CpGfwMethReadsHash = ();
my %CpGrvMethReadsHash = ();

my $strand = "";
my $fwcoordinate = 0;
my $CpGid = "";
my %CpGtotalReadsHash = ();
my %CpGtotalMethReadsHash = ();
my $total_meth = 0;
my $fw_meth = 0;
my $rv_meth = 0;
my $chr = "";

open INFILE, "<:gzip", $input || die "Cannot open $input for reading\n";
open OUTFILE, ">$output" || die "Cannot open $output for writing.\n";
print OUTFILE "#chr\tstart\tend\tnum_C_fw\tnum_total_fw\tmeth_fw\tnum_C_rv\tnum_total_rv\tmeth_rv\ttotal_C\ttotal\ttotal_meth\n";

while ($temp = <INFILE>) {
    chomp $temp;
    @temparray = split("\t", $temp); ## tab separated
    $strand = $temparray[2];
    $chr = $temparray[0];

    if ($strand eq "+") {
        $CpGid = $temparray[0] . "." . $temparray[1];
        $CpGfwReadsHash{$CpGid} = ($temparray[3] + $temparray[4]);
        $CpGfwMethReadsHash{$CpGid} = $temparray[3];

        if ($temparray[3] == 0) {
            $fw_meth = 0;
        } else {
            $fw_meth = sprintf("%.2f", (($CpGfwMethReadsHash{$CpGid}/$CpGfwReadsHash{$CpGid})*100));
        }
    } else {
        $fwcoordinate = $temparray[1] - 1;
        $CpGid = $temparray[0] . "." . $fwcoordinate;
        $CpGrvReadsHash{$CpGid} = ($temparray[3] + $temparray[4]);
        if (exists $CpGfwReadsHash{$CpGid}) {
            $CpGtotalReadsHash{$CpGid} = ($CpGfwReadsHash{$CpGid} + $CpGrvReadsHash{$CpGid});
            $CpGtotalMethReadsHash{$CpGid} = ($CpGfwMethReadsHash{$CpGid} + $temparray[3]);
        } else {
            $CpGtotalReadsHash{$CpGid} = $CpGrvReadsHash{$CpGid};
            $CpGtotalMethReadsHash{$CpGid} = $temparray[3];
            $CpGfwReadsHash{$CpGid} = 0;
            $CpGfwMethReadsHash{$CpGid} = 0;
            $fw_meth = 0;
        }

        if ($temparray[3] == 0) {
            $rv_meth = 0;
        } else {
            $rv_meth = sprintf("%.2f", (($temparray[3]/$CpGrvReadsHash{$CpGid})*100));
        }

        if ($CpGtotalReadsHash{$CpGid} > 0) {
            $total_meth = sprintf("%.2f", (($CpGtotalMethReadsHash{$CpGid}/$CpGtotalReadsHash{$CpGid})*100));
            print OUTFILE "$chr\t$fwcoordinate\t$temparray[1]\t$CpGfwMethReadsHash{$CpGid}\t$CpGfwReadsHash{$CpGid}\t$fw_meth\t$temparray[3]\t$CpGrvReadsHash{$CpGid}\t$rv_meth\t$CpGtotalMethReadsHash{$CpGid}\t$CpGtotalReadsHash{$CpGid}\t$total_meth\n";
            %CpGfwHash = ();
            %CpGrvHash = ();
            %CpGfwReadsHash = ();
            %CpGrvReadsHash = ();
            %CpGfwMethReadsHash = ();
            %CpGrvMethReadsHash = ();
        }
    }
}
close(INFILE);
print "Finish reading file $input\n";
close(OUTFILE);

sub print_helpfile{
    print "\n", '='x100, "\n";
    print ">>> USAGE: perl methylProfile.bismark.pl [Options] <<<\n\n";
    print "Options -i|--input: input file of bismark CpG report gz file.\n";
    print "Options -o|--outfile: output file of combined CpG profile from both forward and reverse strands.\n";
    print "\n>>> Example: perl methylProfile.bismark.pl -i sample.CpG_report.txt.gz -o sample.cpg.profile.strand.combined.csv\n\n\n";
}
