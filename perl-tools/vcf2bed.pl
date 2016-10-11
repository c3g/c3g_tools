#!/usr/bin/perl

use strict;

#if (scalar @ARGV != 1) {
#    print "vcf2bed.pl <vcf_file>\n";
#    exit 1;
#}
    
while (<>) {
    
    if (index($_, '#') == 0) {
        next;
    }
    
    my @fields = split /\t/, $_;
    
    my $chrom = $fields[0];
    my $startPos = $fields[1];
    my $ref = $fields[3];
    my @alts = split (/,/, $fields[4]);
    
    my $endPos = $startPos;
    $startPos = $startPos;
    #$startPos = $startPos - 1;
    
    my $indelCall = "";
        
    foreach my $alt (@alts) {
    
        if (length($ref) > length($alt)) {
            # deletion event
            my $diff = substr($ref, length($alt));
        
            $endPos = $startPos + length($diff);
            #$indelCall = '-' . $diff;
        }
        elsif (length($alt) > length($ref)) {
            # insertion event
            my $diff = substr($alt, length($ref));
        
            #$indelCall = '+' . $diff;
        }
        
        # print out the line

        print "$chrom\t$startPos\t$endPos\n"
        #print "$chrom\t$startPos\t$endPos\t$ref\t$alt\n"
        #print "$chrom\t$startPos\t$endPos\t$ref\t$alt\t$indelCall\n";
    }
    
    
}   

close VCF;
