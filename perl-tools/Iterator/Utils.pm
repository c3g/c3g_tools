#!/usr/bin/env perl

=head1 NAME

Iterator::Utils - Simple object for Fastq sequence

=head1 SYNOPSIS

    my $rec = new Utils($hdr);
    print $rec->output;

=head1 DESCRIPTION

Object for a single read sequence, with methods for basic manipulation and quality control.

=head1 METHODS

=over 5

=cut

package Iterator::Utils;

use warnings;
use strict;

=item new $file
Initialize new sequence object. If the quality scores use Illumina scaling, the $qual_to_sanger flag *must* be set as the object
assumes and requires sanger-scaling.  Quality encoding method is determined by FastqDb class instead.
=cut

sub new
{
	my ($class, $file) = @_;
	die("Missing file\n")           unless defined($file);
	
	# INIT
	my $this = {
		file => $file,
		     # complete sequence without newlines
	};
	bless $this, $class;
	my $qual = check_qual($file);
	#$this->check_qual($file);      # populates id, base, pair, barcode values
	return $qual;
}

sub check_qual
{
	my $file = shift;
	#print "file:\t".$file."\n";

	# INIT VARS
	my $appears_sanger   = 0;
	my $appears_illumina = 0;
	my $phred;

	open my $FASTQ, '<', $file or die $!;

	my $counter = 0;
	
	while( $counter < 5000 ) {
		my @lines = map scalar( <$FASTQ> ), 1 .. 4;
		last unless defined $lines[3];
		chomp @lines;
		my @qual = split(//, $lines[3]);
		foreach my $q (@qual)
		{
			my $ascii = ord($q);
			if    ($ascii < 66) { ++$appears_sanger }
			elsif ($ascii > 74) { ++$appears_illumina }
		}
		$counter++;
	}
	# SAVE RESULT
	if($appears_sanger > $appears_illumina){
		$phred = 33;
	}else{
		$phred = 64;
	}
	die "Can't determine if quality scale is phred 33 or 64...\n" if($appears_sanger == $appears_illumina);
	
	return $phred;
}
1;
