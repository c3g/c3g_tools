#!/usr/bin/env perl

=head1 NAME

Iterator::ValidateFastq - Simple object for Fastq sequence

=head1 SYNOPSIS

    my $rec = new ValidateFastq($hdr);
    print $rec->output;

=head1 DESCRIPTION

Object for a single read sequence, with methods for basic manipulation and quality control.

=head1 METHODS

=over 5

=cut

package Iterator::ValidateFastq;

use warnings;
use strict;

=item new $hdr $seq $qual $qc_params $barcodes $variants
Initialize new sequence object. If the quality scores use Illumina scaling, the $qual_to_sanger flag *must* be set as the object
assumes and requires sanger-scaling.  Quality encoding method is determined by FastqDb class instead.
=cut

sub new
{
	my ($class, $hdr, $qual, $phred) = @_;
	die("Missing hdr\n")           unless defined($hdr);
	
	# INIT
	my $this = {
		hdr => $hdr,
		qual => $qual
		     # complete sequence without newlines
	};
	bless $this, $class;
	$this->header($hdr);      # populates id, base, pair, barcode values
	$this->convert_qual_to_sanger($qual) if($phred == 64);
	return $this;
}

=item convert_qual_to_sanger

Convert the quality string from Illumina to Sanger scaling.

=cut
sub convert_qual_to_sanger
{
	my $this = shift;
	return $this->{qual} = join('', map { chr(ord($_) - 31) } split(//, $this->{qual}));
}

=item header ($new_hdr)
Returns the object's header line.  You can also use this method to give the sequence a new header.
=cut

sub header
{
	my ($this, $hdr) = @_;
	return $this->{hdr} unless defined($hdr) and $hdr;

	$hdr = '@' . $hdr unless $hdr =~ /^@/;
	$this->{hdr}     = $hdr;
	$this->{base}    = undef;    # base ID only (e.g. "A"); always defined
	$this->{pair}    = undef;    # pair ID only (e.g. "1"); only defined if paired
	$this->{barcode} = undef;    # barcode sequence (always upper-case); only defined it barcoded

	# illumina
	if($hdr =~ /^@(\S+:\d+:\d+:\d+:\d+)#([aAtTcCgGnN]+)\/([12])/){	
		#print "1\n";
		# barcoded, paired
		$this->{base}    = $1;
		$this->{barcode} = uc($2);
		$this->{pair}    = $3;
		#print $1
	}elsif ($hdr =~ /^@(\S+:\d+:\d+:\d+:\d+)\/([12])#([aAtTcCgGnN]+)/){                            
		#print "2\n";
		# barcoded, paired
		$this->{base}    = $1;
		$this->{pair}    = $2;
		$this->{barcode} = uc($3);
	}elsif ($hdr =~ /^@(\S+:\d+:\d+:\d+:\d+)\/([12])/){ 
		#print "3\n";
		# paired
		$this->{base} = $1;
		$this->{pair} = $2;
	}elsif ($hdr =~ /^@(\S+:\d+:\d+:\d+:\d+)#([aAtTcCgGnN]+)/){ 
		#print "4\n";
		# barcoded, unpaired
		$this->{base}    = $1;
		$this->{barcode} = uc($2);
	}elsif ($hdr =~ /^@(\S+:\d+:\d+\S+:\d+:\d+:\d+:\d+) ([12]):\S:\d:([aAtTcCgGnN]+)/){
		#print "5\n";
		$this->{base}    = $1;
		$this->{pair}    = $2;	
		$this->{barcode} = uc($3);
	}elsif ($hdr =~ /^@(\S+:\d+:\S+:\d+:\d+:\d+:\d+) ([12]):\S:\S:([aAtTcCgGnN]+)/){
		#print "6\n";
		$this->{base}    = $1;
		$this->{pair}    = $2;
		$this->{barcode} = uc($3);
		#print $2."\n";
	}elsif($hdr =~ /^@(\S+:\d+:\S+:\d+:\d+:\d+:\d+) (\d):\S:\d:/){
			#@HISEQ08:232:C15DJACXX:6:1101:1072:2149 1:N:0:  65  
    		#@HISEQ08:232:C15DJACXX:6:1101:1235:2140 1:N:0:
			$this->{base}    = $1;
			$this->{pair}    = $2;	
	}elsif ($hdr =~ /^@(\S+:\d+:\S+:\d+:\d+:\d+:\d+) ([12]):\S:\S:([aAtTcCgGnN]+)/){
		#print "7\n";
		$this->{base}    = $1;
		$this->{pair}    = $2;	
		$this->{barcode} = uc($3);
	}elsif ($hdr =~ /^@(\S+:\d+:\d+:\d+:\d+) ([12]):\S:\S:([aAtTcCgGnN]+)/){
		#@HWI-ST753_98:3:1101:1590:1994 1:N:0:CGCCTTAAGGA
		$this->{base}    = $1;
		$this->{pair}    = $2;	
		$this->{barcode} = uc($3);
	}elsif ($hdr =~ /^@(\S+:\d+:\d+:\d+:\d+)/){                            
		#print "8\n";
		# unpaired
		$this->{base} = $1;
	}elsif($hdr =~ /^@(\S+) length=\d+ xy=\d+_\d+ region=\d run=.*#([aAtTcCgGnN]+)\/([12])/ ){
		#print "9\n";
		$this->{base}    = $1;
		$this->{pair}    = $3;	
		$this->{barcode} = uc($2);
	}elsif ($hdr =~ /^@(\S+) length=\d+ xy=\d+_\d+ region=\d run=\S/){
		#print "10\n";
		# Roche/454
		$this->{base} = 1;
		#Special fastqs from Dangl lab
	}elsif ($hdr =~ /^@(\d+:[aAtTcCgGnN]+_[aAtTcCgGnN]+)#([aAtTcCgGnN]+)/){
		#print "11\n";
		$this->{base} = $1;
		$this->{barcode} = uc($2);
	}elsif ($hdr =~ /^@(\S+)#([aAtTcCgGnN]+)\/([12])/){
		#print "12\n";
		$this->{base}    = $1;
		$this->{pair}    = $3;	
		$this->{barcode} = uc($2);	
	}elsif ($hdr =~ /^@(\S+)#([aAtTcCgGnN]+)/){
		#print "13\n";
        # Roche/454 with barcode
        $this->{base} = $1;
        $this->{barcode} = uc($2);
	}elsif ($hdr =~ /^@(\S+)/){
		#print "14\n";
		# generic fastq (not an Illumina or Roche read header)
		$this->{base} = $1;
	}else{
		die("Unable to parse sequence header: $hdr\n");
	}

	return $hdr;
}

=item id
Returns the object's ID, which is the sequence's unique identifier without comments which may be present in the header.
It cannot be changed directly; set header, base, barcode, or pair instead.
=cut
sub id
{
	my $this = shift;
	my $id   = $this->{base};
	$id .= "#" . $this->{barcode} if $this->{barcode};
	$id .= "/" . $this->{pair}    if $this->{pair};
	return $id;
}

=item base
Returns the base of the ID; both paired reads will have the same base.
=cut
sub base { return shift->{base} }

=item qual
Returns the qual in sanger format; both paired reads will have the same base.
=cut
sub qual { return shift->{qual} }

=item barcode
If the read contains an Illumina molecular barcode ID, it will be returned; otherwise returns undef.
Supplying an optional argument will set the barcode.
=cut
sub barcode
{
	my ($this, $barcode) = @_;
	if (defined($barcode))
	{
		$this->{barcode} = $barcode ? $barcode : undef;
	}
	return $this->{barcode};
}

=item pair
Returns the read's ord in the pair; undef otherwise.
=cut
sub pair { return shift->{pair} }
1;
