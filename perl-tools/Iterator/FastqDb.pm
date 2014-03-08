#!/usr/bin/env perl

=head1 NAME

Iterator::Fastq - Simple object for Fastq sequence

=head1 SYNOPSIS

    my $rec=new Fastq( $hdr, $id, $base, $barcode, $pair, $seq, $qual, $qc_params, $barcodes, $barcode_variants);
    print $rec->output;

=head1 DESCRIPTION

Object for a single read sequence, with methods for basic manipulation and quality control.

=head1 METHODS

=over 5

=cut

package Iterator::Fastq;

use warnings;
use strict;
use constant {
	CHARACTERS_PER_LINE    => 80,    # for formatting Fasta/Qual output only
	CLOSE_ENOUGH_TO_END    => 6,     # hits of adapters this close to end will result in trimming to be done to end
	DEFAULT_MINLEN         => 20,
	DEFAULT_MEANQ          => 20,
	DEFAULT_WINSIZE        => 5,
	DEFAULT_LOW_COMPLEXITY => 0.8,
	DEFAULT_MAXN           => 3,};

=item new $hdr $seq $qual $qc_params $barcodes $variants

Initialize new sequence object. If the quality scores use Illumina scaling, the $qual_to_sanger flag *must* be set as the object
assumes and requires sanger-scaling.  Quality encoding method is determined by FastqDb class instead.

=cut

sub new
{
	my ($class, $hdr, $seq, $qual, $qc_params, $barcodes, $barcode_variants) = @_;
	die("Missing hdr\n")           unless defined($hdr);
	die("Missing seq for $hdr\n")  unless defined($seq);
	die("Missing qual for $hdr\n") unless defined($qual);
	if (length($seq) != length($qual))
	{
		warn("Seq and qual don't match for: $hdr\n");
		return undef;
	}
	# INIT
	my $this = {
		seq      => $seq,     # complete sequence without newlines
		qual     => $qual,    # complete quality string without newlines
		filtered => undef     # reason why filtered
	};
	bless $this, $class;
	$this->header($hdr);      # populates id, base, pair, barcode values
	$this->convert_qual_to_sanger if $qc_params->{qual_to_sanger};
	$this->trim_roche_mid($qc_params->{roche_mid_len}) if defined($qc_params->{roche_mid_len});
	$this->check_barcode($barcodes, $barcode_variants) if defined($barcodes);
	$this->qc($qc_params);
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

=item trim_roche_mid

Trim the first x bases which are the 454 MID.

=cut

# TRIM THE ROCHE MOLECULAR ID FROM THE 5' END OF THE SEQUENCE
sub trim_roche_mid
{
	my ($this, $len) = @_;
	return unless defined($len);
	die("Invalid MID length, $len\n") unless $len > 0;
	if (length($this->{seq}) <= $len)
	{
		$this->{seq} = $this->{qual} = undef;
		return;
	}
	$this->{barcode} = uc(substr($this->{seq}, 0, $len));
	$this->{seq}  = substr($this->{seq},  $len);
	$this->{qual} = substr($this->{qual}, $len);
}

=item check_barcode

Check if barcode is valid, perform 1-base error correction, or filter read.

=cut

sub check_barcode
{
	my ($this, $barcodes, $barcode_variants) = @_;
	return unless defined($barcodes) and defined($barcode_variants);
	my $barcode = $this->barcode;
	return unless $barcode;
	return if exists($barcodes->{$barcode});
	if (exists($barcode_variants->{$barcode}))
	{
		$this->barcode($barcode_variants->{$barcode});
	} else
	{
		$this->del('invalid_barcode');
	}
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

=item unpair

Method to clear pairing of a read (e.g. singleton).

=cut

sub unpair { shift->{pair} = undef }

=item seq

Returns the read's complete sequence, without newlines.

=cut

sub seq
{
	my $this = shift;
	return $this->{filtered} ? '' : $this->{seq};
}

=item len

Returns the length of the sequence (returns 0 if filtered).

=cut

sub len
{
	my $this = shift;
	return $this->{filtered} ? 0 : length($this->{seq});
}

=item revcomp

Reverse-complements a sequence and quality scores.

=cut

sub revcomp
{
	my $this = shift;
	return unless $this->{seq};
	$this->{seq} =~ tr/ATCGatcg/TAGCtagc/;
	my @seq = reverse split(//, $this->{seq});
	$this->{seq} = join('', @seq);
	my @qual = reverse split(//, $this->{qual});
	$this->{qual} = join('', @qual);
}

=item qual ($new_qual)

Returns the read's quality string, without newlines.

=cut

sub qual
{
	my $this = shift;
	return $this->{filtered} ? '' : $this->{qual};
}

=item qual_arrayref

Returns an arrayref of sanger-scaled quality scores.

=cut

sub qual_arrayref
{
	my $this = shift;
	return [] unless $this->{qual};
	return [] if $this->{filtered};
	my @qual = map { ord($_) - 33 } split(//, $this->{qual});
	return \@qual;
}

=item del

Mark this record as filtered.

=cut

sub del
{
	my ($this, $reason) = @_;
	die("Missing filter reason\n") unless $reason;
	$this->{filtered} = $reason;
}

=item filtered

Returns the reason why the read was filtered; undef otherwise.

=cut

sub filtered { return shift->{filtered} }

=item output

Returns a multiline string of the sequence in Fastq format.  Filtered unpaired reads return an empty string,
while filtered paired reads return an empty record, so pairing won't be broken.

=cut

sub output
{
	my $this = shift;
	if ($this->{filtered})
	{
		if (defined($this->{pair}))
		{
			return '@' . $this->id . "\n\n+\n\n";
		} else
		{
			return '';
		}
	} else
	{
		return '@' . $this->id . "\n" . $this->{seq} . "\n+\n" . $this->{qual} . "\n";
	}
}

=item output_fasta

Returns a multiline string of the sequence in Fasta format.  Returns an empty string if record was filtered and is
unpaired, returns an empty record if was filtered (header but no seq).

=cut

sub output_fasta
{
	my $this = shift;
	if ($this->{filtered})
	{
		if (defined($this->{pair}))
		{
			return '>' . $this->id . "\n\n";
		} else
		{
			return '';
		}
	} else
	{
		my $fasta = '>' . $this->id . "\n";
		my $seq   = $this->{seq};
		while (length($seq) > CHARACTERS_PER_LINE)
		{
			$fasta .= substr($seq, 0, CHARACTERS_PER_LINE) . "\n";
			$seq = substr($seq, CHARACTERS_PER_LINE);
		}
		$fasta .= $seq . "\n";
		return $fasta;
	}
}

=item output_qual

Returns a multiline string of the sequence's quality scores in phred format.  Returns empty string if filtered
and unpaired; returns empty record if paired.

=cut

sub output_qual
{
	my $this = shift;
	if ($this->{filtered})
	{
		if (defined($this->{pair}))
		{
			return '>' . $this->id . "\n\n";
		} else
		{
			return '';
		}
	} else
	{
		my @qual   = map { ord($_) - 33 } split(//, $this->{qual});
		my $output = ">" . $this->id . "\n";
		my $i      = CHARACTERS_PER_LINE - 1;
		while (scalar(@qual) > CHARACTERS_PER_LINE)
		{
			$output .= join(' ', @qual[ 0 .. $i ]) . "\n";
			@qual = @qual[ CHARACTERS_PER_LINE .. $#qual ];
		}
		$output .= join(' ', @qual) . "\n";
		return $output;
	}
}

###############################################################################
## QC METHODS

=back

=head2 QC Methods

=over 5

=item qc $winsize $meanq $minlen $maxn

To perform quality end trimming, minimum length filtering, and filtering reads with too many Ns.

=cut

sub qc
{
	my ($this, $qc_params) = @_;
	$this->trim3            if $qc_params->{trim3};
	$this->trim_terminal_Ns if $qc_params->{trimN};
	$this->qual_end_trim($qc_params->{winsize}, $qc_params->{meanq}) if defined($qc_params->{winsize}) and defined($qc_params->{meanq});
	$this->length_filter($qc_params->{minlen})                 if defined($qc_params->{minlen});
	$this->N_filter($qc_params->{maxn})                        if defined($qc_params->{maxn});
	$this->low_complexity_filter($qc_params->{low_complexity}) if defined($qc_params->{low_complexity});
	$this->low_qual_filter($qc_params->{low_q}, $qc_params->{max_num_low_q}) if defined($qc_params->{low_q}) and defined($qc_params->{max_num_low_q});
	$this->mean_qual_filter($qc_params->{min_mean_q}) if defined($this->{min_mean_q});
}

=item trim3

Trim Q2 bases from 3' end, which correspond to Illumina "Read Segment Quality Control Indicator"

=cut

sub trim3
{
	my ($this) = @_;
	return unless $this->{qual} =~ /(#+)$/;
	my $len = length($this->{qual}) - length($1);
	$this->{qual} = substr($this->{qual}, 0, $len);
	$this->{seq}  = substr($this->{seq},  0, $len);
	$this->del('poor_quality') unless $len > 0;
}

=item trim_terminal_Ns

Removes uninformative Ns at the sequence ends.

=cut

sub trim_terminal_Ns
{
	my $this = shift;
	my $len  = length($this->{seq});
	return unless $len > 0;
	if ($this->{seq} =~ /^(N+)(.*)$/)
	{
		$this->{seq} = $2;
		$this->{qual} = substr($this->{qual}, length($1));
	}
	if ($this->{seq} =~ /(N+)$/)
	{
		$len = length($this->{seq}) - length($1);
		$this->{seq}  = substr($this->{seq},  0, $len);
		$this->{qual} = substr($this->{qual}, 0, $len);
	}
	$len = length($this->{seq});
	$this->del('poor_quality') unless $len > 0;
}

=item qual_end_trim

To trim both ends of the read using a sliding window until the mean quality score surpasses the threshold.
The sequence can be trimmed away to nothing so using the length filter is recommended.

=cut

sub qual_end_trim
{
	my ($this, $winsize, $meanq) = @_;
	$winsize = DEFAULT_WINSIZE unless $winsize;
	$meanq   = DEFAULT_MEANQ   unless $meanq;
	return unless $this->{seq};

	my @qual = map { ord($_) - 33 } split(//, $this->{qual});
	if (@qual < $winsize)
	{
		# FAIL
		$this->del('too_short');
		return;
	}

	# SLIDING WINDOW TRIM LEFT
	my $w      = $winsize - 1;
	my $winsum = 0;
	foreach (@qual[ 0 .. $w ]) { $winsum += $_ }
	while (@qual > $winsize and ($winsum / $winsize) < $meanq)
	{
		$this->{seq}  = substr($this->{seq},  1);
		$this->{qual} = substr($this->{qual}, 1);
		my $q = shift @qual;
		$winsum -= $q;
		$q = $qual[$w];
		$winsum += $q;
	}
	if (($winsum / $winsize) < $meanq)
	{
		# FAIL
		$this->del('too_short');
		return;
	}

	# SLIDING WINDOW TRIM RIGHT
	$w      = @qual - $winsize;
	$winsum = 0;
	foreach (@qual[ $w .. $#qual ]) { $winsum += $_ }
	while (@qual > $winsize and ($winsum / $winsize) < $meanq)
	{
		$this->{seq}  = substr($this->{seq},  0, length($this->{seq}) - 1);
		$this->{qual} = substr($this->{qual}, 0, length($this->{qual}) - 1);
		my $q = pop @qual;
		$winsum -= $q;
		$w = @qual - $winsize;
		$q = $qual[$w];
		$winsum += $q;
	}
	if (($winsum / $winsize) < $meanq)
	{
		# FAIL
		$this->del('too_short');
		return;
	}
	if (length($this->{seq}) != length($this->{qual}))
	{
		warn("End-trim error\n");
		$this->del('error');
	}
}

=item length_filter

If the sequence is shorter than the minimum length, the record is filtered.

=cut

sub length_filter
{
	my ($this, $minlen) = @_;
	$minlen = DEFAULT_MINLEN unless defined($minlen);
	$this->del('too_short') if length($this->{seq}) < $minlen;
}

=item N_filter

If the sequence contains more than the allowed number of Ns, the record is filtered.

=cut

sub N_filter
{
	my ($this, $maxn) = @_;
	$maxn = DEFAULT_MAXN unless defined($maxn);
	my $n = $this->{seq} =~ s/N/N/gi;
	$this->del('too_many_N') if $n > $maxn;
}

=item low_complexity_filter $pct_len

If the sequence is >= $pct_len mono- or di-nucleotide repeats, the record is filtered.

=cut

sub low_complexity_filter
{
	my ($this, $pct_len) = @_;
	$pct_len = DEFAULT_LOW_COMPLEXITY unless defined($pct_len);
	my $seq = $this->{seq};
	my $len = length($seq);
	return unless $len;
	foreach my $nn (qw/AA TT CC GG CA GT CT GA AT CG/)
	{
		my $n = $seq =~ s/$nn/$nn/gi;
		if ($n * 2 / $len >= $pct_len)
		{
			$this->del('low_complexity');
			return;
		}
	}
}

=item low_qual_filter

Filter record if there are too many low-quality bases.

=cut

sub low_qual_filter
{
	my ($this, $minq, $max) = @_;
	my $n = 0;
	foreach my $q (@{$this->qual_arrayref}) { ++$n if $q < $minq }
	$this->del('too_many_low_qual') if $n > $max;
}

=item mean_qual_filter

Filter record if mean quality score is below threshold.

=cut

sub mean_qual_filter
{
	my ($this, $minq) = @_;
	my $meanq = sum(@{$this->qual_arrayref}) / $this->len;
	$this->del('low_mean_qual') if $meanq < $minq;
}

=back

=head1 BUGS AND LIMITATIONS

Reads must be named in accordance with Illumina naming conventions.

=head1 COPYRIGHT

Copyright (c) 2010 U.S. Department of Energy Joint Genome Institute

All right reserved. This program is free software; you can redistribute it
and/or modify it under the same terms as Perl itself.

=head1 AUTHOR

Edward Kirton <ESKirton@LBL.gov>

=cut

###############################################################################
## FASTQDB OBJECT

=head1 NAME

FastqDb - Simple object for handling Fastq files

=head1 SYNOPSIS

    # simple script to QC reads (outputs to stdout)
    use FastqDb;
    my $file=shift or die("Infile required\n");
    my $fastqFile = new FastqFile($file);
    while (my $fastq=$fastqFile->next_seq) {
        $fastq->qc; # see Fastq package for details
        print $fastq->output; # no output if read failed QC
    }

=head1 DESCRIPTION

This module provides an iterator for Fastq files, plus basic QC functions

=head1 METHODS

=over 5

=cut

package Iterator::FastqDb;

use warnings;
use strict;
use IO::File;
use File::Basename;
use Cwd;
use constant {SAMPLE_SIZE => 10_000};

=item new $file

Constructor accepts a path to an existing Fastq file; if undefined, empty, or 'STDIN', 
it will instead read from standard input.

The quality format will be automatically detected using the first sequences.

The output always uses sanger-scaled quality scores.

Optional parameters are:

    format => illumina|sanger; convert if illumina
    barcodes_file => path or 'AUTO' to learn barcode sequences from subset of data
    barcodes_col => column (starting with 1) containing barcode sequence (in file above)

=cut

sub new
{
	my ($class, $file, $params) = @_;

	# OPEN INPUT
	my $fh = new IO::File;
	if (!defined($file) or $file eq '' or $file =~ /^STDIN$/i)
	{
		$file = 'STDIN';
		$fh   = *STDIN;
	} else
	{
		# USE FULL PATH
		my $cwd = getcwd . '/';
		my ($base, $dir, $suffix) = fileparse($file, qr/\.[^.]*/);
		$dir = $cwd if $dir eq './';
		$file = "$dir$base$suffix";
		# OPEN FILE
		if ($file =~ /\.gz$/)
		{
			open($fh, "gunzip -c $file|") or die("Unable to open Fastq file, $file: $!\n");
		} elsif ($file =~ /\.bz2$/)
		{
			open($fh, "bzcat $file|") or die("Unable to open Fastq file, $file: $!\n");
		} else
		{
			open($fh, "<$file") or die("Unable to open Fastq file, $file: $!\n");
		}
	}

	# INIT OBJECT
	my $this = {
		# INPUT:
		file   => $file,
		fhr    => \$fh,
		line   => 0,
		paired => undef,
		# SUMMARY STATS:
		counter => {pass => 0}    # "pass" or reason failed => count
	};
	bless $this, $class;

	# INIT QC PARAMS
	$this->_qc_params($params);

	# SAMPLE FILE FOR TRAINING
	$this->_load_buffer;

	# DETERMINE QUAL
	$this->_check_qual;

	# MOLECULAR BARCODES
	$this->_init_barcodes;

	# DONE
	$this->_process_buffer;

	return $this;
}

# DESTRUCTOR - closes file before going away
sub DESTROY
{
	my $this = shift;
	return unless exists($this->{fhr}) and defined($this->{fhr});
	my $fh = ${$this->{fhr}};
	close($fh) unless $this->{file} eq 'STDIN';
}

##############################################
## INITIALIZE QC PARAMETERS: PRIVATE METHOD ##
##############################################

sub _qc_params
{
	my ($this, $params) = @_;

	# INIT
	$this->{qc_params} = {
		'trim3'          => 0,
		'trimN'          => 0,
		'qual_to_sanger' => undef,
		'roche_mid_len'  => undef,
		'meanq'          => undef,
		'winsize'        => undef,
		'minlen'         => undef,
		'maxn'           => undef,
		'low'            => undef,
		'low_q'          => undef,
		'max_num_low_q'  => undef,
		'min_mean_q'     => undef};
	$this->{barcodes_file}      = undef;
	$this->{barcodes_seq_col}   = 0;
	$this->{barcodes_label_col} = undef;
	return unless defined($params);

	# PARSE USER-SUPPLIED OPTIONS
	foreach my $param (keys %$params)
	{
		my $val = $params->{$param};
		if ($param eq 'roche_mid_len')
		{
			die("Invalid roche_mid_len parameter, $val\n") unless $val > 0;
			$this->{qc_params}->{roche_mid_len} = $val;
		} elsif ($param eq 'barcodes_file')
		{
			$this->{barcodes_file} = $val =~ /^auto$/i ? 'AUTO' : $val;
		} elsif ($param eq 'barcodes_seq_col')
		{
			die("Invalid barcodes_seq_col, $val\n") unless $val >= 1;
			$this->{barcodes_seq_col} = $val - 1;    # change to 0-based index
		} elsif ($param eq 'barcodes_label_col')
		{
			die("Invalid barcodes_label_col, $val\n") unless $val >= 1;
			$this->{barcodes_label_col} = $val - 1;    # change to 0-based index
		} elsif ($param eq 'winsize')
		{
			die("Invalid winsize, $val\n") unless $val >= 1;
			$this->{qc_params}->{'winsize'} = $val;
		} elsif ($param eq 'meanq')
		{
			die("Invalid meanq, $val\n") unless $val >= 1;
			$this->{qc_params}->{'meanq'} = $val;
		} elsif ($param eq 'minlen')
		{
			die("Invalid minlen, $val\n") unless $val >= 0;
			$this->{qc_params}->{'minlen'} = $val;
		} elsif ($param eq 'maxn')
		{
			die("Invalid maxn, $val\n") unless $val >= 0;
			$this->{qc_params}->{'maxn'} = $val;
		} elsif ($param eq 'low')
		{
			die("Invalid low, $val\n") unless $val >= 0 and $val <= 1;
			$this->{qc_params}->{'low'} = $val;
		} elsif ($param eq 'trim3')
		{
			die("Invalid trim3, $val\n") unless $val == 0 or $val == 1;
			$this->{qc_params}->{'trim3'} = $val;
		} elsif ($param eq 'trimN')
		{
			die("Invalid trimN, $val\n") unless $val == 0 or $val == 1;
			$this->{qc_params}->{'trimN'} = $val;
		} elsif ($param eq 'low_q')
		{
			die("Invalid low_q, $val\n") unless $val >= 0 and $val <= 41;
			$this->{qc_params}->{'low_q'} = $val;
		} elsif ($param eq 'max_num_low_q')
		{
			die("Invalid max_num_low_q, $val\n") unless $val >= 0;
			$this->{qc_params}->{'max_num_low_q'} = $val;
		} elsif ($param eq 'min_mean_q')
		{
			die("Invalid min_mean_q, $val\n") unless $val >= 0 and $val <= 41;
			$this->{qc_params}->{'min_mean_q'} = $val;
		} else
		{
			$val = '' unless defined($val);
			die("Invalid parameter, $param => $val\n");
		}
	}
}

##############
## ACCESSORS #
##############

=item file ($new_file)

Returns the Fastq db complete path or 'STDIN'.

=cut

sub file { return shift->{file} }

=item paired

Returns true if paired, false otherwise.

=cut

sub paired { return shift->{paired} }

##############
## ITERATOR ##
##############

=item next_seq

An iterator for the fastq database; returns a single read object.

=cut

sub next_seq
{
	my $this = shift;

	# GET NEXT SEQ
	my $seq = $this->_next_from_buffer;
	$seq = $this->_next_from_file unless defined($seq);
	return undef unless defined($seq);

	# UPDATE COUNTER
	if (my $reason = $seq->filtered)
	{
		if (exists($this->{counter}->{$reason}))
		{
			$this->{counter}->{$reason} += 1;
		} else
		{
			$this->{counter}->{$reason} = 1;
		}
	} else
	{
		$this->{counter}->{pass} += 1;

		# ONLY COUNT BARCODES FOR GOOD SEQS
		$this->{barcodes}->{$seq->barcode} += 1 if defined($this->{barcodes});
	}
	return $seq;

	##################
	# PRIVATE METHODS:

	sub _next_from_buffer
	{
		my $this = shift;
		return undef unless exists($this->{buffer});
		my $rec = shift @{$this->{buffer}};
		delete($this->{buffer}) unless defined($rec);
		return $rec;
	}

	sub _next_from_file
	{
		my $this = shift;
		return undef unless exists($this->{fhr}) and defined($this->{fhr});
		my $fh = ${$this->{fhr}};
		my ($hdr, $seq, $sep, $qual);
		my $error = 0;
		while (<$fh>)
		{
			chomp;
			$this->{line} += 1;
			if (!defined($hdr))
			{
				if (/^@/)
				{
					$hdr   = $_;
					$error = 0;
				} elsif (!$error)
				{
					print "[Line " . $this->{line} . "] Expected header, got:\n$_\n";
				}
				# else: after an error (below), don't complain about every line before the next heaer
			} elsif (!defined($seq))
			{
				if ($_ eq '' or /^[ATCGNKMRYSWBVHDX-]+$/i)#Support IUPAC alphated (ambiguous letters...)
				{
					$seq = $_;
				} else
				{
					print "[Line " . $this->{line} . "] Expected seq, got:\n$hdr\n$_\n";
					$hdr   = undef;
					$error = 1;
				}
			} elsif (!defined($sep))
			{
				if (/^\+/)
				{
					$sep = $_;
				} else
				{
					print "[Line " . $this->{line} . "] Expected '+', got:\n$hdr\n$seq\n$_\n";
					$hdr = $seq = undef;
					$error = 1;
				}
			} elsif (!defined($qual))
			{
				if (length($_) == length($seq))
				{
					$qual = $_;
					last;
				} else
				{
					print "[Line " . $this->{line} . "] Qual doesn't match seq, got:\n$hdr\n$seq\n$sep\n$_\n";
					$hdr = $seq = $sep = undef;
					$error = 1;
				}
			}
		}
		unless (defined($hdr))
		{
			$this->close_file;
			return undef;
		}
		return new Iterator::Fastq($hdr, $seq, $qual, $this->{qc_params}, $this->{barcodes}, $this->{barcode_variants});
	}
}

=item summary

Return summary of good and filtered reads.

=cut

sub summary
{
	my $this    = shift;
	my %counter = %{$this->{counter}};
	# DEEP COPY
	my %summary = %counter;
	return \%summary;
}

=item close_file

This will close the file and delete the file handle.  Normally this is called by the iterator upon reaching the end
of the file, so this method only needs to be called if you wish to close the file prematurely.

=cut

sub close_file
{
	my $this = shift;
	if (exists($this->{fhr}) and defined($this->{fhr}) and $this->{file} ne 'STDIN')
	{
		my $fh = ${$this->{fhr}};
		close($fh);
	}
	delete $this->{fhr};
}

##############################
## TRAINING: PUBLIC METHODS ##
##############################

# WRITE THE BUFFER TO A FILE FOR PARAMETER TRAINING BY AN EXTERNAL PROGRAM.
# THE BUFFER WOULD HAVE BEEN PREPROCESSED (QUAL CONVERSION, TRIM/FILTER, ETC.)

=item write_sample C<$file> C<$fasta>

Write the first sequences to C<$file>, for parameter training.
If C<$file> is undefined, empty, or "STDOUT", the records are written to standard output.
This method does not disrupt pipes when the object is reading from STDIN (or writing to STDOUT) as the sequences were buffered.
The records written will be processed with the same parameters as the remainder of the file.
Records are always in FastqSanger format unless the C<$fasta> flag is true, which generates a Fasta file.
Returns number of records written, mean read length, standard deviation read length, mean pair length, standard deviation of
pair length (pair stats are undefined if unpaired).

=cut

sub write_sample
{
	my ($this, $file, $fasta) = @_;
	my $out = new IO::File;
	if (defined($file) and $file and $file !~ /^STDOUT$/i)
	{
		open($out, ">$file") or die("Unable to open outfile, $file: $!\n");
	} else
	{
		$file = 'STDOUT';
		$out  = *STDOUT;
	}
	foreach my $seq (@{$this->{buffer}})
	{
		print $out $fasta ? $seq->output_fasta : $seq->output;
	}
	close($out) unless $file eq 'STDOUT';
	return $this->{sample_size};
}

=item next_sample_seq

Iterator for sample data, returns seq object.

=cut

sub next_sample_seq
{
	my $this = shift;
	if (!exists($this->{sample_index}))
	{
		$this->{sample_index} = 0;
		return $this->{buffer}->[0];
	} elsif (++$this->{sample_index} == $this->{sample_size})
	{
		delete($this->{sample_index});
		return undef;
	} else
	{
		return $this->{buffer}->[ $this->{sample_index} ];
	}
}

=item sample_array

Returns an array of all sampled sequence objects.  Iterating (above) saves RAM and is the preferred method.

=cut

sub sample
{
	my $this = shift;
	return @{$this->{buffer}};    # deep copy
}

###############################
## TRAINING: PRIVATE METHODS ##
###############################
#
# The begining of the file is buffered and used by one or more methods to train the parameters
# to be used on the entire database.
#
# The following may be learned:
#
# - whether reads are paired or not (output empty reads if paired)
# - whether input is in older FastqIllumina format (and must be converted to Sanger-scale)
# - the molecular barcodes used (if barcoded)
# - which sequencing adapters were used (adapter db must be provided)
# - distributions of read, insert, and fragment sizes (for pairwise assembly)
#
# There are no public methods.

# LOAD BEGINING OF FILE INTO RAM, TO BE USED BY OTHER METHODS FOR TRAINING
sub _load_buffer
{
	my $this = shift;

	# INIT VARS
	my @buffer = ();

	# PRELOAD WITHOUT ANY READ PREPROCESSING
	my $qc_params_user = $this->{qc_params};
	$this->{qc_params} = {
		'trim3'          => 0,
		'trimN'          => 0,
		'qual_to_sanger' => undef,
		'roche_mid_len'  => undef,
		'meanq'          => undef,
		'winsize'        => undef,
		'minlen'         => undef,
		'maxn'           => undef,
		'low'            => undef,
		'low_q'          => undef,
		'max_num_low_q'  => undef,
		'min_mean_q'     => undef};
	my $paired = undef;

	# READ INPUT
	for (my $i = 1; $i <= SAMPLE_SIZE; $i++)
	{
		my $seq = $this->next_seq or last;
		push @buffer, $seq;
		if ($seq->pair)
		{
			if (!defined($paired)) { $paired = 1 }
			elsif ($paired == 0) { 
				#die("Mix of paired/unpaired reads\n") 
			}
		} else
		{
			if (!defined($paired)) { $paired = 0 }
			elsif ($paired == 1) { 
				#die("Mix of paired/unpaired reads\n") 
			}
		}
	}
	$this->{buffer}      = \@buffer;
	$this->{paired}      = $paired;
	$this->{sample_size} = scalar(@buffer) or die("No valid records in " . $this->{file} . "\n");

	# RESTORE QC PARAMETERS TO USER-SPECIFIED VALUES
	$this->{qc_params} = $qc_params_user;
}

# AFTER TRAINING, THE PREVIOUSLY BUFFERED READS MUST BE PROCESSED USING THE TRAINED PARAMETERS
sub _process_buffer
{
	my $this = shift;

	foreach my $seq (@{$this->{buffer}})
	{
		# QUALITY ENCODING
		$seq->convert_qual_to_sanger if $this->{qc_params}->{qual_to_sanger};

		# MOLECULAR BARCODES
		if (exists($this->{barcodes}))
		{
			my $barcode = $seq->barcode;
			if (exists($this->{barcodes}->{$barcode}))
			{
				$this->{barcodes}->{$barcode} += 1;
			} elsif (exists($this->{barcode_variants}->{$barcode}))
			{
				my $corrected_barcode = $this->{barcode_variants}->{$barcode};
				$seq->barcode($corrected_barcode);
				$this->{barcodes}->{$corrected_barcode} += 1;
			} else
			{
				$seq->del('bad_barcode');
			}
		}
		# TRIM/FILTER
		$seq->qc($this->{qc_params});
	}

	# RESET COUNTER (READS WILL BE TALLIED DURING RETRIEVAL)
	$this->{counter} = {};
}

########################################################
## AUTO-DETECT QUALITY SCORE ENCODING: PRIVATE METHOD ##
########################################################

sub _check_qual
{
	my $this = shift;

	# INIT VARS
	my $appears_sanger   = 0;
	my $appears_illumina = 0;

	# COUNT QUAL CHARACTERS IN NON-OVERLAPPING ASCII RANGES
	foreach my $seq (@{$this->{buffer}})
	{
		my @qual = split(//, $seq->qual);
		foreach my $q (@qual)
		{
			my $ascii = ord($q);
			if    ($ascii < 66) { ++$appears_sanger }
			elsif ($ascii > 74) { ++$appears_illumina }
		}
	}

	# SAVE RESULT
	#print "SANGER:".$appears_sanger."\n";
	#print "ILLUMINA:".$appears_illumina."\n";
	$this->{qc_params}->{qual_to_sanger} = $appears_sanger > $appears_illumina ? 0 : 1;
}

########################################
## MOLECULAR BARCODES: PUBLIC METHODS ##
########################################

=item barcode_counts

Returns a hash of barcodes and counts.  Returns undef if reads were not barcoded.
Normally this is only run after iterating through all reads in the database.

=cut

sub barcode_counts
{
	my $this = shift;
	return undef unless defined($this->{barcodes});
	my %barcodes = %{$this->{barcodes}};    # deep copy
	return \%barcodes;
}

=item barcode_labels

Returns a hash of barcodes and labels.  Returns undef if reads were not barcoded or labels not defined.

=cut

sub barcode_labels
{
	my $this = shift;
	return undef unless defined($this->{barcode_labels});
	my %barcode_labels = %{$this->{barcode_labels}};    # deep copy
	return \%barcode_labels;
}

########################################
## MOLECULAR BARCODES: PRIVATE METHOD ##
########################################

sub _init_barcodes
{
	my $this = shift;
	return unless defined($this->{barcodes_file});

	if ($this->{barcodes_file} eq 'AUTO')
	{
		$this->_learn_barcodes;
	} else
	{
		$this->_load_barcodes;
	}
}

# load the barcodes provided in a tabular or Fasta file
sub _load_barcodes
{
	my $this = shift;
	my $file = $this->{barcodes_file} or die("No barcodes file provided\n");

	# INIT VARS
	my $seq_col   = $this->{barcodes_seq_col};
	my $label_col = $this->{barcodes_label_col};
	my $barcode;
	my %barcodes         = ();
	my %barcode_variants = ();
	my %barcode_labels   = ();

	# READ FILE
	open(IN, "<$file") or die("Unable to open molecular barcodes file, $file: $!\n");
	my @data = <IN>;
	close(IN);

	# DETERMINE FORMAT AND PARSE ACCORDINGLY
	if ($data[0] =~ /^>\w+/)
	{
		# FASTA
		my $label;
		while (@data)
		{
			my $line = shift @data;
			chomp $line;
			if ($line =~ /^>(.+)$/)
			{
				$label = $1;
			} elsif ($line =~ /^[ATCG]+$/i)
			{
				$barcode .= uc($line);
				$barcodes{$barcode} = 0;
				my @variants = _generate_barcode_variants($barcode);
				foreach my $variant (@variants)
				{
					$barcode_variants{$variant} = $barcode;
				}
				$barcode_labels{$barcode} = $label;
			} elsif ($line)
			{
				die("Error reading molecular barcodes file\n");
			}
			# else blank line
		}
	} else
	{
		# TABULAR
		while (@data)
		{
			my $line = shift @data;
			chomp $line;
			my @row = split(/\t/, $line);
			$barcode = uc($row[$seq_col]);
			die("Invalid barcode, $barcode\n") unless $barcode =~ /^[ATCGN]+$/;
			$barcodes{$barcode} = 0;
			my @variants = _generate_barcode_variants($barcode);
			foreach my $variant (@variants)
			{
				$barcode_variants{$variant} = $barcode;
			}
			if (defined($label_col))
			{
				$barcode_labels{$barcode} = $row[$label_col];
			}
		}
	}

	# SAVE
	$this->{barcodes}         = \%barcodes;
	$this->{barcode_variants} = \%barcode_variants;
	$this->{barcode_labels}   = \%barcode_labels;
}

# DISCOVER BARCODES USING ABUNDANT SEQUENCES AND ASSUMING ALL REAL BARCODES HAVE COUNTS
# WITHIN THREE STANDARD DEVIATIONS OF THE MEAN ABUNDANCE

sub _learn_barcodes
{
	my $this = shift;

	# INIT VARS
	my %putative_barcodes = ();    # seq => count
	my %barcodes          = ();    # sequence of real barcode or variant => number of times observed
	my %barcode_variants  = ();    # variant barcode sequence => real barcode sequence
	my $barcode;
	my $n;
	my %counts = ();

	# COUNT ALL BARCODE SEQUENCES IN SAMPLE
	foreach my $seq (@{$this->{buffer}})
	{
		if ($barcode = $seq->barcode)
		{
			if (exists($putative_barcodes{$barcode}))
			{
				$putative_barcodes{$barcode} += 1;
			} else
			{
				$putative_barcodes{$barcode} = 1;
			}
		}
	}

	# SORT PUTATIVE BARCODES IN DESCENDING ORDER OF COUNTS
	foreach $barcode (keys %putative_barcodes)
	{
		$n = $putative_barcodes{$barcode};
		if (exists($counts{$n}))
		{
			push @{$counts{$n}}, $barcode;
		} else
		{
			$counts{$n} = [$barcode];
		}
	}
	my @counts = sort { $b <=> $a } keys %counts;

	# PICK REAL BARCODES (ASSUMES ABUNDANCE OF ALL REAL BARCODES IS WITHIN 3SD OF MEAN)
	my @real_barcode_counts = ();
	my $mean;
	my $stddev;
	while ($n = shift @counts)
	{
		if (@real_barcode_counts == 0)
		{
			1;
		} elsif (@real_barcode_counts < 3)
		{
			($mean, $stddev) = _mean_stddev(\@real_barcode_counts);
			last if $n * 100 < $mean;
		} else
		{
			($mean, $stddev) = _mean_stddev(\@real_barcode_counts);
			last if $n < ($mean - 3 * $stddev);
		}

		# REAL BARCODE(S)
		foreach $barcode (@{$counts{$n}})
		{
			next unless exists($putative_barcodes{$barcode});
			push @real_barcode_counts, $n;
			$barcodes{$barcode} = 0;
			# DELETE VARIANTS
			my @variants = _generate_barcode_variants($barcode);
			foreach my $variant (@variants)
			{
				if (exists($putative_barcodes{$variant}))
				{
					delete($putative_barcodes{$variant});
				}
				$barcode_variants{$variant} = $barcode;
			}
		}
		delete($counts{$n});
	}

	# GENERATE LABELS HASH (LABELS ARE JUST THE SEQ IN THIS CASE)
	my %labels = ();
	foreach $barcode (keys %barcodes) { $labels{$barcode} = $barcode }

	# SAVE RESULTS
	$this->{barcodes}         = \%barcodes;
	$this->{barcode_variants} = \%barcode_variants;
	$this->{barcode_labels}   = \%labels;
	return;
}

############################
# NONMEMBER HELPER FUNCTION

# calculate mean and standard deviation
sub _mean_stddev
{
	my $ar  = shift;
	my $n   = scalar(@$ar);
	my $tot = 0;
	foreach my $x (@$ar) { $tot += $x }
	my $mean = $tot / $n;
	my $d    = 0;
	foreach my $x (@$ar) { $d += ($mean - $x)**2 }
	my $stddev = sqrt($d / $n);
	return ($mean, $stddev);
}

# generate all 1-base mismatch sequences of molecular barcode provided
sub _generate_barcode_variants
{
	my $seq      = shift;
	my @seq      = split(//, $seq);
	my @variants = ();
	for (my $i = 0; $i <= $#seq; $i++)
	{
		foreach my $nt (qw/A T C G N/)
		{
			next if $nt eq $seq[$i];
			my @variant = @seq;
			$variant[$i] = $nt;
			push @variants, join('', @variant);
		}
	}
	return @variants;
}

1;

=back

=head1 BUGS AND LIMITATIONS

Uses a modified Illumina naming convention: base#barcode/pair

=head1 COPYRIGHT

Copyright (c) 2010 U.S. Department of Energy Joint Genome Institute

All right reserved.

=head1 AUTHOR

Edward Kirton <ESKirton@LBL.gov>
Julien Tremblay <julien.tremblay@mail.mcgill.ca>

=cut

__END__
