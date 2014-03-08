#!/usr/bin/env perl

=head1 NAME

Fasta - Simple object for Fasta sequence

=head1 SYNOPSIS

    my $seq=new Fasta( $hdr, $seq );
    $seq->qc;
    print $seq->output;

=head1 DESCRIPTION

Object for a single read sequence, with methods for basic manipulation and quality control.

=head1 METHODS

=over 5

=cut

package Iterator::Fasta;

use warnings;
use strict;
use constant {
    CHARACTERS_PER_LINE => 80, # for formatting Fasta output
    CLOSE_ENOUGH_TO_END => 6,   # hits of adapters this close to end will result in trimming to be done to end
};

=item new $hdr $seq

Initialize new sequence object.

=cut

sub new {
    my ($class,$hdr,$seq)=@_;
    chomp $hdr;
    chomp $seq;
    die("Incomplete Fasta record: hdr=$hdr, seq=$seq\n") unless $hdr and $seq;
    # INIT
    my $this={
        hdr=>$hdr, # complete header (e.g. ">A/1 description")
        seq=>$seq # complete sequence without newlines
    };
    bless $this,$class;
    $this->_parse_header;
    return $this;
}

sub _parse_header {
    my $this=shift;

    my $hdr=$this->{hdr};
    $this->{id}=undef; # complete ID (e.g. "A/1#GATTACA"); always defined
    $this->{base}=undef; # base ID only (e.g. "A"); always defined
    $this->{pair}=undef; # pair ID only (e.g. "1"); only defined if paired
    $this->{barcode}=undef; # barcode sequence (upper-case); only defined it barcoded

    if ($hdr =~ /^>(\S+)#([aAtTcCgGnN]+)\/([12])/) { # barcoded, paired
        $this->{base}=$1;
        $this->{barcode}=uc($2);
        $this->{pair}=$3;
    } elsif ($hdr =~ /^>(\S+)\/([12])#([aAtTcCgGnN]+)/) { # barcoded, paired
        $this->{base}=$1;
        $this->{pair}=$2;
        $this->{barcode}=uc($3);
    } elsif ($hdr =~ /^>(\S+)\/([12])/) { # paired
        $this->{base}=$1;
        $this->{pair}=$2;
    } elsif ($hdr =~ /^>(\S+)#([aAtTcCgGnN]+)/) { # barcoded, unpaired
        $this->{base}=$1;
        $this->{barcode}=uc($2);
    } elsif ($hdr =~ /^>(\S+)/) { # unpaired or other
        $this->{base}=$1;
    } else {
        die "Unable to parse header: $hdr\n";
    }

    $this->{id}=$this->{base};
    $this->{id} .= "/".$this->{pair} if $this->{pair};
    $this->{id} .= "#".$this->{barcode} if $this->{barcode}
}

=item header ($new_hdr)

Returns the object's header line.  You can also use this method to give the sequence a new header.

=cut

sub header {
    my ($this,$new_hdr)=@_;
    if (defined($new_hdr) and $new_hdr) {
        $new_hdr = $1 if $new_hdr =~ /^[>@](.+)$/;
        $new_hdr = '>'.$new_hdr unless $new_hdr =~ /^>/;
        $this->{hdr} = $new_hdr;
        $this->_parse_header();
    }
    return $this->{hdr};
}

=item id

Returns the object's ID, which is the sequence's unique identifier without comments which may be present in the header.
It cannot be changed directly, but will be updated whenever the header, base, or pair is changed.

=cut

sub id { return shift->{id} }

=item base

If the read is paired, returns it's base ID (same for both members of pair); returns undefined if not paired.
Providing the optional argument will change the base (and the id and header).

=cut

sub base {
    my ($this,$base)=@_;
    if ($base) {
        $this->{id}=$this->{base}=$base;
        $this->{id} .= "/".$this->{pair} if $this->{pair};
        $this->{id} .= "#".$this->{barcode} if $this->{barcode}
    }
    return $this->{base};
}

=item barcode

If the read contains an Illumina molecular barcode ID, it will be returned; otherwise returns undef.
Supplying an optional argument will set the barcode.

=cut

sub barcode {
    my ($this,$barcode)=@_;
    if ($barcode) {
        $this->{barcode}=$barcode;
        $this->{id}=$this->{base};
        $this->{id} .= "/".$this->{pair} if $this->{pair};
        $this->{id} .= "#".$barcode;
        $this->_recreate_header;
    }
    return $this->{barcode};
}

=item pair

Returns the read's ord in the pair; undef otherwise.  It may be changed by supplying the optional extra argument.
To clear the pairing information, pass it an empty string, not NULL.

=cut

sub pair {
    my ($this,$pair)=@_;
    if (defined($pair)) {
        if ($pair) {
            $this->{pair}=$pair;
            $this->{id}=$this->{base};
            $this->{id} .= "/".$pair;
            $this->{id} .= "#".$this->{barcode} if $this->{barcode}
        } else {
            # delete pairing information (e.g. singleton)
            $this->{pair}=undef;
            $this->{id}=$this->{base};
            $this->{id} .= "#".$this->{barcode} if $this->{barcode}
        }
    }
    return $this->{pair};
}

=item seq ($new_seq)

Returns the read's complete sequence, without newlines.  Optional argument changes it.

=cut

sub seq {
    my ($this,$new_seq)=@_;
    if ($new_seq) {
        chomp $new_seq;
        $this->{seq}=$new_seq;
    }
    return $this->{seq};
}

=item revcomp

Reverse-complements a sequence.

=cut

sub revcomp {
    my $this=shift;
    return unless $this->{seq};
    $this->{seq} =~ tr/ATCGatcg/TAGCtagc/;
    my @seq=reverse split(//, $this->{seq});
    $this->{seq}=join('', @seq);
}

=item output

Returns a multiline string of the sequence in Fasta format.  Returns no output if sequence is empty.

=cut

sub output {
    my $this=shift;
    return '' unless $this->{seq}; # will be empty if failed QC
    return $this->{hdr}."\n"._format($this->{seq});
}

# PRIVATE METHOD TO ADD NEWLINE CHARACTERS 
sub _format {
    my $old=shift;
    return '' unless $old;
    return $old unless CHARACTERS_PER_LINE;
    my $new='';
    while (length($old)> CHARACTERS_PER_LINE) {
        $new .= substr($old,0, CHARACTERS_PER_LINE)."\n";
        $old = substr($old, CHARACTERS_PER_LINE);
    }
    $new .= $old."\n";
    return $new;
}

###############################################################################
## QC METHODS

=back

=head2 QC Methods

=over 5

=item qc $winsize $meanq $minlen $maxn

To perform minimum length filtering, and filtering reads with too many Ns.

=cut

sub qc {
    my ($this, $minlen, $maxn)=@_;
    $this->trim_terminal_Ns;
    $this->length_filter($minlen);
    $this->N_filter($maxn);
}

=item trim_terminal_Ns

Discard uninformative Ns from the ends of the sequence.

=cut

sub trim_terminal_Ns {
    my $this=shift;
    if ($this->seq =~ /^(N+)(.*)$/) {
        $this->seq($2);
    }
    if ($this->seq =~ /^(.*)N+$/) {
        $this->seq($1);
    }
}

=item length_filter

If the sequence is shorter than the minimum length, the sequence string is emptied so they will not be
returned by the output method.  Returns true if sequence was filtered.

=cut

sub length_filter {
    my ($this,$minlen)=@_;
    my $len=length($this->{seq});
    if ( !defined($minlen) or $len >= $minlen) { 
        return 0;
    } else {
        $this->{seq}='';
        return 1;
    }
}

=item N_filter

If the sequence contains more than the allowed number of Ns, the sequence string is emptied.
Returns true if sequence was filtered.

=cut

sub N_filter {
    my ($this,$maxn)=@_;
    return 1 unless defined($maxn);
    # COUNT NUMBER OF 'N' BASES IN TRIMMED SEQ
    my $tmpseq=$this->{seq};
    my $n= $tmpseq=~ s/N//gi;
    if ($n <= $maxn) {
        return 0;
    } else {
        # FAIL
        $this->{seq}='';
        return 1;
    }
}

=item low_complexity_filter $pct_len

If the sequence is >= $pct_len mono- or di-nucleotide repeats, clears the sequence string and returns true.

=cut

sub low_complexity_filter {
    my ($this,$pct_len)=@_;
    my $seq=$this->{seq};    
    $seq =~ s/\n//g;
    my $len = length($seq);
    my $filter=0;
    foreach my $nn (qw/AA TT CC GG CA GT CT GA AT CG/) {
        my $n = $seq =~ s/$nn/$nn/g;
        if ($n >= $pct_len/200*$len) {
            $filter = 1;
            last;
        }
    }
    if ($filter) {
        $this->{seq}='';
    }
    return $filter;
}

#=item trim
#
#Given start and end coordinates of adapter/primer sequence, trims sequence string and returns final length.
#
#=cut
#
#sub trim {
#    my ($this, $start, $end)=@_;
#    return unless defined($start) and defined($end);
#    ($start,$end)=($end,$start) if $end<$start; 
#    my $len=length($this->{seq});
#    #print STDERR "- trim read, ", $this->id, " ($len bp) from $start-$end\n"; # DEBUG
#
#    if ($start <= CLOSE_ENOUGH_TO_END and $len - $end <= CLOSE_ENOUGH_TO_END ) {
#        # fail read
#        $this->{seq}='';
#    } elsif ($start <= CLOSE_ENOUGH_TO_END ) {
#        # trim left
#        $this->{seq}=substr($this->{seq},$end);
#    } elsif ($len - $end <= CLOSE_ENOUGH_TO_END ) {
#        # trim right
#        $this->{seq}=substr($this->{seq},0,$start);
#    }
#    return length($this->{seq});
#    #print STDERR "\t+ length after trimming is ", length($this->{seq})," bp\n"; # DEBUG
#}

=back

=head1 BUGS AND LIMITATIONS

No support for paired reads.

=head1 COPYRIGHT

Copyright (c) 2010 U.S. Department of Energy Joint Genome Institute

All right reserved. This program is free software; you can redistribute it
and/or modify it under the same terms as Perl itself.

=head1 AUTHOR

Edward Kirton <ESKirton@LBL.gov>

=cut

###############################################################################
## FASTADB OBJECT

=head1 NAME

FastaDb - Simple object for handling Fasta files

=head1 SYNOPSIS

    # simple script to QC reads (outputs to stdout)
    use FastaDb;
    my $file=shift or die("Infile required\n");
    my $fastaFile = new FastaFile($file);
    while (my $fasta=$fastaFile->next_seq) {
        $fasta->qc; # see Fasta package for details
        print $fasta->output; # no output if read failed QC
    }

=head1 DESCRIPTION

This module provides an iterator for Fasta files, plus basic QC functions

=head1 METHODS

=over 5

=cut

package Iterator::FastaDb;

use warnings;
use strict;
use IO::File;
use File::Basename;
use Cwd;

sub new {
    my ($class,$file)=@_;
    die("Fasta file required\n") unless $file;

    my $fh = new IO::File;
    if ($file =~ /^stdin$/i) {
        $fh = *STDIN;
    } else {
        # USE FULL PATH
        my $cwd=getcwd.'/';
        my ($base,$dir,$suffix)=fileparse($file, qr/\.[^.]*/);
        $dir=$cwd if $dir eq './';
        $file="$dir$base$suffix";
        open($fh, "<$file") or die("Unable to open fasta file, $file: $!\n");
    }
 
    my $this = { file => $file, next_header => undef, fhr => \$fh };
    bless $this,$class;
    return $this;
}

# DESTRUCTOR - closes file before going away
sub DESTROY {
    my $this=shift;
    return unless exists($this->{fhr}) and defined($this->{fhr});
    my $fh=${$this->{fhr}};
    close $fh;
}

=item file ($new_file)

Returns the Fasta db complete path.  You may also use this method to open a new file for reading.

=cut

sub file {
    my ($this,$new_file)=@_;
    if ($new_file) {
        $this->close_file;
        my $fh=new IO::File('<'.$new_file) or die("Unable to open Fasta file, $new_file: $!\n");
        $this->{fhr}=\$fh;
    }
    return $this->{file};
}

=item next_seq

An iterator for the fasta database; returns a single sequence as multiline text.

=cut

sub next_seq {
    my $this=shift;
    return undef unless exists($this->{fhr}) and defined($this->{fhr});

    # INIT VARS
    my $fh = ${$this->{fhr}};
    my $hdr=undef;
    my $seq='';
    if (defined($this->{next_header})) {
        $hdr=$this->{next_header};
        $this->{next_header}=undef;
    } else {
        # GET FIRST HEADER
        while (<$fh>) {
            chomp;
            if (/^>\S+/) {
                $hdr=$_;
                last;
            } elsif (/^#/ or ! $_) {
                next;
            } else {
                die("Invalid Fasta input file; choked on line: $_\n");
            }
        }
    }
    unless (defined($hdr)) {
        $this->close_file;
        return undef;
    }

    # GET SEQ
    while (<$fh>) {
        chomp;
        if (/^>\S+/) {
            $this->{next_header}=$_;
            last;
        } elsif($_) {
            $seq .= $_;
        }
    }

    # CHECK IF END OF FILE OR PART
    if ( !defined($this->{next_header}) ) {
        $this->close_file;
    }

    # RETURN ONE FASTA SEQ
    return new Iterator::Fasta($hdr,$seq);
}

=item close_file

This will close the file and delete the file handle.  Normally this is called by the iterator upon reaching the end
of the file, so this method only needs to be called if you wish to close the file prematurely.

=cut

sub close_file {
    my $this=shift;
    return unless exists($this->{fhr}) and defined($this->{fhr});
    my $fh=${$this->{fhr}};
    close $fh;
    delete $this->{fhr};

}

1;

=back

=head1 BUGS AND LIMITATIONS

No support for paired reads.

=head1 COPYRIGHT

Copyright (c) 2010 U.S. Department of Energy Joint Genome Institute

All right reserved. This program is free software; you can redistribute it
and/or modify it under the same terms as Perl itself.

=head1 AUTHOR

Edward Kirton <ESKirton@LBL.gov>
Julien Tremblay <julien.tremblay@mail.mcgill.ca>

=cut
