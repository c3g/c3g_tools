#!/usr/bin/env python

### Edouard Henrion (2021/02/25)
### Convert headers of raw/undemultiplexed fastq files from MGI to Illumina format

import os
import sys
import re
import getopt
import itertools
from functools import partial
from multiprocessing import Pool
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO, bgzf
import gzip

instrument_ID = ""
run_number = 0
index1_length = 0
index2_length = 0

def getarg(argument):
    rlist=""
    fastq1=""
    fastq2="."
    gz=0
    op="."
    instID="."
    run=0
    i1_len=10
    i2_len=0
    ncpus=4
    optli,arg = getopt.getopt(argument[1:],"i:m:o:g:x:y:s:r:t:h",['input','mate','output_prefix','gzip','I1_length','I2_length','instrument','run_number','threads','help'])
    if len(optli) < 1 :
        usage()
        sys.exit("Error : Missing argument(s)")
    for option, value in optli:
        if option in ("-i","--input"):
            fastq1=str(value)
        if option in ("-m","--mate"):
            fastq2=str(value)
        if option in ("-g","--gzip"):
            gz=int(value)
        if option in ("-o","--output_prefix"):
            op=str(value)
        if option in ("-x","--I1_length"):
            i1_len=int(value)
        if option in ("-y","--I2_length"):
            i2_len=int(value)
        if option in ("-s","--instrument"):
            instID=str(value)
        if option in ("-r","--run_number"):
            run=int(value)
        if option in ("-t","--threads"):
            ncpus=int(value)
        if option in ("-h","--help"):
            usage()
            sys.exit()
    if not os.path.exists(fastq1) :
        sys.exit("Error - fastq file not found:\n"+fastq1)
    if not os.path.exists(fastq2) :
        sys.exit("Error - fastq file not found:\n"+fastq2)
    if instID == "." :
        sys.exit("Error - instrument ID is mandatory...\n")
    if run == 0 :
        sys.exit("Error - Run number is mandatory...\n")
    if op == "." :
        output_f1 = os.path.splitext(fastq1)[0]
        output_f2 = os.path.splitext(fastq2)[0]
        if gz == 0 :
            output_f1 = os.path.splitext(output_f1)[0]
            output_f2 = os.path.splitext(output_f2)[0]
        output_f1 = output_f1 + ".converted.fastq.gz"
        output_f2 = output_f2 + ".converted.fastq.gz"
    else :
        output_f1 = op + "_R1.fastq.gz"
        output_f2 = op + "_R2.fastq.gz"

    global instrument_ID
    global run_number
    global index1_length
    global index2_length
    instrument_ID = instID
    run_number = run
    index1_length = i1_len
    index2_length = i2_len

#    global f1_out
#    global f2_out

    return gz, fastq1, fastq2, output_f1, output_f2, ncpus

def usage():
    print "\n------------------------------------------------------------------------"
    print "convertFastqHeader.py will convert fastq files with MGI-formated headers"
    print " into fastq files with Illumina-formated headers."
    print "Since the barcode sequence is part of the Illumina header and because"
    print " MGI places the barcode sequence(s) at the end of the 2nd read,"
    print " providing the read2 fastq file is mandatory to run this script."
    print "The instrument ID and the run number have to be passed in parameter"
    print " since they are needed for the Illumina header but not provided by the"
    print " MGI header."
    print "Size of the barcode sequences can be customized especially to handle"
    print " both single-index and dual-index datasets."
    print "Converted outputs are always written in compressed fastq.qz format."
    print "This program was written by Edouard HENRION."
    print "For more information, contact: edouard.henrion@computationalgenomics.ca"
    print "------------------------------------------------------------------------\n"
    print "USAGE : fastqIcounter.py [option] "
    print "       -i :        read1 input fastq file"
    print "       -m :        mate (read2) input fastq file"
    print "       -g :        gziped input yes: 0 ; no: 1 (Default 0)"
    print "       -o :        output file prefix (default: output in the same directory as input, adding '.converted' suffix at the input base name)"
    print "       -x :        size of the first barcode (default 10)"
    print "       -y :        size of the second barcode (default 0 i.e. single-end)"
    print "       -r :        run number"
    print "       -s :        instrument/sequencer ID"
    print "       -t :        number of threads/cpus for parallelized processing (default 4)"
    print "       -h :        this help \n"

# https://stackoverflow.com/a/43922107/6198494
def chunker_list(seq, size):
    print "Start chunking the input fastq files..."
    return (seq[i::size] for i in range(size))

def parseMGI(
    header_string,
    instrument_ID,
    run_number,
    sample_barcode
    ):

    m = re.search("(?P<fcid>V\w+)L(?P<lane>\d{1})C(?P<col>\d{3})R(?P<row>\d{3})(?P<tile>\d+)/(?P<read>\d{1})", header_string)
    if m:
        # Pack into dictionary
        header_dict = {
            "instrument": instrument_ID,
            "run_number": run_number,
            "sample_barcode": sample_barcode,
            "lane": int(m.group('lane')),
            "x-pos": int(m.group('col')),
            "y-pos": int(m.group('row')),
            "tile": int(m.group('tile')),
            "flowcell-id": m.group('fcid'),
            "read_no": int(m.group('read'))
        }
    else:
        sys.exit("Could not parse MGI header properly : " + header_string)
    return header_dict

def convertIllumina(
    header_dict
    ):

    template_str = "{instrument}:{run_number}:{flowcell_id}:{lane}:{tile}:{x_pos}:{y_pos} {read}:N:0:{sample_barcode}"
    formatted_string = template_str.format(
        instrument = header_dict["instrument"],
        run_number = header_dict["run_number"],
        flowcell_id = header_dict["flowcell-id"],
        lane = header_dict["lane"],
        tile = header_dict["tile"],
        x_pos = header_dict["x-pos"],
        y_pos = header_dict["y-pos"],
        read = header_dict["read_no"],
        sample_barcode = header_dict["sample_barcode"]
    )
    return formatted_string

def processMGIfastq(
    data,
    ):
    mgi_head1, seq1, qual1 = data[0]
    mgi_head2, seq2, qual2 = data[1]

    if index2_length == 0:
        # Single index
        sample_barcode = str(seq2)[-index1_length:]
    else:
        # Dual Index
        sample_barcode = str(seq2)[-(index2_length+index1_length):-index1_length]+"+"+str(seq2)[-index1_length:]
#    print sample_barcode
    illumina_head1 = convertIllumina(
        parseMGI(
            mgi_head1,
            instrument_ID,
            run_number,
            sample_barcode
        )
    )

    illumina_head2 = convertIllumina(
        parseMGI(
            mgi_head2,
            instrument_ID,
            run_number,
            sample_barcode
        )
    )

    # Parallelization attempt
#    f1_out.write(illumina_records[0])
#    f2_out.write(illumina_records[1])
#    return sample_barcode
    return ["@%s\n%s\n+\n%s\n" % (illumina_head1, seq1, qual1), "@%s\n%s\n+\n%s\n" % (illumina_head2, seq2, qual2)]


def main():
    gz, fastq1, fastq2, output_fastq1, output_fastq2, ncpus = getarg(sys.argv)

    # Parallelization attempt
#    # create a pool of processing nodes
#    p = Pool(ncpus)
#
#    # use partial to create a function needing only one argument
#    fastqfunc = partial(
#        processMGIfastq
#    )

    if gz == 0 :
        open_input = compile("gzip.open", "<string>", "eval")
    else :
        open_input = compile("open", "<string>", "eval")

    with eval(open_input)(fastq1, 'rb') as f1_in, eval(open_input)(fastq2, 'rb') as f2_in:
        print "fastqs opened for reading :\n" + "\n".join([fastq1, fastq2])

        with gzip.open(output_fastq1, "wt") as f1_out, gzip.open(output_fastq2, "wt") as f2_out:
            print "fastqs opened for writing :\n" + "\n".join([output_fastq1, output_fastq2])

            # Serial mode
            for mgi_r1_record, mgi_r2_record in itertools.izip(FastqGeneralIterator(f1_in), FastqGeneralIterator(f2_in)):
                illumina_records = processMGIfastq([mgi_r1_record, mgi_r2_record])
                f1_out.write(illumina_records[0])
                f2_out.write(illumina_records[1])

            # Parallelization attempt
#            chunk_count = 0
#            for i in range(10000):
#            for chunked_fastqs in chunker_list(list(itertools.izip(FastqGeneralIterator(f1_in), FastqGeneralIterator(f2_in))), ncpus):
#                chunk_count += 1
#                print "Processing chunk #" + str(chunk_count) + " out of " + str(10000)
#                illumina_fastqs = p.imap(fastqfunc, list(itertools.izip(FastqGeneralIterator(f1_in), FastqGeneralIterator(f2_in)))[i::10000])

#                for illumina_fastq in illumina_fastqs:
#                    f1_out.write(illumina_fastq[0])
#                    f2_out.write(illumina_fastq[1])

    print "Fastq headers converted successfully !"

main()
