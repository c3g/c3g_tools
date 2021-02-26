#!/usr/bin/python

### Edouard Henrion (2021/02/25)
### Convert headers of raw/undemultiplexed fastq files from MGI to Illumina format

import os
import sys
import string
import getopt
import re
import itertools
from Bio import SeqIO, bgzf
import gzip
import operator


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
    optli,arg = getopt.getopt(argument[1:],"i:m:o:g:x:y:s:r:h",['input','mate','output_prefix','gzip','I1_length','I2_lenght','instrument','run_number','help'])
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
        output_f1 = os.path.splitext(fastq1)
        output_f2 = os.path.splitext(fastq2)
        if gz == 0 :
            output_f1 = os.path.splitext(output_f1)
            output_f2 = os.path.splitext(output_f2)
        output_f1 = output_f1 + ".converted.fastq.gz"
        output_f1 = output_f2 + ".converted.fastq.gz"
    else :
        output_f1 = op + "_R1.fastq.gz"
        output_f2 = op + "_R2.fastq.gz"

    return gz, fastq1, fastq2, output_f1, output_f2, instID, run, i1_len, i2_len

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
    print "       -h :        this help \n"

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

def main():
    gz, fastq1, fastq2, output_fastq1, output_fastq2, instrument_ID, run_number, barcode1_length, barcode2_length = getarg(sys.argv)

    if gz == 0 :
        f1 = gzip.open(fastq1, 'rb')
        f2 = gzip.open(fastq2, 'rb')
    else :
        f1=fastq1
        f2=fastq2
    print "fastqs opened for reading :\n" + "\n".join([fastq1, fastq2])

    handle_out_R1 = gzip.open(output_fastq1, "wt")
    handle_out_R2 = gzip.open(output_fastq2, "wt")
    print "fastqs opened for writing :\n" + "\n".join([output_fastq1, output_fastq2])

    # Loop over both R1 and R2 fastq in parallel
    # R2 is needed to fech the barcode sequcence(s)
    # One coupled-record at a time, both R1 and R2 headers have to be converted before passing to the next record
    for mgi_record1, mgi_record2 in itertools.izip(SeqIO.parse(f1, "fastq"), SeqIO.parse(f2, "fastq")):
        # Retrieve the barcode sequence(s) from the R2 fastq
        if barcode2_length == 0:
            # Single index
            sample_barcode = str(mgi_record2.seq)[-barcode1_length:]
        else:
            # Dual Index
            sample_barcode = str(mgi_record2.seq)[-(barcode2_length+barcode1_length):-barcode1_length]+"+"+str(mgi_record2.seq)[-barcode1_length:]

        # Parse the MGI header of R1 e.g. V300057102L1C001R00400000262/1
        header_dict_R1 = parseMGI(
            mgi_record1.id,
            instrument_ID,
            run_number,
            sample_barcode
        )
        # Parse the MGI header of R2 e.g. V300057102L1C001R00400000262/2
        header_dict_R2 = parseMGI(
            mgi_record2.id,
            instrument_ID,
            run_number,
            sample_barcode
        )

        # Create Illumina header for R1
        illumina_header_R1 = convertIllumina(header_dict_R1)
        # Create Illumina header for R2
        illumina_header_R2 = convertIllumina(header_dict_R2)

        # Create the Illumina output record for R1
        illumina_record1 = mgi_record1
        illumina_record1.description = illumina_header_R1
        illumina_record1.id = illumina_header_R1.split(' ')[0]
        illumina_record1.name = illumina_header_R1.split(' ')[0]
        # Create the Illumina output record for R2
        illumina_record2 = mgi_record2
        illumina_record2.description = illumina_header_R2
        illumina_record2.id = illumina_header_R2.split(' ')[0]
        illumina_record2.name = illumina_header_R2.split(' ')[0]

        # write R1 record into output
        handle_out_R1.write(illumina_record1.format('fastq'))
        # write R2 record into output
        handle_out_R2.write(illumina_record2.format('fastq'))
        
    print "Fastq headers converted successfully !"

main()
