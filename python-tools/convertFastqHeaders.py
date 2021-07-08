#!/usr/bin/env python

### Edouard Henrion (2021/02/25)
### Convert headers of fastq file from MGI to Illumina format

# Update - 2021/03/15
### Now takes the path of the folder which contains the fastq files, and process all the fastq files present in the folder
### Input fastq files has to be demultpixed prior conversion i.e. needs to have the sample barcode in header
### Input fastq files have to be MGI-formated otherwise result fastq files will be empty...

import os
import sys
import re
import getopt
from glob import glob
from itertools import izip, islice
from functools import partial
from multiprocessing import Pool
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip

instrument_ID = ""
run_number = 0
index1_length = 0
index2_length = 0

def getarg(argument):
    rlist=""
    input_path=""
    output_path="."
    gz=0
    instID="."
    run=0
    ncpus=4
#    optli,arg = getopt.getopt(argument[1:],"i:o:g:s:r:t:h",['input_path','output_path','gzip','instrument','run_number','threads','help'])
    optli,arg = getopt.getopt(argument[1:],"i:o:g:s:r:h",['input_path','output_path','gzip','instrument','run_number','help'])
    if len(optli) < 1:
        usage()
        sys.exit("Error : Missing argument(s)")
    for option, value in optli:
        if option in ("-i","--input_path"):
            input_path=str(value)
        if option in ("-o","--output_path"):
            output_path=str(value)
        if option in ("-g","--gzip"):
            gz=int(value)
        if option in ("-s","--instrument"):
            instID=str(value)
        if option in ("-r","--run_number"):
            run=int(value)
        if option in ("-t","--threads"):
            ncpus=int(value)
        if option in ("-h","--help"):
            usage()
            sys.exit()
    if not os.path.exists(input_path):
        sys.exit("Error - input fastq files path not found:\n"+input_path)
    if instID == ".":
        sys.exit("Error - instrument ID is mandatory...\n")
    if run == 0:
        sys.exit("Error - Run number is mandatory...\n")
    if output_path == ".":
        output_path = input_path
    elif not os.path.exists(output_path):
        sys.exit("Error - output path not found:\n"+output_path)

#    return gz, input_path, output_path, instID, run, ncpus
    return gz, input_path, output_path, instID, run

def usage():
    print "\n------------------------------------------------------------------------"
    print "convertFastqHeader.py will convert fastq files with MGI-formated headers"
    print " into fastq files with Illumina-formated headers."
    print "Since the barcode sequence is part of the Illumina header, the input fastq"
    print " files have to be demultpixed prior conversion i.e. needs to have the"
    print " sample barcode in the header."
    print "The instrument ID and the run number have to be passed in parameter"
    print " since they are needed for the Illumina header but not provided by the"
    print " MGI header."
    print "Converted output is always written in compressed fastq.qz format."
    print "This program was written by Edouard HENRION."
    print "For more information, contact: edouard.henrion@computationalgenomics.ca"
    print "------------------------------------------------------------------------\n"
    print "USAGE : fastqIcounter.py [option] "
    print "       -i :        path of the folder containing the fastq files to convert"
    print "       -o :        path of the folder where to write the converted fastq files,"
    print "                   named by adding '.converted' suffix at the input base name)"
    print "                   (default: output in the same directory as input)"
    print "       -g :        gziped input yes: 0 ; no: 1 (Default 0)"
    print "       -r :        run number"
    print "       -s :        instrument/sequencer ID"
#    print "       -t :        number of threads (default 4)"
    print "       -h :        this help \n"

def get_output_file(
    in_file,
    out_path,
    gz
    ):
    out_file = os.path.splitext(os.path.basename(in_file))[0]
    if gz == 0 :
        out_file = os.path.splitext(out_file)[0]
    return os.path.join(out_path, out_file + ".converted.fastq.gz")

def chunked_iterable(
    iterable,
    size
    ):
    it = iter(iterable)
    while True:
        chunk = tuple(islice(it, size))
        if not chunk:
            break
        yield chunk

def parse_MGI(
    header_string,
    instrument_ID,
    run_number
    ):

    m = re.search("(?P<fcid>V\w+)L(?P<lane>\d{1})C(?P<col>\d{3})R(?P<row>\d{3})(?P<tile>\d+):(?P<sample_barcode>[ATCG]+-?[ATCG]*)/(?P<read>\d{1})", header_string)
    if m:
        # Pack into dictionary
        header_dict = {
            "instrument": instrument_ID,
            "run_number": run_number,
            "sample_barcode": m.group('sample_barcode'),
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

def convert_Illumina(
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

def process_MGI_fastq(
    mgi_records,
    instrument_ID,
    run_number
    ):

    illumina_records = []
    for mgi_record in mgi_records:
        mgi_head, seq, qual = mgi_record

        illumina_head = convert_Illumina(
            parse_MGI(
                mgi_head,
                instrument_ID,
                run_number
            )
        )
        illumina_records.append("@%s\n%s\n+\n%s\n" % (illumina_head, seq, qual))

    return illumina_records

def main():
#    gz, input_path, output_path, instrument_ID, run_number, ncpus = getarg(sys.argv)
    gz, input_path, output_path, instrument_ID, run_number, = getarg(sys.argv)

    # the script looks for the following file types
    file_types = ('*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq')
    input_files = []
    for file_type in file_types:
        input_files.extend(glob(input_path + file_type))

    output_files = [get_output_file(in_file, output_path, gz) for in_file in input_files]
    
    if gz == 0 :
        open_input = compile("gzip.open", "<string>", "eval")
    else :
        open_input = compile("open", "<string>", "eval")

    in_fastq_handles = [(eval(open_input)(i, 'rb')) for i in input_files]
    print "fastqs opened for reading :\n" + "\n".join(input_files)

    out_fastq_handles = [gzip.open(i, 'wt') for i in output_files]
    print "fastqs opened for writing :\n" + "\n".join(output_files)

#    if ncpus > 1:
#        # Parallelization
#        # create a pool of processing nodes
#        p = Pool(ncpus)
#
#        # Loop over chunks of 40000 lines
#        for mgi_record_chunk in list(chunked_iterable(izip(*(FastqGeneralIterator(fh) for fh in in_fastq_handles)), 40000)):
#            # use partial to create a function needing only one argument
#            illumina_record_chunk = list(p.imap(
#                partial(
#                    process_MGI_fastq,
#                    instrument_ID,
#                    run_number
#                ),
#                mgi_record_chunk
#            ))
#            for illumina_records in illumina_record_chunk:
#                [fh_out.write(illumina_records[i]) for i, fh_out in enumerate(out_fastq_handles)]
#
#    else:
#        # Serial mode
    for mgi_records in list(izip(*(FastqGeneralIterator(fh) for fh in in_fastq_handles))):
        illumina_records = process_MGI_fastq(
            mgi_records,
            instrument_ID=instrument_ID,
            run_number=run_number
        )
        [fh_out.write(illumina_records[i]) for i, fh_out in enumerate(out_fastq_handles)]

    print "Fastq headers converted successfully !"

main()
