#!/usr/bin/env python

### Mareike Janiak (2023/08/08) - mareike.janiak@computationalgenomics.ca
### based on previous work by Edouard Henrion

import errno
import os
import sys
import json
import getopt
import shutil
import signal
import time
import random

def getarg(argument):
    json_file=""
    inputs=[]
    readset=""
    optli,arg = getopt.getopt(argument[1:],"j:i:r:h",['json_file','inputs','readset','help'])
    if len(optli) == 0 :
        usage()
        sys.exit("Error : No argument given")
    for option, value in optli:
        if option in ("-j","--json_file"):
            json_file=str(value)
        if option in ("-i","--inputs"):
            inputs.append(str(value))
        if option in ("-r","--readset"):
            readset=str(value)
        if option in ("-h","--help"):
            usage()
            sys.exit()

    if not json_file:
        sys.exit("Error - json_file parameter is required")
    elif json_file and not os.path.exists(json_file):
        sys.exit("Error - JSON file not found:\n" + str(json_file))
    if len(inputs) == 0:
        sys.exit("Error - input parameter is required")
    if not readset:
        sys.exit("Error - readset parameter is required")
    else:
        not_found = []
        for in_file in inputs:
            if not os.path.exists(in_file):
                not_found.append(in_file)
            if not_found:
                sys.exit("Error - input file(s) not found:\n " + "\n ".join(not_found))

    return json_file, inputs, readset

def getFileSize(
    input_file
    ):
    stat = os.stat(input_file)
    size = stat.st_size
    return size

def getFileSizeHash(
    inputs,
    dict_to_update
    ):

    fastq1 = ""
    fastq2 = ""
    bam = ""
    bai = ""

    for in_file in inputs:
        if "_R1_001.fastq.gz" in in_file:
            fastq1 = in_file
            fastq1_size = getFileSize(fastq1)
            dict_to_update['fastq_1']['size'] = fastq1_size
        elif "_R2_001.fastq.gz" in in_file:
            fastq2 = in_file
            fastq2_size = getFileSize(fastq2)
            dict_to_update['fastq_2']['size'] = fastq2_size
        elif ".bam" in in_file:
            bam = in_file
            bam_size = getFileSize(bam)
            dict_to_update['bam']['size'] = bam_size
        elif ".bai" in in_file:
            bai = in_file
            bai_size = getFileSize(bai)
            dict_to_update['bai']['size'] = bai_size
        else:
            sys.exit("Error - Unexpected input file found")

    return dict_to_update


def report(
    json_file,
    inputs,
    readset
    ):

    with open(json_file, 'r') as json_fh:
        run_report_json = json.load(json_fh)

    readsets = [readset]
    for readset in readsets:
        dict_to_update = run_report_json["readsets"][readset]
        new_dict = getFileSizeHash(inputs, dict_to_update)

        run_report_json["readsets"][readset] = new_dict

    file_size_str = json.dumps(run_report_json, indent=4)

    #Print to file
    with open(json_file, 'w') as out_json:
        out_json.write(file_size_str)

def lock(filepath):
    unlocked = True
    while unlocked:
        try:
            os.makedirs(filepath + '.lock')
        except OSError as exception:
            if exception.errno == errno.EEXIST and os.path.isdir(filepath + '.lock'):
                # The lock folder already exists, we need to wait for it to be deleted
                sleep_time = random.randint(1, 100)
                time.sleep(sleep_time)
                pass
            else:
                # An unexpected error has occured : let's stop the program and raise the error"
                raise exception
        else:
            # The lock folder was successfully created !"
            unlocked = False

def unlock(filepath):
    shutil.rmtree(filepath + '.lock', ignore_errors=True)

def usage():
    print("\n-------------------------------------------------------------------------------------")
    print("runProcessingFileSizesToJson.py" )
    print("This program was written by Mareike Janiak, based on work by Edouard Henrion")
    print("For more information, contact: mareike.janiak@computationalgenomics.ca")
    print("----------------------------------------------------------------------------------\n")
    print("USAGE : runProcessingMetricsToJson.py")
    print("    -j    JSON file to be updated")
    print("    -r    readset for which to fetch metrics and update the JSON")
    print("    -i    input metrics file being parsed to update the JSON (could be multiple input metrics files)")
    print("    -h    this help\n")

def main():
    json_file, inputs, readset = getarg(sys.argv)
    # finally (unlock) will execute even if exceptions occur
    try:

        # Make sure the json_file is unlock if process receives SIGTERM too (not python exception)
        def sigterm_handler(_signo, _stack_frame):
            unlock(json_file)
            sys.exit(0)
        signal.signal(signal.SIGTERM, sigterm_handler)

        # First lock the file to avoid multiple and synchronous writing attemps
        lock(json_file)
        print("Updating run processing report (" + json_file + ") for readset " + readset + " file sizes ")

        report(
            json_file,
            inputs,
            readset
        )

    finally:
        # Finally unlock the file
        unlock(json_file)

if __name__ == '__main__':
    main()
