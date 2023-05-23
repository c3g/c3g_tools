#!/usr/bin/env python

### Edouard Henrion (2020/09/01) - edouard.henrion@computationalgenomics.ca

import os
import sys
import getopt
import re
import json
import csv
from collections import namedtuple
import xml.etree.ElementTree as ET
from bs4 import BeautifulSoup

def getarg(argument):
        readset_json=""
        stats_json}=""
        output_json=""
        optli,arg = getopt.getopt(argument[1:],"r:s:o:h",['readset_json','stats_json','output_json','help'])
        if len(optli) == 0 :
            usage()
            sys.exit("Error : No argument given")
        for option, value in optli:
            if option in ("-r","--readset_json"):
                readset_json=str(value)
            if option in ("-s","--stats_json"):
                stats_json=str(value)
            if option in ("-o","--output_json"):
                output_json=str(value)
            if option in ("-h","--help"):
                usage()
                sys.exit()

        if not os.path.exists(readset_json):
            sys.exit("Error - main JSON file not found:\n" + str(readset_json))
        if not os.path.exists(stats_json):
            sys.exit("Error - index JSON file not found:\n" + str(stats_json))

        return readset_json, stats_json, output_json

def usage():
        print "\n----------------------------------------------------------------"
        print "runProcessingEditIndexStatsFile.py extracts the information from"
        print "the BCL2fastq Stats.json output file, adds the barcode names along"
        print "the barcode sequences using the provided readset JSON file and"
        print "and prints the result in JSON format to output file"
        print "This program was written by Edouard Henrion"
        print "For more information, contact: edouard.henrion@computationalgenomics.ca"
        print "----------------------------------------------------------------------------------\n"
        print "USAGE : runProcessingReport.py"
        print "    -r   JSON file containing index per readset information, built by the pipeline"
        print "    -s   JSON file outputed by BCL2fastq (generally in <...>/Stats/Stats.json)"
        print "    -o   output file (JSON format)"
        print "    -h   this help\n"

def main():

        readset_json_file, stats_json_file, output_json_file = getarg(sys.argv)

        with open(readset_json_file, 'r') as rjf:
            index_per_sample_json = json.load(rjf)

        with open(stats_json_file, 'r') as sjf:
            stats_json = json.load(sjf)

        for i, sample in enumerate(stats_json['ConversionResults'][0]['DemuxResults']):
            if index_per_sample_json[sample['SampleName']]:
                sample['IndexMetrics'][0]['IndexName'] = index_per_sample_json[sample['SampleName']][0]['INDEX_NAME']

        # Print to file
        with open(output_json_file, 'w') as ojf:
            json.dump(stats_json, ojf, indent=4)

main()

