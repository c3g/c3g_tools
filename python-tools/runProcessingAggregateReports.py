#!/usr/bin/env python

### Edouard Henrion (2020/03/18) - edouard.henrion@computationalgenomics.ca

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
        main_json=""
        sample_json_directory=""
        output_file=""
        optli,arg = getopt.getopt(argument[1:],"m:s:o:h",['main_json','sample_json_directory','output_file','help'])
        if len(optli) == 0 :
            usage()
            sys.exit("Error : No argument given")
        for option, value in optli:
            if option in ("-m", "--main_json"):
                main_json=str(value)
            if option in ("-s","--sample_json_directory"):
                sample_json_directory=str(value)
            if option in ("-o","--output_file"):
                output_file=str(value)
            if option in ("-h","--help"):
                usage()
                sys.exit()

        if not os.path.exists(main_json):
            sys.exit("Error - main JSON file not found:\n" + str(main_json))
        if not os.path.exists(sample_json_directory):
            sys.exit("Error - sample JSON report folder not found:\n" + str(sample_json_directory))

        return main_json, sample_json_directory, output_file

def usage():
        print "\n-------------------------------------------------------------------------------------"
        print "runProcessingAggregateReport.py takes all the individual JSON sample reports"
        print "and aggregates them into one single JSON sun report"
        print "This program was written by Edouard Henrion"
        print "For more information, contact: edouard.henrion@computationalgenomics.ca"
        print "----------------------------------------------------------------------------------\n"
        print "USAGE : runProcessingReport.py"
        print "    -m   JSON file containing general information to be iupdated with sample json reports"
        print "    -s   path of the folder containing all the sample JSON reports"
        print "    -o   output_file"
        print "    -h   this help\n"

def main():

        main_json_file, sample_json_directory, output_file = getarg(sys.argv)

        run_validation_hash = {}
        with open(main_json_file, 'r') as main:
            run_validation_json = json.load(main)
            run_validation_json['run_validation'] = []

        for filename in os.listdir(sample_json_directory):
            filepath = os.path.join(sample_json_directory, filename)
            sample = re.sub(".report.json", "", filename)

            with open(filepath, 'r') as json_file:
                sample_stats_json = json.load(json_file)

            run_validation_json['run_validation'].append(sample_stats_json)

        # Compute some lane-wise metrics
        pf_clusters_count_list = [int(sample['index']['PF Clusters']) for sample in run_validation_json['run_validation']]
        run_validation_json['total_pf_clusters'] = sum(pf_clusters_count_list)
        max_pf_cluster = max(pf_clusters_count_list)
        min_pf_cluster = min(pf_clusters_count_list)
        run_validation_json['spread'] = float(max_pf_cluster)/min_pf_cluster

        # Print to file
        with open(output_file, 'w') as out_json:
            json.dump(run_validation_json, out_json, indent=4)

main()

