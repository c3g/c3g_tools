#!/usr/bin/env python

### Edouard Henrion (2017/06/06)
### CpG_stats

import os
import sys
import getopt
import numpy
import csv

def getarg(argument):
    profile=""
    optli,arg = getopt.getopt(argument[1:], "i:o:h", ['input', 'output', 'help'])
    for option, value in optli:
        if option in ("-i", "--input"):
            profile = str(value)
        if option in ("-o", "--output"):
            out = str(value)
        if option in ("-h", "--help"):
            usage()
            sys.exit()
    if not os.path.exists(profile) :
        sys.exit("Error - input CpG profile file not found:\n" + profile)
    return profile, out

def usage():
    print("USAGE : CpG_coverageStats.py [option] ")
    print("       -i :        CpG combined profile input file")
    print("       -o :        output text file")
    print("       -h :        this help \n")

def main():
    print("\n---------------------------------------------------------------------------------")
    print("CpG_coverageStats.py will generate the mean & median CpG coverage from a given")
    print("combined CpG profile i.e. a bismark CpG report for which the strand information")
    print("would have been combined using the perl script 'methylProfile.bismark.pl' which can")
    print("be found within the mugqic_tools, as part of the perl tools.")
    print("This program was written by Edouard HENRION")
    print("For more information, contact: edouard.henrion@computationalgenomics.ca")
    print("----------------------------------------------------------------------------------\n")

    profile, outfile = getarg(sys.argv)

    coverage_list = []
    with open(profile, "rb") as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        for row in reader:
            coverage_list.append(float(row['total']))

    coverage_list.sort()

    median_coverage=numpy.median(numpy.array(coverage_list))
    mean_coverage=numpy.mean(numpy.array(coverage_list))

    out=open(outfile, 'w')
    out.write("""Median coverage\t{median}\nMean coverage\t{mean}\n""".format(
        median=median_coverage,
        mean=mean_coverage
    ))

main()
