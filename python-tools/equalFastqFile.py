#!/usr/bin/env python

### Eric Audemard (2013/02/20) - eric.audemard@mail.mcgill.ca

import os
import sys
import string
import getopt
import gzip

def getarg(argument):
	refF=""
	fastqF=""
	optli,arg = getopt.getopt(argument[1:],"r:f:h",['ref','fastq','help'])
	if len(optli) == 0 :
		usage()
		sys.exit("Error : No argument given")
	for option, value in optli:
		if option in ("-r","--ref"):
			refF=str(value)
		if option in ("-f","--fastq"):
			fastqF=str(value)
		if option in ("-h","--help"):
			usage()
			sys.exit()
	if not os.path.exists(refF) :
                sys.exit("Error - ref file not found:\n"+str(refF))
	if not os.path.exists(fastqF) :
                sys.exit("Error - fastq file not found:\n"+str(fastqF))
        return refF, fastqF

def deleteRead(rf, ff): 
	refF = gzip.open(rf, "rb" )
	fastqF = gzip.open(ff, "rb" )
	
	of = ff.split(".fastq.gz")[0] + ".equal.fastq.gz"
	outF = gzip.open(of, "wb")
	
	ref=refF.readline()
	fastq=fastqF.readline()
	while ref != "" :
		while fastq != "" and ref != fastq :
			fastq=fastqF.readline()
			fastq=fastqF.readline()
			fastq=fastqF.readline()
			fastq=fastqF.readline()
		outF.writelines(fastq)
		fastq=fastqF.readline()
		outF.writelines(fastq)
		fastq=fastqF.readline()
		outF.writelines(fastq)
		fastq=fastqF.readline()
		outF.writelines(fastq)
		fastq=fastqF.readline()
		
		ref=refF.readline()
		ref=refF.readline()
		ref=refF.readline()
		ref=refF.readline()
	
	refF.close()
	fastqF.close()
	outF.close()
	


def usage():
	print "USAGE : equalfastqFile.py [option] "
	print "       -r :        reference fastq file (sort by name needed)"
	print "       -f :        fastq file (sort by name needed)"
	print "       -h :        this help \n"


def main():
	print "\n----------------------------------------------------------------------------------"
	print "equalfastqFile.py delete read in fastqFile not presente in refFile "
	print "This program was written by Eric AUDEMARD"
	print "For more information, contact: eric.audemard@mail.mcgill.ca"
	print "----------------------------------------------------------------------------------\n"
	rf, ff = getarg(sys.argv)
	deleteRead(rf, ff)

main()
