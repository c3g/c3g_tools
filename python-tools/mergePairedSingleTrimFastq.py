#!/usr/bin/env python

### Mathieu Bourgey (2015/10/01) - mathieu.bourgey@mcgill.ca

import os
import sys
import string
import getopt
import re
import gzip
import itertools
from Bio import SeqIO


def getarg(argument):
	optli,arg = getopt.getopt(argument[1:],"p:q:s:t:o:h",['paired1','paired2','single1','single2','output','help'])
	if len(optli) == 0 :
		usage()
		sys.exit("Error : No argument given")
	for option, value in optli:
		if option in ("-p","--paired1"):
			if os.path.exists(str(value)) :
				p1=str(value)
			else :
				sys.exit("Error - paired 1 fastq file not found:\n"+str(value))
		if option in ("-q","--paired2"):
			if os.path.exists(str(value)) :
				p2=str(value)
			else :
				sys.exit("Error - paired 2 fastq file not found:\n"+str(value))
		if option in ("-s","--single1"):
			if os.path.exists(str(value)) :
				s1=str(value)
			else :
				sys.exit("Error - single1 file not found:\n"+str(value))
		if option in ("-t","--single2"):
			if os.path.exists(str(value)) :
				s2=str(value)
			else :
				sys.exit("Error - single2 file not found:\n"+str(value))
		if option in ("-o","--output"):
			out=str(value)
		if option in ("-h","--help"):
			usage()
			sys.exit()
	return p1, p2, s1, s2, out

def usage():
	print "\n-----------------------------------------------------------------------------------------------------"
	print "mergePairedSingleTrimFastq.py will merge paired and single fastq outputed to one main fatsq that could be use by bwa mem -p"
	print "This program was written by Mathieu BOURGEY"
	print "For more information, contact: mathieu.bourgey@mcgill.ca"
	print "------------------------------------------------------------------------------------------------------\n"
	print "USAGE : mergePairedSingleTrimFastq.py [options] "
	print "mandatory options:"
	print "       -p :        paired 1 fastq file"
	print "       -q :        paired 2 fastq file"
	print "       -s :        single 1 fastq file"
	print "       -t :        single 2 fastq file"
	print "       -o :        output file"
	print "other option:"
	print "       -h :        this help \n"
	


def interleacePaired(iter1, iter2) :
	while True :
		yield iter1.next()
		yield iter2.next()

def concatSingle(iter1) :
	while True:
		yield iter1.next()



def main():
	p1, p2, s1, s2, outf = getarg(sys.argv)    
	fp1, fp2, fs1, fs2 = gzip.open(p1, 'rb'), gzip.open(p2, 'rb'),  gzip.open(s1, 'rb'), gzip.open(s2, 'rb') 
	records_pair = interleacePaired(SeqIO.parse(fp1, "fastq"), SeqIO.parse(fp2, "fastq"))
	records_single1=concatSingle(SeqIO.parse(fs1, "fastq"))
	records_single2=concatSingle(SeqIO.parse(fs2, "fastq"))
	outfile = gzip.open(outf, 'wb')
	count = SeqIO.write(itertools.chain(records_pair,records_single1,records_single2), outfile,  "fastq")
	print count
	outfile.close()
	fp1.close()
	fp2.close()
	fs1.close() 
	fs2.close()

main()
