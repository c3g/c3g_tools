#!/usr/bin/python

### Mathieu Bourgey (2011/03/21)
### extract sequence based on a given list

import os
import sys
import string
import getopt
import re
import itertools
from Bio import SeqIO
import gzip
import operator


def getarg(argument):
	rlist=""
	fastq1=""
	fastq2="."
	gz=0
	of="indexCount.tsv"
	optli,arg = getopt.getopt(argument[1:],"i:m:o:g:h",['input','mate','out','gzip','help'])
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
                if option in ("-o","--out"):
			of=str(value)
		if option in ("-h","--help"):
			usage()
			sys.exit()
	if not os.path.exists(fastq1) :
		sys.exit("Error - fastq file not found:\n"+fastq1)
	if fastq2 != ".":
		if not os.path.exists(fastq2) :
			sys.exit("Error - fastq file not found:\n"+fastq2)
	return  gz, fastq1, fastq2, of



def usage():
	print "USAGE : fastqIcounter.py [option] "
	print "       -i :        fastq input file"
	print "       -m :        mate fastq file (optional)"
	print "       -g :        gziped input yes: 0 ; no: 1 (Default 0)"
	print "       -o :        output file name (indexCount.tsv)"
	print "       -h :        this help \n"


def main():
	print "\n------------------------------------------------------------"
	print "fastqIcounter.py will extract the reads in a fastq index file" 
	print " and count the number of each index"
	print "This program was written by Mathieu BOURGEY"
	print "For more information, contact: mathieu.bourgey@mcgill.ca"
	print "------------------------------------------------------------\n"
	gz, fastq1, fastq2, of = getarg(sys.argv)
	idxfound={}
	if gz == 0 :
                f1=gzip.open(fastq1, 'rb')
		if fastq2 != ".":
			f2=gzip.open(fastq2, 'rb')
	else :
		f1=fastq1
		if fastq2 != ".":
			f2=fastq2
        print "fastq(s) opened"
	if fastq2 != ".":
		for  record1,record2 in itertools.izip(SeqIO.parse(f1, "fastq"),SeqIO.parse(f2, "fastq")):
                    if idxfound.has_key(str(record1.seq)+"-"+str(record2.seq)) :
                        idxfound[str(record1.seq)+"-"+str(record2.seq)]+=1
                    else :
                        idxfound[str(record1.seq)+"-"+str(record2.seq)]=1
	else :
		for  record1 in SeqIO.parse(f1, "fastq"):
                    if idxfound.has_key(str(record1.seq)) :
                        idxfound[str(record1.seq)]+=1
                    else :
                        idxfound[str(record1.seq)]=1
        print "indexes readed and counted"
        out=open(of,'w')
        k=idxfound.keys()
        for i in k :
            out.write(str(i) + "\t" + str(idxfound[i]) + "\n")
        out.close()
        print "rev-sorted index outpued in:"
        print of


main()

	    
	
	





