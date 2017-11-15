#!/usr/bin/env python

### Eric Audemard (2012/07/21) - eric.audemard@mail.mcgill.ca

import os
import sys
import string
import getopt
import re


def getarg(argument):
	scaF=""
	lenF=""
	optli,arg = getopt.getopt(argument[1:],"s:l:h",['scaffold','length','help'])
	if len(optli) == 0 :
		usage()
		sys.exit("Error : No argument given")
	for option, value in optli:
		if option in ("-s","--scaffold"):
			scaF=str(value)
		if option in ("-l","--length"):
			lenF=str(value)
		if option in ("-h","--help"):
			usage()
			sys.exit()
	if not os.path.exists(scaF) :
                sys.exit("Error - scaffold result file not found:\n"+str(scaF))
	if not os.path.exists(lenF) :
		sys.exit("Error - lenght file not found:\n"+str(lenF))
	return scaF, lenF
	
def addLegth(sF, lF):
	s=open(sF,'r')
	l=open(lF,'r')
	out=open(sF+".length",'w')
	ls=s.readline()
	ll=l.readline()
	while ls != "" :
		if ls[0] == ">" :
			cl=ll.split()
			cs=ls.split()
			out.write(cs[0]+" "+cl[1]+"\n")
			ll=l.readline()
		else :
			out.write(ls)
		ls=s.readline()
	out.write("\n")
	s.close()
	l.close()
	out.close()


def usage():
	print "USAGE : addLength.py [option] "
	print "       -s :        scaffold file"
	print "       -l :        lenght file"
	print "       -h :        this help \n"


def main():
	print "\n---------------------------------------------------------------------------------"
	print "addLength.py add lenght of each scaffold build by Ray. (fasta file like abyss)"
	print "This program was written by Eric Audemard"
	print "For more information, contact: eric.audemard@mail.mcgill.ca"
	print "----------------------------------------------------------------------------------\n"
	sf, lf = getarg(sys.argv)
	addLegth(sf, lf)
	
main()
	
	
