#!/usr/bin/env python

### Eric Audemard (2011/10/21) - eric.audemard@mail.mcgill.ca

import os
import sys
import string
import getopt
import re


def getarg(argument):
	sclipF=""
	length=21
	optli,arg = getopt.getopt(argument[1:],"s:l:h",['sclip','length','help'])
	if len(optli) == 0 :
		usage()
		sys.exit("Error : No argument given")
	for option, value in optli:
		if option in ("-s","--sclip"):
			sclipF=str(value)
		if option in ("-l","--length"):
			length=int(value)
		if option in ("-h","--help"):
			usage()
			sys.exit()
	if not os.path.exists(str(sclipF)) :
		sys.exit("Error - softclip file not found:\n"+str(sclipF))
	return sclipF, length

def sclipToFasta(sclipF, length):
  
	s=open(sclipF,'r')
	
	of = sclipF.split(".txt")[0] + ".fasta"
	ofName = sclipF.split(".txt")[0] + ".name.fasta"
	out=open(of,'w')
	outName=open(ofName,'w')
	
	ls=s.readline() # read header 
	ls=s.readline() # read first sclip
	
	if len(ls)==0 :
		s.close()
		out.close()
		return
	
	tabS = ls.split("\t")
	lChrom = tabS[0]
	lPos = tabS[1]
	seqKeep = [tabS[4]]
	name=[tabS[3]]
	ls=s.readline()
	
	it = 0
	boo = 1
	
	while len(ls)!=0 :
		tabS = ls.split("\t")
		
		chrom = tabS[0]
		pos = tabS[1]
		
		if (len(tabS[4])>=length) :
			if lChrom == chrom and lPos == pos :
				it = 0
				boo = 1		# =0 when we find an exact substring
				while (boo==1 and it<len(seqKeep)) :
					if (len(seqKeep[it])>len(tabS[4])) : 
						if (seqKeep[it].find(tabS[4])!=-1) :
							boo = 0
					else :
						if (tabS[4].find(seqKeep[it])!=-1) :
							boo = 0
					it += 1
				
				if (boo == 1) : # new sequence for this position
					seqKeep.append(tabS[4])
					name.append(tabS[3])
				else : # add the substring
					it -= 1
					name[it] = name[it]+"|"+tabS[3]
					if (len(seqKeep[it])<len(tabS[4])) :
						seqKeep[it] = tabS[4]
			else :	#Change sclip
				it = 0
				while (it<len(seqKeep)) :
					out.write(">"+name[it]+"\n")
					out.write(seqKeep[it]+"\n")
					outName.write(">"+name[it]+"\n")
					it +=1
					
				lChrom = tabS[0]
				lPos = tabS[1]
				seqKeep = [tabS[4]]
				name=[tabS[3]]
		
		ls=s.readline()
		
	
	#print the last one
	it = 0
	while (it<len(seqKeep)) :
		out.write(">"+name[it]+"\n")
		out.write(seqKeep[it]+"\n")
		outName.write(">"+name[it]+"\n")
		it += 1
		
	s.close()
	out.close()
	
	
def usage():
	print "USAGE : sclipToFasta.py [option] "
	print "       -s :        softclip file (with sequence)"
	print "       -l :        minimun length (default : 21)"
	print "       -h :        this help \n"


def main():
	print "\n----------------------------------------------------------------------------------"
	print "sclipToFasta.py create fasta file from a softclip file."
	print "When a softclip is a substring at an other, this softclip is delete of the fasta file"
	print "This program was written by Eric AUDEMARD"
	print "For more information, contact: eric.audemard@mail.mcgill.ca"
	print "----------------------------------------------------------------------------------\n"
	sf,l = getarg(sys.argv)
	print "Arguments : ok"
	sclipToFasta(sf,l)
	print "Fasta : ok"

main()






