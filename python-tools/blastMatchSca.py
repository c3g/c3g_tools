#!/usr/bin/env python

### Eric Audemard (2011/10/21) - eric.audemard@mail.mcgill.ca

import os
import sys
import string
import getopt
import re

from Bio import SeqIO
from Bio.Blast import NCBIXML



def getarg(argument):
	blastF=""
	fastaF=""
	optli,arg = getopt.getopt(argument[1:],"f:b:h",['fasta','blast','help'])
	if len(optli) == 0 :
		usage()
		sys.exit("Error : No argument given")
	for option, value in optli:
		if option in ("-f","--fasta"):
			fastaF=str(value)
		if option in ("-b","--blast"):
			blastF=str(value)
		if option in ("-s","--sclip"):
			sclipF=str(value)
		if option in ("-h","--help"):
			usage()
			sys.exit()
	if not os.path.exists(str(fastaF+".fasta")) :
		sys.exit("Error - fasta file not found:\n"+str(fastaF+".fasta"))
	if not os.path.exists(blastF) :
                sys.exit("Error - base of blast file not found:\n"+blastF)
	return fastaF,blastF

def extractMatchBis(fastaF,blastF):
  

	f=open(fastaF+".fasta",'r')
  
	print "fasta : "+fastaF+".fasta"
	
	out=open(fastaF+".fasta.blast",'w')
	outDetail=open(fastaF+".fasta.blast.detail",'w')
	
	nbfile = 0
	records = list(NCBIXML.parse(open(blastF)))
	print "xml : "+blastF
	nbfile = nbfile + 1
	r=0
	lf=f.readline()
	
	lengthA="0"
	lengthF=0
	db="None"
	evalue=["-1"]
	numM="0"
	strMatch=""
	ligne = ""
	typeM = "None"
	
	out.write("ID\tL_ID\tT_MATCH\tMATCH\tNUM_M\tEVALUE\n")
	outDetail.write("ID\tMATCH\tSCORE\tEVALUE\n")
	
	while len(lf)!=0 :
		title = ["None"]
		score = [-1]
		
		if records[r].alignments :
			listAl = list(records[r].alignments)
			itAl=0
			
			title = [listAl[itAl].title]
			db = title[itAl].split(",")
			db = db[0]
			lengthA = str(listAl[itAl].length)
			hsp = list(listAl[itAl].hsps)
			evalue = [str(hsp[0].expect)]
			strMatch = hsp[0].match
			numM = str(len(strMatch.split("|"))-1)
			score = [hsp[0].score]
			  
			itAl += 1
			if (itAl < len(listAl)) :
				hsp = list(listAl[itAl].hsps)
				
			while itAl < len(listAl) and score[0]==hsp[0].score : 
				title.append(listAl[itAl].title)
				score.append(hsp[0].score)
				evalue.append(str(hsp[0].expect))
				itAl += 1
				if (itAl < len(listAl)) :
					hsp = list(listAl[itAl].hsps)
			
			#count number hit match on virus
			nbVirus=0
			nbHuman=0
			nbOther=0
			
			name = title[0]
			if (name.find("virus")!=-1 or name.find("Virus")!=-1) :
				nbVirus += 1
			else :
				if (name.find("Homo sapiens")!=-1 or name.find("Human")!=-1 or name.find("homo sapiens")!=-1 or name.find("human")!=-1 or name.find("H.sapien")!=-1) :
					nbHuman += 1
				else :
					nbOther += 1
			
			if (nbVirus>0 and nbHuman==0 and nbOther==0) :
				typeM="Virus"
			if (nbHuman>0 and nbVirus==0 and nbOther==0) :
				typeM="Human"
			if (nbOther>0 and nbVirus==0 and nbHuman==0) :
				typeM="Other"
			if (nbVirus>0 and nbHuman>0 and nbOther==0) :
				typeM="VirusHuman"
			if (nbVirus>0 and nbOther>0 and nbHuman==0) :
				typeM="VirusOther"
			if (nbVirus>0 and nbOther>0 and nbHuman>0) :
				typeM="VirusHumanOther"
			if (nbOther>0 and nbHuman>0 and nbVirus==0) :
				typeM="HumanOther"
		else :
			db = "None"
			length = "0"
			evalue = ["-1"]
			numM="0"
			strMatch=""
			typeM = "None"
		
		ligne=lf.split("\n")
		reads = ligne[0].split(">")
		reads = reads[1].split(" ")
		
		lf=f.readline()
		while len(lf)!=0 and lf[0] != ">" :
			lengthF = lengthF + len(lf) - 1
			lf=f.readline()
		
		
		out.write(reads[0] + "\t" + str(lengthF) + "\t" + typeM + "\t" + db + "\t" + numM + "\t" + evalue[0] + "\n")
		it=0
		while it != len(title) :
			outDetail.write(reads[0] + "\t"+ title[it] + "\t" + str(score[it]) + "\t"+ evalue[it] + "\n")
			it += 1
		lengthF = 0
		r = r + 1
		
		

	f.close()
	out.close()
	outDetail.close()
	
	
def usage():
	print "USAGE : blastMatchSca.py [option] "
	print "       -f :        base of fasta file (...contigs.HQ for ...contigs.HQ.fasta)"
	print "       -b :        base of blast xml file (...blast for ...blast.0.xml)"
	print "       -h :        this help \n"


def main():
	print "\n----------------------------------------------------------------------------------"
	print "blastMatchSca.py extract first hit from xml to tab format (for PUURe)."
	print "This program was written by Eric AUDEMARD"
	print "For more information, contact: eric.audemard@mail.mcgill.ca"
	print "----------------------------------------------------------------------------------\n"
	ff,bf = getarg(sys.argv)
	print "Arguments : ok"
	extractMatchBis(ff,bf)
	print "Extract Match : ok"

main()

	    
	
	





