#!/usr/bin/env python

### Mathieu bourgey (2017/09/21) - mathieu.bourgey@.mcgill.ca

import os
import sys
import string
import getopt
import re
import csv


def getarg(argument):
	readSet_file=""
	project_file=""
	output_file="readset_recoded.tsv"
	optli,arg = getopt.getopt(argument[1:],"r:p:o:h",['readset','project','output','help'])
	if len(optli) == 0 :
		usage()
		sys.exit("Error : No argument given")
	for option, value in optli:
		if option in ("-r","--readset"):
			readSet_file=str(value)
		if option in ("-p","--project"):
			project_file=str(value)
                if option in ("-o","--output"):
                        output_file=str(value)
		if option in ("-h","--help"):
			usage()
			sys.exit()
	if not os.path.exists(readSet_file) :
                sys.exit("Error - readset file not found:\n"+str(readSet_file))
	if not os.path.exists(project_file) :
		sys.exit("Error - project file not found:\n"+str(project_file))
	return readSet_file, project_file,  output_file
	

def extractAlias(pf):
    dico={}
    f=open(pf,'rb')
    l=csv.reader(f)
    col_alias=None
    col_name=None
    col_barcode=None
    col_run=None
    col_lane=None
    col_status=None
    ct=0
    for row in l:
      if ct == 0 :
          #work on header: Catch the alias column and the path collumn (fastqR1 or BAM)
          for i in range(0,len(row)):
              if row[i] == "Library Name" :
                col_alias=i
              elif row[i] == "Name" :
                col_name=i
              elif row[i] == "Library Barcode" :
                col_barcode=i
              elif row[i] == "Status" :
                col_status=i
              elif row[i] == "Run" :
                col_run=i
              elif row[i] == "Region" :
                col_lane=i
      else :
          #work on csv body: extract Alias and the corresponding basename of R1 and/or bam 
          if row[col_status] == "Data is valid":
            dico[".".join([row[col_name],row[col_barcode],row[col_run],row[col_lane]])]=row[col_alias]
      ct+=1
    return dico

def usage():
	print "USAGE : readset_recodeBy_Alias.py [option] "
	print "       -r :        readset file"
	print "       -p :        project file"
	print "       -o :        output file (Default readset_recoded.tsv)"
	print "       -h :        this help \n"


def main():
	print "\n---------------------------------------------------------------------------------"
	print "readset_recodeBy_Alias.py extract the Alias information in the nanuq porject file" 
	print "and replace the sample name by the alias in the readset file"
	print "This program was written by Mathieu Bourgey"
	print "For more information, contact: mathieu.bourgey@mcgill.ca"
	print "----------------------------------------------------------------------------------\n"
	rf, pf, output = getarg(sys.argv)
	alias=extractAlias(pf)
	f=open(rf,'rb')
	o=open(output,'wb')
	l=f.readline()
	o.write(l)
	col_Sample=None
	col_readset=None
        c=l.split()
        for i in range(0,len(c)):
            if c[i] == "Sample" :
                col_Sample=i
            elif c[i] == "Readset" :
                col_readset=i
        l=f.readline()
        ct=1
        while l != "" :
            ct+=1
            c=l.split("\t")
            c_term=c[-1:][0].split()
            c[-1:]=c_term
            c[col_Sample]=alias[c[col_readset]]
            o.write("\t".join(c)+"\n")
            l=f.readline()
        f.close()
        o.close()
	
main()
	
	
