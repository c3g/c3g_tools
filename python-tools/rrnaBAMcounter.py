#!/usr/bin/env python

import os
import sys
import string
import getopt
import re
import pysam
import gzip


def getarg(argument):
	typ="transcript"
	optli,arg = getopt.getopt(argument[1:],"i:g:o:t:h",['input','gtf','output','type','help'])
	if len(optli) < 3 :
		usage()
		sys.exit("Error : Missing argument(s)")
	for option, value in optli:
		if option in ("-i","--input"):
			inF=str(value)
		if option in ("-g","--gtf"):
			gtF=str(value)
		if option in ("-o","--output"):
			outF=str(value)
		if option in ("-t","--type"):
			typ=str(value).lower()
		if option in ("-h","--help"):
			usage()
			sys.exit()
	if not os.path.exists(inF) :
		sys.exit("Error - input file not found:\n"+inF)
	if typ != "transcript" and typ != "gene" :
		sys.exit("Error - wrong type:\n"+typ)
	
	return  inF, gtF, outF, typ



def usage():
	print "USAGE : rrnaBAMCounter.py [option] "
	print "       -i :        input file"
	print "       -o :        output file"
	print "       -g :        gtf file"
	print "       -t :        type of the instance rRNA: transcript / gene (default: transcript)"
	print "       -h :        this help \n"


## generate a mapping table detween rRNA instance id and rRNA name
def id_to_name(gf,t,n):
	idTag=t+"_id"
	nameTag=t+"_name"
	g=open(gf,'r')
	l=g.readline()
	map_dict={}
	while l != "" :
		if l[0] != "#" :
			c=l.replace(' "','\t').replace('"; ','\t').split("\t")
			if idTag in c and nameTag in c :
                                if c[c.index(idTag)+1] in n :
                                        if not map_dict.has_key(c[c.index(idTag)+1]) :
                                                map_dict[c[c.index(idTag)+1]] = c[c.index(nameTag)+1]
		l=g.readline()
	g.close()
	return map_dict



def main():
	print "\n------------------------------------------------------------"
	print "rrnaBAMCounter.py will estimate rRNA abundancy after "
	print "mapping reads on fasta containing rRNA sequeunces"
	print "This program was written by Mathieu BOURGEY"
	print "For more information, contact: mbourgey@mail.mcgill.ca"
	print "------------------------------------------------------------\n"
	inF, gtf, outF, typ = getarg(sys.argv)
	
	### open bam file
	samfile = pysam.AlignmentFile(inF,'rb')
	mappingTable=id_to_name(gtf,typ,samfile.references)
	## open output
	o=open(outF,'w')
	o2=open(outF+"_short.tsv",'w')
	## count in a hash :total reads; mapped read to rRNA; mapped read by rRNA instance in the fasta
	ct={}
	ct["total"]=0
	ct["rRNA_mapped"]=0
	for x in samfile:
 		ct["total"]+=1
 		if not x.is_unmapped :
  			ct["rRNA_mapped"]+=1
                        if not samfile.getrname(x.reference_id) in ct :
                                ct[samfile.getrname(x.reference_id)]=0
                        ct[samfile.getrname(x.reference_id)]+=1
        ##associate rRNA instance to count (include instance with no read mapped
	##and go find the Gene nama in the gtf
	rname=["Mapped","rRNA_mapped"]
	srname=["Mapped","rRNA_mapped","rRNA_mapped_%"]
	rcount=[str(ct["total"]),str(ct["rRNA_mapped"])]
	srcount=[str(ct["total"]),str(ct["rRNA_mapped"]),str(round((float(ct["rRNA_mapped"])/float(ct["total"]))*100,2))]
	for r in samfile.references:
                if str(r) in mappingTable :
                        rname.append(mappingTable[str(r)])
                else :
                        rname.append(str(r))
		if r in ct :
			rcount.append(str(ct[r]))
		else:
			rcount.append("0")	
	o.write("\t".join(rname)+"\n")
	o.write("\t".join(rcount)+"\n")
	o2.write("\t".join(srname)+"\n")
        o2.write("\t".join(srcount)+"\n")
	o.close()
	o2.close()
	samfile.close()

main()
