#!/usr/bin/env python

### Henrion Edouard (2017/07/20)
### get GC count by bin

import os
import sys
import string
import getopt
import re
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

def getarg(argument):
    refF=""
    optli,arg = getopt.getopt(argument[1:], "s:r:m:o:h", ['siz', 'ref', 'out', 'help'])
    if len(optli) < 3 :
        usage()
        sys.exit("Error : Missing argument(s)")
    for option, value in optli:
        if option in ("-s","--siz"):
            bsi=int(value)
        if option in ("-r","--ref"):
            refF=str(value)
        if option in ("-o","--output"):
            out=str(value)
        if option in ("-h","--help"):
            usage()
            sys.exit()
    if not os.path.exists(refF) :
        sys.exit("Error - reference file not found:\n"+refF)
    return  bsi, refF, out

def usage():
    print "USAGE : getFastaBinedGC.py [option] "
    print "       -s :        bin size (bp)"
    print "       -r :        reference genome file "
    print "       -o :        output basename"
    print "       -h :        this help \n"

def main():
    print "\n----------------------------------------------------------------------------------"
    print "getFastaBinedGC.py will generate a bed files with the GC content of each bin from "
    print "the given genome file."
    print "This program was written by Edouard HENRION"
    print "For more information, contact: edouard.henrion@computationalgenomics.ca"
    print "----------------------------------------------------------------------------------\n"
    siz, refF, outF = getarg(sys.argv)
    out=open(outF,'w')
    out.write("#Chrom\tStart\tEnd\tGCcontent\tUnMappContent\n")
    for seq_record in SeqIO.parse(refF, "fasta"):
        chroS=len(seq_record)
        chroID=str(seq_record.id)
        start=0
        end=start+siz-1
        while start <= chroS :
            if end > chroS:
                end=chroS
            refSeq=seq_record.seq[start:end]
            GCcont=int(GC(refSeq))
            ## write Output
            out.write(chroID +"\t"+str(start)+"\t"+str(end+1)+"\t"+str(GCcont)+"\n")
            start=start+siz
            end=start+siz-1
    out.close()

main()
