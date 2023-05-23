#!/usr/bin/env python3

### Henrion Edouard (2017/07/20)
### updated on 2023/04/11 for Python3 compatibility
### get GC count by bin

import os
import sys
import string
import getopt
import re
from Bio import SeqIO
from Bio.SeqUtils import GC

def getarg(argument):
    ref_file = ""
    optli, arg = getopt.getopt(argument[1:], "s:r:o:h", ['size', 'ref_file', 'out_file', 'help'])
    if len(optli) < 3 :
        usage()
        sys.exit("Error : Missing argument(s)")
    for option, value in optli:
        if option in ("-s", "--size"):
            bin_size = int(value)
        if option in ("-r", "--ref_file"):
            ref_file = str(value)
        if option in ("-o", "--out_file"):
            out_file = str(value)
        if option in ("-h", "--help"):
            usage()
            sys.exit()
    if not os.path.exists(ref_file) :
        sys.exit(f"Error - reference file not found:\n{ref_file}")
    return bin_size, ref_file, out_file

def usage():
    print("USAGE : getFastaBinedGC.py [option] ")
    print("       -s, --size     :        bin size (bp)")
    print("       -r, --ref_file :        reference genome file ")
    print("       -o, --out_file :        output file name")
    print("       -h, --help     :        this help \n")

def main():
    print("\n----------------------------------------------------------------------------------")
    print("getFastaBinedGC.py will generate a bed files with the GC content of each bin from ")
    print("the given genome file.")
    print("For more information, contact: edouard.henrion@computationalgenomics.ca")
    print("----------------------------------------------------------------------------------\n")
    bin_size, ref_file, out_file = getarg(sys.argv)
    out = open(out_file, 'w', encoding="utf-8")
    out.write("#Chrom\tStart\tEnd\tGCcontent\tUnMappContent\n")
    for seq_record in SeqIO.parse(ref_file, "fasta"):
        chrom_size = len(seq_record)
        chrom_ID = str(seq_record.id)
        start = 0
        end = start + bin_size - 1
        while start <= chrom_size :
            if end > chrom_size:
                end = chrom_size
            ref_seq = seq_record.seq[start:end]
            GC_content = int(GC(ref_seq))
            ## write Output
            out.write(f"{chrom_ID }\t{start}\t{end+1}\t{GC_content}\n")
            start = start + bin_size
            end = start + bin_size - 1
    out.close()

main()
