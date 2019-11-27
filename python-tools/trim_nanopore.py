#!/usr/bin/env python

########################################################################################
# Script to select a certain portion of a long-read and save to a separate file.
# For use with the GenPipes Nanopore pipeline.
#
# Will trim the long-reads a user defined number of bp (-s) starting from position 100
# The idea is to avoid any potential adaptors at the beginning of the read.
# Author: H. Galvez
########################################################################################

from Bio import SeqIO
import sys
import argparse


def main(argv):
    """
    Function to perform the trimming and write output to file.
    Trimming size is determined by the -s parameter.
    Trimming will begin at position 100 in the read to avoid
    potential adapters at the beginning of the read.
    Trimmed reads can be used for BLAST QC.
    """
    parser = argparse.ArgumentParser(prog='trim_nanopore.py',
                                     usage="trim_nanopore.py -i <input> -o <output> -s <size>")

    parser.add_argument('-i', action='store', dest='infile',
                        help='Input fastq')

    parser.add_argument('-o', action='store', dest='outfile',
                        help='Output fastq')

    parser.add_argument('-s', action='store', dest='size',
                        help='Trim size')

    results = parser.parse_args()

    trimmed_reads = []

    for seq_record in SeqIO.parse(results.infile, "fastq"):
        trimmed_reads.append(seq_record[100:(100 + int(results.size))])

    SeqIO.write(trimmed_reads, results.outfile, "fastq")

if __name__ == "__main__":
    main(sys.argv[1:])

