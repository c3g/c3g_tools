#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import argparse
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--fasta_file", help="Original fasta to filter" , type=str, required=True)
    parser.add_argument("-l", "--filter_file", help="Filtered fasta headers file" , type=str,required=True)    
    parser.add_argument("-c", "--fasta_id_column", help="Column of fasta ID in the filtered fasta headers file" , type=int, required=False, default=0)        
    parser.add_argument("-o", "--output", help="Output files prefix" , type=str,required=True)    
    args = parser.parse_args()
    
    # Get list (hash table) of components to keep
    filter_file=open(args.filter_file,'r')
    line=filter_file.readline()
    
    components={}
    while line != "" :
        words=line.split()        
        if len(words) > args.fasta_id_column:
            components[words[args.fasta_id_column]]=1
        line=filter_file.readline()

    filter_file.close()

    # Filter contigs fasta 
    seq_records_FF=[]
    for seq_record in SeqIO.parse(args.fasta_file, "fasta"):
            if seq_record.id in components:
                    seq_records_FF.append(seq_record)

    # Write Fasta and tabular output
    output_file=open(args.output + ".fasta",'w')
    SeqIO.write(seq_records_FF, output_file, "fasta")
    output_file.close()
    output_file=open(args.output + ".tsv",'w')
    SeqIO.write(seq_records_FF, output_file, "tab")
    output_file.close()
