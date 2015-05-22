#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines tools.
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

import csv
import argparse
import string
import sys
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-r", "--report", help="Input Files", type=argparse.FileType('rb'), required=True)
    parser.add_argument("-i", "--item_column", help="Column name for the item to select (\"#gene_id\" for genes, \"transcript_id\" for transcripts)" , type=str, required=False, default= "#gene_id")
    parser.add_argument("-b", "--Top_BLASTX_hit", help="Column name for top blast hit, default is sprot_Top_BLASTX_hit", type=str , required=False, default= "sprot_Top_BLASTX_hit")
    parser.add_argument("-g", "--gene_ontology", help="Name of column for gene_ontology, default is gene_ontology_blast", type=str, required=False, default= "gene_ontology_blast")
    parser.add_argument("-o", "--output", help="Output File prefix", type=str , required=True)
    args = parser.parse_args()
    
    csv.field_size_limit(sys.maxsize)
    
    # Parameters
    infile=args.report
    blast=args.Top_BLASTX_hit
    go=args.gene_ontology
    outfile=args.output
    key=args.item_column
    
    #Open output files
    outfile_blast=open(os.path.abspath(outfile) + "_blast.tsv" , "wb")
    outfile_go=open(os.path.abspath(outfile) + "_go.tsv" , "wb")
    
    # Read file    
    reader = csv.DictReader(infile,  delimiter='\t', quoting=csv.QUOTE_NONE)    
    field_names_blast= [key] + [blast + "_ext"] + [vals for vals in reader.fieldnames if vals not in [key] ]
    field_names_go=[key] + [go]
    csvwriter = csv.DictWriter(outfile_blast, field_names_blast, extrasaction='ignore',  delimiter='\t', quoting=csv.QUOTE_NONE)
    csvwriter_go = csv.DictWriter(outfile_go, field_names_go, extrasaction='ignore',  delimiter='\t', quoting=csv.QUOTE_NONE)
    
    # Write blast annotated results
    csvwriter.writeheader()
    csvwriter_go.writeheader()
    for line in reader:
        best_blast=line[blast].split("^")[0]           
        line[blast + "_ext" ] = best_blast        
        csvwriter.writerow(line)
        go_lines=string.split(line[go], "`")            
        for go_line in go_lines:
            go_table=go_line.split("^")
            line[go]=go_table[0]
            if go_table[0] != ".":
                csvwriter_go.writerow(line)
    del reader
    del csvwriter
    del csvwriter_go