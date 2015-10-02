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
import collections
import itertools
import operator
import re
import warnings

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='Parse trinotate output and produce genes and isoforms annotations files: extract blast symbol, Go terms and/or filter for any field in trinotate output')
    parser.add_argument("-r", "--report", help="Trinotate report (tsv format)", type=argparse.FileType('rb'), required=True)
    parser.add_argument("-i", "--item_column", help="Column name for the item to select (\"#gene_id\" for genes, \"transcript_id\" for transcripts)" , type=str, required=False, default= "#gene_id")
    parser.add_argument("-b", "--Top_BLASTX_hit", help="Column name for top blast hit, default is sprot_Top_BLASTX_hit", type=str , required=False, default= "sprot_Top_BLASTX_hit")
    parser.add_argument("-g", "--gene_ontology", help="Column name for gene_ontology, default is gene_ontology_blast", type=str, required=False, default= "gene_ontology_blast")
    parser.add_argument("-o", "--output", help="Output File prefix", type=str , required=True)
    parser.add_argument("-l", "--length_file", help="Path to the gene/transcript length file. If declared, only the longest transcript per gene is reported. (\t separated file with columns gene_id\tlength)", type=argparse.FileType('rb') , required=False)
    parser.add_argument("-f", "--filter", help="Filter annotations using a python expression, generate a filtered annotations file", type=str, nargs="+", required=False)
    args = parser.parse_args()
    
    csv.field_size_limit(sys.maxsize)
    
    # Parameters
    infile=args.report
    blast=args.Top_BLASTX_hit
    go=args.gene_ontology
    outfile=args.output
    key=args.item_column
    length_fname=args.length_file
    
    # Open output files
    outfile_blast=open(os.path.abspath(outfile) + "_blast.tsv" , "wb")
    outfile_go=open(os.path.abspath(outfile) + "_go.tsv" , "wb")
      
    # Read file    
    trinotatereader = csv.DictReader(infile,  delimiter='\t', quoting=csv.QUOTE_NONE, restval=".")    
    field_names_blast= [key] + ["Symbol"] + [vals for vals in trinotatereader.fieldnames if vals not in [key] ]
    if length_fname:
        field_names_blast= field_names_blast + ["longest_transcript_length", "longest_transcript_id"]
    field_names_go=[key] + [go]
    csvwriter = csv.DictWriter(outfile_blast, field_names_blast, extrasaction='ignore',  delimiter='\t', quoting=csv.QUOTE_NONE)
    csvwriter_go = csv.DictWriter(outfile_go, field_names_go, extrasaction='ignore',  delimiter='\t', quoting=csv.QUOTE_NONE)
    
    # Force to print only one line per transcript id
    transcript_id = trinotatereader.fieldnames[1]

    # Read transcript length file. If this file is defined, the output will be
    # a gene/isoform subset of Trinotate report with only the longest isoform per gene.
    longest_item = {}
    length_id = "length"
    
    if length_fname :
        length_map = collections.defaultdict(dict)
        lengthreader = csv.DictReader(length_fname,  delimiter='\t', quoting=csv.QUOTE_NONE)
        item_id = lengthreader.fieldnames[0]
        length_id = lengthreader.fieldnames[1]        
        for row in lengthreader:
            length_map[row[item_id]] = row
            length_map[row[item_id]].update(dict( gene_group = "_".join(row[item_id].split("_")[:2]) ))
        group_by='gene_group'
        # Painfully sorting by gene ID
        rows = sorted(length_map.values(), key=lambda d: d[group_by])
        # Grouping and obtaining the max length
        groups = itertools.groupby(rows, lambda d: d[group_by])
        for k,g in groups:
            longest_item[k]=max(g, key=lambda d: float(d[length_id]))
        #longest_item = dict( k, max(g, key=lambda d: float(d[length_id])) for k, g in groups)
        # Free 
        del lengthreader
        del groups
        del rows
        length_fname.close()
    #print "NOTICE: length of longest_item: " + str(len(longest_item))
    count =1
    # If filter is defined, filter elements given a set of filter definitions          
    if args.filter:
        transcript_filtered={}
        filter_updated=" ".join(args.filter)        
        outfile_filtered=open(os.path.abspath(outfile) + "_filtered.tsv" , "wb")
        filtered_writer = csv.DictWriter(outfile_filtered, [key], extrasaction='ignore',  delimiter='\t', quoting=csv.QUOTE_NONE)
        re_fields=dict(( fn, re.compile(fn)) for fn in trinotatereader.fieldnames)
    
    # This function will replace all occurrences of a variable (a column name) in the filter string 
    def multiple_replace(text, adict):
        rx = re.compile('|'.join(map(re.escape, adict)))
        def one_xlat(match):
            try:
                value = float(adict[match.group(0)])
                return adict[match.group(0)]
            except ValueError:        
                return "\"" + adict[match.group(0)] + "\""
                pass
        return rx.sub(one_xlat, text)

    # Write blast / go annotated results
    item_done={}
    csvwriter.writeheader()
    csvwriter_go.writeheader()
    for line in iter(trinotatereader):
        # If filter is defined, filter elements given a set of filter definitions 
        # Search for variables in trinotate reader and create a script that must be evaluated
        if args.filter:                        
            # Controls for duplicated lines (one transcript per gene)            
            if not transcript_filtered.has_key(line[transcript_id]):
                # Search for variables in reader
                if all(x is None for x in [re_fields[k].search(filter_updated) for k in line.keys()]):
                    warnings.warn("The expression doesn't include any of the queried fields: " + ",".join(trinotatereader.fieldnames) + " in " + filename )
                commands=["def validate(): "]                    
                filter_all=multiple_replace(filter_updated, line)
                commands.append( "return(" + filter_all + ")"  + '\n')                    
                exec '\n\t'.join(commands)
                #print filter_all 
                #print validate()
                # If expression is true, write gene/transcript to the output file
                if validate():
                    filtered_writer.writerow(line)
                    transcript_filtered[line[transcript_id]] = None
        # Controls for duplicated lines (one transcript per gene)
        if item_done.has_key(line[transcript_id]):
            continue
        else:
            item_done[line[transcript_id]] = None
        line["gene_group"] = "_".join(line[key].split("_")[:2])
        if len(longest_item) == 0 or ( len(longest_item) > 0  and longest_item.has_key(line["gene_group"]) and longest_item[line["gene_group"]][transcript_id] == line[transcript_id] ) :
            best_blast=line[blast].split("^")[0]           
            line["Symbol"] = best_blast        
            line["longest_transcript_length"] = longest_item[line["gene_group"]][length_id] if longest_item.has_key(line["gene_group"]) else None
            line["longest_transcript_id"] = longest_item[line["gene_group"]][transcript_id] if longest_item.has_key(line["gene_group"]) else None
            csvwriter.writerow(line)
            go_lines=string.split(line[go], "`")
            for go_line in go_lines:
                go_table=go_line.split("^")
                line[go]=go_table[0]
                if go_table[0] != ".":
                    csvwriter_go.writerow(line)
    del trinotatereader
    del csvwriter
    del csvwriter_go    