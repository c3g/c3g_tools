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
import warnings
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input_files", help="List of files to merge" , type=str, nargs="+", required=True)
    parser.add_argument("-d", "--delimiter", help="Delmiter for input files" , type=str, nargs="+", required=False, default='\t')    
    parser.add_argument("-o", "--output", help="Output File prefix", type=str , required=True)
    parser.add_argument("-c", "--common", help="Common columns used to join, if more than one, use comma separated column names, i.e. gene,transcript. By default, the first column is used as key", nargs="+",  type=str , required=False)    
    parser.add_argument("-s", "--subset", help="A subset of column names to print in the output file", nargs="+",  type=str , required=False, default=None)
    parser.add_argument("-x", "--exclude", help="A subset of column names to exclude from  the output file", nargs="+",  type=str , required=False, default=None)
    parser.add_argument("-l", "--left", help="If selected, left outer join", action="store_true")
    parser.add_argument("-t", "--sort", help="A column name to sort by", type=str , required=False)
    parser.add_argument("-f", "--filter", help="An python expression using any of the fields in the input files", nargs="+",  type=str , required=False)
    parser.add_argument("-n", "--make_names", help="If TRUE then the names of the columns in the input files that are duplicated are adjusted so that they are all preserved. If false, values are updated every time an input file is loaded", action="store_true")
    args = parser.parse_args()
    
    #error found when testing, don't ask me if system size is exceeded
    csv.field_size_limit(sys.maxsize)
    
    # Parameters
    infiles=args.input_files
    outfile=args.output
    subset=args.subset
    exclude=args.exclude
    
    # Separator for output file is \t
    out_sep='\t'
    in_sep=[]
    
    # Separator , common must have the same length as infiles
    if len(args.delimiter) == len(infiles) and len(args.delimiter) > 1:
        in_sep=args.delimiter
    elif len(args.delimiter) == 1 :
        in_sep=[args.delimiter[0]] * len(infiles)
    else:
        raise Exception("Error: Two or more input files are required " + args.delimiter )
 
    if not args.common is None and len(args.common) == len(infiles) and len(args.common) > 1:        
        key=args.common
    elif not args.common is None  and len(args.common) == 1:
        key=[args.common[0]] * len(infiles)
    else:
        key=[None] * len(infiles)
    
    # Treat the None case
    if "None" in key:
        key=map(lambda x:x if x is x!= 'None' else None,key)
    
    print "NOTICE: merging " + " ".join(infiles) + ", delimited by: " + str(in_sep) + ", with keys: " +  str(key) 
    #print str(key)
    #print str(args.left)
    # Open output files    
    data = collections.OrderedDict()
    fieldnames = []
    replacement_names={}
    
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
    i=0
    # Filter 
    if args.filter:
        filter_updated = " ".join(args.filter)
        print "NOTICE: filtering input files " + str(infiles) + " using the expression " + filter_updated

    for filename in infiles:
        with open(filename, "rb") as fp: # python 2
            if in_sep[i] in ['\t',"\\t"]:                
                reader = csv.DictReader(fp, quoting=csv.QUOTE_NONE, dialect=csv.excel_tab)
            else:
                reader = csv.DictReader(fp, delimiter=in_sep[i], quoting=csv.QUOTE_NONE)
            if args.make_names:
                new_names = dict((rn,rn + "_" + str(i)) if rn in fieldnames else (rn,rn) for rn in reader.fieldnames)
                #print str([new_names.values()])
                #print str([new_names.keys()])
                if i==0 :
                    replacement_names = { v:k for k, v in new_names.items() if k in fieldnames }
                else:
                    replacement_names.update( { v:k for k, v in new_names.items() if k in fieldnames } )
                reader.fieldnames = [new_names[fn] for fn in reader.fieldnames]
            fieldnames.extend(reader.fieldnames)
            for row in reader:
                # Filter data         
                if args.filter:                                        
                    # Search for variables in reader 
                    re_fields=dict(( fn, re.compile(fn)) for fn in reader.fieldnames)
                    if all(x is None for x in [re_fields[k].search(filter_updated) for k in row.keys()]):
                        warnings.warn("The expression doesn't include any of the queried fields: " + ",".join(reader.fieldnames) + " in " + filename )
                    commands=["def validate(): "]                    
                    filter_all=multiple_replace(filter_updated, row)
                    commands.append( "return(" + filter_all + ")"  + '\n')                    
                    exec '\n\t'.join(commands)
                    #print filter_all 
                    #print validate()
                    # If expression is not true, the next line is read
                    if not validate():
                        continue                             
                key_field = out_sep.join(row[k] for k in key[i].split(",") ) if key[i] else row[reader.fieldnames[0]]
                # Default behavior is cross join (common elements of all tables are added)
                # if left outer join, preserve only keys of the first file
                if (args.left and i == 0 ) or not args.left:
                    data.setdefault(key_field , {}).update(row)
                elif data.has_key(key_field) :
                    data[key_field].update(row)
            del reader
        i+=1
        
    if not subset is None:
        diff = [x for x in subset if x not in fieldnames]
        if len(diff) > 0 :
            error_subset = " ".join(diff)        
            warnings.warn("WARNING: Columns to include are not in input files: all columns will be included " + error_subset )
        else:    
            fieldnames = subset
            if args.make_names:
                fieldnames.extend([k for k,v in replacement_names.items() if v in fieldnames])
    else:
      fieldnames = list(data.fromkeys(fieldnames))
    if exclude is None:
        diff = []        
    else:
        diff = [x for x in exclude if x not in fieldnames]
    
    #print str(exclude)
    #print str(diff)
    #print str(replacement_names)
    
    if exclude and len(diff) > 0 :
        error_exclude = " ".join(diff)        
        warnings.warn("WARNING: Columns to exclude are not in input files: " + error_exclude )
    elif exclude is None:
        pass
    else:        
        if args.make_names:
            exclude.extend([k for k,v in replacement_names.items() if v in exclude])
        fieldnames = [x for x in fieldnames if x not in exclude ]
    print "NOTICE: will preserve the following headers: " + str(fieldnames)
    # Sort by a selected field
    if args.sort and args.sort in fieldnames:
        tmpdata=collections.OrderedDict()
        try:
            # Try numeric conversion
            tmpdata=sorted(data.items(), key=lambda d: float(d[1][args.sort]) if d[1].has_key(args.sort) else None)            
            data=collections.OrderedDict(tmpdata)
        except ValueError:            
            tmpdata=sorted(data.items(), key=lambda d: d[1][args.sort] if d[1].has_key(args.sort) else None )
            data=collections.OrderedDict(tmpdata)

    # Print output file  
    with open(outfile, "wb") as fp:
        writer = csv.writer(fp, delimiter=out_sep)
        writer.writerow(fieldnames)
        for row in data.itervalues():
            writer.writerow([row.get(field, '') for field in fieldnames])
