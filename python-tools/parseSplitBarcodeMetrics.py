#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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
#
# Author : Mareike Janiak
# Contact : mareike.janiak@computationalgenomics.ca
#

import pandas as pd
import argparse

def join_tables(
        metrics_table,
        samplesheet
        ):
    metrics = parse_inputs(metrics_table, "\t")
    metrics['BarcodeName'] = metrics['#Barcode'].str.replace('barcode','')
            
    samplesheet = parse_inputs(samplesheet, ",")
    samplesheet['BarcodeName'] = samplesheet[['Sample_Name','Sample_ID']].apply(clean_barcode_name, axis=1)
                        
    joined_table = pd.merge(samplesheet,metrics, how="outer", on='BarcodeName')
    
    return joined_table

def clean_barcode_name(val):
    barcode = val[0].replace(val[1],'').strip()
    return barcode.lstrip('_')

def parse_inputs(
    input_file,
    delimiter
    ):
    
    output = pd.read_csv(input_file, sep=delimiter)
    return(output)

def print_metrics(
        joined_table, 
        output
        ):
                                                        
        final_table = joined_table[['BarcodeName','Sample_ID', ' Correct', ' Corrected', ' Total',' Percentage(%)']]
        final_table = final_table.rename(columns={' Correct': 'Correct', ' Corrected': 'Corrected', ' Total': 'Total', ' Percentage(%)': 'Percentage'})
        final_table.to_csv(output, index=False, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Combine splitBarcode metrics file with sample sheet information.""")
    parser.add_argument("-i", "--input_file", help="Input BarcodeStat file from splitBarcode", type=argparse.FileType('r'), required=True)
    parser.add_argument("-s", "--sample_sheet", help="Samplesheet for run, containing sample IDs and barcode information.", type=argparse.FileType('r'), required=True)
    parser.add_argument("-o", "--output", help="Output metric file, result of the combination of all the input metrics files.", type=str, required=True)
    args = parser.parse_args()
                                                                                        
    # combine tables by barcode name
    metrics = join_tables(args.input_file, args.sample_sheet)
                                                                                                
    # Print combined metrics to output file
    print_metrics(
            metrics,
            args.output)
