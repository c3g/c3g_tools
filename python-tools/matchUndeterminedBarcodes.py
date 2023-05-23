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

import argparse
import re
import csv

#####
## Parse undetermined barcode file and look for matches in standard barcodes file. 
## Create table with read count, barcode sequence, and any match to standard barcodes
## For all unmatched barcodes with more than 5000 reads

def parse_undetermined_barcodes(undetermined_counts):
    barcodes = dict()

    with open(undetermined_counts, 'r') as f:
        for line in f:
            row = line.strip()
            read_count = row.split()[0]
            if int(read_count) > 50000:
                barcode = row.split()[1]
                barcodes[barcode] = read_count
    return barcodes


def match_undetermined_barcodes(barcodes, barcodes_by_sequence):
    rows = [["ReadCount", "Sequence", "Matches"]]
    with open(barcodes_by_sequence, 'r') as db:
        db = db.read()
        for barcode in barcodes:
            p = "{barcode}[ACTG]?[ACTG]?\t.*".format(barcode=barcode)
            match = re.findall(p, db)
            if match:
                match = [m.replace('\t', ' ') for m in match]
                rows.append([barcodes[barcode], barcode, '; '.join(match)])
            if not match:
                rows.append([barcodes[barcode], barcode, ""])
    return rows
    
def print_undetermined_barcodes_with_matches(rows, output_file):
    with open(output_file, 'w') as out:
        writer = csv.writer(out, delimiter = '\t')
        writer.writerows(rows)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Check common undetermined barcodes for matches in database of sequencing barcodes.")
    parser.add_argument("-i", "--input_file", help="File containing the counts and sequences of undetermined barcodes for the run.", type=str, required=True)
    parser.add_argument("-b", "--barcode_sequence_file", help="File containing database of sequencing barcode sequences and names.", type=str, required=True)
    parser.add_argument("-o", "--output", help="Output file with common undetermined barcodes and their matches to barcode database.", type=str, required=True)
    args = parser.parse_args()

    # parse undetermined barcodes count file for run
    barcodes = parse_undetermined_barcodes(
            args.input_file
            )

    # look for matches in barcode database for most common barcodes
    rows = match_undetermined_barcodes(
            barcodes,
            args.barcode_sequence_file
            )

    # print most common barcodes and any matches to output file
    print_undetermined_barcodes_with_matches(
            rows, 
            args.output
            )
