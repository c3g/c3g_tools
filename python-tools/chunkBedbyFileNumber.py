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

import argparse
import csv
import glob
import logging
import os
import re
import sys
import math
import pandas as pd

if sys.version_info[0] < 3:
    import httplib
else:
    import http.client

log = logging.getLogger(__name__)

def main():
    # Parse options
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", help="input bed file", type=argparse.FileType('r'), required=True)
    parser.add_argument("-o", "--output_dir", help="output_dir", required=True)
    parser.add_argument("-c", "--chunk", help="chunk size by number of files (optional ; Default no chunk)", type=int, default=0)
    args = parser.parse_args()

    input_file = args.input
    output_dir = args.output_dir
    
    bed_df = pd.read_csv(input_file, sep="\t", header=None)

    # Extract the chromosome, start position, and end position columns
    chrom_col = bed_df[0]
    start_col = bed_df[1]
    end_col = bed_df[2]

    # Calculate the total number of base pairs in the BED file
    total_bp = sum(end_col - start_col)

    # Calculate the desired number of base pairs per file
    num_files = args.chunk
    bp_per_file = math.ceil(total_bp // num_files)

    adjusted = total_bp - (bp_per_file * num_files)

    print(f"{total_bp} total bases spread of {num_files} files is {bp_per_file} per file")

    print(f"{adjusted} leftover ")

    # Initialize variables for tracking the current file and running total of base pairs
    current_file_num = 1
    current_bp = 0

    # Iterate over the intervals and write them to the appropriate output file
    total_current_bp = 0
    for i in range(len(chrom_col)):
        # Calculate the length of the current interval
        interval_bp = end_col[i] - start_col[i]

        # If the running total of base pairs exceeds the desired number of base pairs per file, start a new output file
        
        if current_bp + interval_bp > bp_per_file:
            total_current_bp += current_bp
            adjusted_bp_per_file = (total_bp - total_current_bp) // (num_files - current_file_num)
            print(f"Iterated through {total_current_bp} in {current_file_num} of files, adjusting bp per file to {adjusted_bp_per_file}")
            current_file_num += 1
            bp_per_file = adjusted_bp_per_file
            current_bp = 0

        # Write the interval to the current output file
        with open(f"{output_dir}/{current_file_num:04d}-scattered.bed", "a") as f:
            f.write(f"{chrom_col[i]}\t{start_col[i]}\t{end_col[i]}\n")
            current_bp += interval_bp

if __name__ == '__main__':
    main()
