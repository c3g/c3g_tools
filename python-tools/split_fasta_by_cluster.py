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

from argparse import ArgumentParser


def parse_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--input-fasta', help='fasta file to split')
    arg_parser.add_argument('--min-lines', help='minimum number of lines per file split into')
    arg_parser.add_argument('--padding', help='padding for suffix of split files') 

    return arg_parser.parse_args()


def split_fasta_by_cluster(fasta_file, min_lines, padding):
    minimum_lines_per_file = min_lines
    # keeps track of previous and current sequence to compare duplicates
    previous_seq = ""
    current_seq = ""
    # stores a cluster of duplicates for writing to file
    write_buffer = ""
    # keep track of the header of current_seq to add to write_buffer
    last_header = ""

    split_index = 1
    file_to_write = open(fasta_file + ".split." + str(split_index).zfill(padding), 'w')

    with open(fasta_file) as file_to_read:
        line_count = 0
        for line in file_to_read:
            if line[0] == '>':
                if previous_seq != current_seq:
                    # When minimum lines exceeded open a new output file
                    if line_count > minimum_lines_per_file:
                        split_index += 1
                        line_count = len(write_buffer.splitlines())
                        file_to_write.close()
                        file_to_write = open(fasta_file + ".split." + str(split_index).zfill(padding), 'w')

                    file_to_write.write(write_buffer)
                    previous_seq = current_seq

                write_buffer = ""
                write_buffer += last_header + current_seq
                current_seq = ""
                last_header = line
            else:
                current_seq += line

            line_count += 1

        # Write remaining sequence once end of input file is reached
        write_buffer += last_header + current_seq
        file_to_write.write(write_buffer)
        file_to_write.close()


args = parse_args()
split_fasta_by_cluster(args.input_fasta, int(args.min_lines), int(args.padding))
