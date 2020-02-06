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
    arg_parser.add_argument('--file-one', help='file containing rRNA IDs from readset one')
    arg_parser.add_argument('--file-two', help='file containing rRNA IDs from readset two')
    arg_parser.add_argument('--output-file', help='name of file with combined ids') 

    return arg_parser.parse_args()


def strip_id_of_paired_end(input_file, output_file, paired_end):
    """
        Create a single file containing all ids of rrna sequences found in both end. Remove read id.
    """
    write_mode = 'w' if paired_end == '/1' else 'a' 
    file_to_write = open(output_file, write_mode)

    with open(input_file) as f:
        for line in f:
            if line.startswith("@") or not line:
                continue

            rrna_id = line
            if line[-3:-1] == paired_end:
                rrna_id = line[:-3] + '\n'

            file_to_write.write(rrna_id)


args = parse_args()

strip_id_of_paired_end(args.file_one, args.output_file, '/1')
strip_id_of_paired_end(args.file_two, args.output_file, '/2')
