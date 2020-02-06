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

"""
Split fasta or fastq into 2 output files, based on a JSON file with IDs
"""
import argparse
import textwrap


# Parse args
def parse_args():
    arg_parser = argparse.ArgumentParser(description=__doc__)

    original_data = arg_parser.add_mutually_exclusive_group()
    original_data.add_argument('--fastq')
    original_data.add_argument('--fasta')

    arg_parser.add_argument('--id-file')

    arg_parser.add_argument('--included', help='Output filename for reads that are in the IDs')
    arg_parser.add_argument('--excluded', help='Output filename for reads that are not in the IDs')
    arg_parser.add_argument('--out-format', help='Output format; fasta or fastq')
    arg_parser.add_argument('--readset-name', help='Name of readset')

    return arg_parser.parse_args()


def get_ids(id_file):
    """
    Get a set of ids from a JSON file

    Example JSON file:
    {
        "rows": [
                    {"id": "@SRR5"},
                    {"id": "@SRR7"}
                ]
    }

    Example output:
    {'@SRR5', '@SRR7'}

    :param id_file: JSON filename
    :return: set of str
    """
    ids = {}
    with open(id_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("@") or not line:
                continue

            _id = line.split('\t')[0]
            if '/' in _id:
                _id = _id.split('/')[0]
            ids[_id] = 1

    return ids


def partition_reads(all_reads, ids, in_file, out_file, out_format):
    """
    Write all_reads into 2 files based on where the read's id is in ids

    :param all_reads: dictionary of all reads
    :param ids: set of IDs
    :param in_file: the file to write reads that are in ids to
    :param out_file: the file to write reads that aren't in ids to
    :param out_format: fastq or fasta
    """
    included = open(in_file, 'w')
    excluded = open(out_file, 'w')
    
    for read in all_reads:
        data = all_reads[read]

        if out_format == 'fastq':
            if read.strip().split('/')[0] in ids:
                included.write('@' + data[0] + '\n' + data[1] + '\n+\n' + data[2] + '\n')
            else:
                excluded.write('@' + data[0] + '\n' + data[1] + '\n+\n' + data[2] + '\n')
        else:
            if read.strip().split('/')[0] in ids:
                included.write('>' + data[0] + '\n' + textwrap.fill(data[1], width=60) + '\n')
            else:
                excluded.write('>' + data[0] + '\n' + textwrap.fill(data[1], width=60) + '\n')


def get_all_reads(input_file, in_format, readset_name):
    seqs = {}

    # stores the last read seen
    last_read = ""

    # tells whether or not it's quality or sequence line
    is_sequence_line = True
 
    with open(input_file) as f:
        if in_format == 'fastq':
            # fastq
            # stores sequence, +, and quality data associated with last_read
            seq_info = ["", "", ""]

            for line in f:
                if line.startswith('@' + readset_name):
                    if last_read:
                        seq_info[0] = last_read
                        seqs[last_read] = seq_info
                        seq_info = ["", "", ""]

                    last_read = (line[1:]).split('\n')[0]
                    is_sequence_line = True

                elif line.startswith('+') and is_sequence_line:
                    is_sequence_line = False

                elif is_sequence_line:
                    # add sequence
                    seq_info[1] += line.strip()

                else:
                    # add quality data
                    seq_info[2] += line.strip()
 
            # add last sequence
            if last_read:
                seq_info[0] = last_read
                seqs[last_read] = seq_info

        else:
            # fasta
            # stores sequence header and sequence
            seq_info = ["", ""]

            for line in f:
                if line.startswith('>'):
                    if last_read:
                        seq_info[0] = last_read
                        seqs[last_read] = seq_info
                        seq_info = ["", ""]

                    last_read = (line[1:]).split('\n')[0]

                else:
                    seq_info[1] += line.strip()

            # add last sequence
            if last_read:
                seq_info[0] = last_read
                seqs[last_read] = seq_info

    return seqs


args = parse_args()
input_file, in_format = (args.fastq, 'fastq') if args.fastq else (args.fasta, 'fasta')

all_reads = get_all_reads(input_file, in_format, args.readset_name)
partition_reads(all_reads, get_ids(args.id_file), args.included, args.excluded, args.out_format)
