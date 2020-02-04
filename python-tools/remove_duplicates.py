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
Remove duplicates from a fastq file

Record the IDs, cluster size, and the length of each sequence in a *.json file
"""
import re
import textwrap
from argparse import ArgumentParser


def parse_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--input-fastq', help='Original fastq before de-duplication')
    arg_parser.add_argument('--unique-fasta', help='Fasta that has been de-duplicated')
    arg_parser.add_argument('--unique-uc', help='UC output from duplicate clustering step')

    arg_parser.add_argument('--output-read-description', help='File to store the IDs of reads that have been removed')
    arg_parser.add_argument('--output-fastq', help='Fastq that has been de-duplicated')
    arg_parser.add_argument('--output-fasta', help='Fasta that has been de-duplicated')
    arg_parser.add_argument('--readset-name', help='Name of readset')

    return arg_parser.parse_args()


def get_unique_ids(uc_file):
    """
    Get the fastq IDs of reads that are unique
    """
    def get_id(line):
        return re.match('^(\S+\t){8}(\S+)', line).group(2)

    # Lines that start with 'S' are not duplicates
    return {get_id(line) for line in open(uc_file) if line.startswith('S')}


def get_cluster_id_to_size(fasta):
    """
    ID of a cluster -> # of reads in cluster

    :param fasta: fasta that has been clustered with usearch
    :return: dict
    """
    def get_id(line):
        return re.search('Cluster(\d+);', line).group(1)

    def get_size(line):
        return int(re.search('size=(\d+)', line).group(1))

    return {get_id(line): get_size(line) for line in open(fasta) if line.startswith('>')}


def get_fastq_id_to_cluster_id(uc_file):
    """
    Fastq ID -> ID of cluster that read is part of

    :param uc_file: *.uc from usearch
    :return: dict
    """
    def fastq_id(line):
        return re.match('^(\S+\t){8}(\S+)', line).group(2)

    def cluster_id(line):
        return re.match('^(\S+\t){1}(\S+)', line).group(2)

    return {fastq_id(line): cluster_id(line) for line in open(uc_file) if line.startswith('S')}


def get_unique_reads(fastq, uc_file, readset_name):
    """
    Subset the reads from fastq to remove duplicates

    :param fastq: original fastq with duplicates
    :param uc_file: *.uc file
    :return: list
    """
    unique_ids = get_unique_ids(uc_file)
    all_reads = get_all_reads(fastq, readset_name)

    return [all_reads[read] for read in all_reads if read in unique_ids]


def assign_cluster_size(unique_reads, uc_file, unique_fasta):
    """
    Set the 'cluster_size' attribute on all unique reads

    :param unique_reads: iterable of reads
    :param uc_file: *.uc file
    :param unique_fasta: output from usearch that contains the number of duplicates for each read
    """
    fastq_id_to_cluster_id = get_fastq_id_to_cluster_id(uc_file)
    cluster_id_to_size = get_cluster_id_to_size(unique_fasta)

    for read in unique_reads:
        read[3] = cluster_id_to_size[fastq_id_to_cluster_id[read[0]]]

def write_id_to_cluster_size(reads, id_file):
    """
    Write a table in json containing fastq ID and cluster size
    Eg:
    {"rows": [
        {"id": "SRR1", "cluster_size": 3}
        {"id": "SRR23", "cluster_size": 1}
    ]}
    """
    with open(id_file, 'w+') as f:
        # Write header
        f.write("@id\tcluster_size\tlength\n")
        for read in reads:
            f.write("{_id}\t{cluster_size}\t{length}\n".format(_id=read[0],
                                                               cluster_size=read[3],
                                                               length=len(read[1])))


def write_unique_reads(reads, out_file, out_format):
    """
    Write reads to file according to formatt

    :param reads: iterable of reads
    :param out_file: *.fastq file or *.fasta
    :param out_format: 'fastq' or 'fasta'
    :return:
    """
    with open(out_file, 'w+') as out:
        if out_format == 'fastq':
            for data in reads:
                out.write('@' + data[0] + '\n' + data[1] + '\n+\n' + data[2] + '\n')
        else:
            for data in reads:
                out.write('>' + data[0] + '\n' + textwrap.fill(data[1], width=60) + '\n')


def get_all_reads(input_file, readset_name):
    seqs = {}
 
    # stores the last read seen
    last_read = ""
    
    # tells whether or not it's quality or sequence line
    is_sequence_line = True

    with open(input_file) as f:
        # stores sequence header, sequence, quality data, and cluster size
        seq_info = ["", "", "", 0]
        
        for line in f:
            if line.startswith('@' + readset_name):
                if last_read:
                    seq_info[0] = last_read
                    seqs[last_read] = seq_info
                    seq_info = ["", "", "", 0]

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

    return seqs


args = parse_args()

# Get the unique reads and determine how many there are in each cluster
unique_reads = get_unique_reads(args.input_fastq, args.unique_uc, args.readset_name)
assign_cluster_size(unique_reads, args.unique_uc, args.unique_fasta)

# Keep track of how many duplicates there are of each read
write_id_to_cluster_size(unique_reads, args.output_read_description)

# Write out a *.fasta and *.fastq file containing only the unique reads
write_unique_reads(unique_reads, args.output_fastq, 'fastq')
write_unique_reads(unique_reads, args.output_fasta, 'fasta')
