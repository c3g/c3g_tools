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

import argparse
import csv
import glob
import httplib
import logging
import os
import re

log = logging.getLogger(__name__)

class Sequence:
    def __init__(self, name, size):
        self.name = name
        self.size = size

    @property
    def name():
        return self._name
    @property
    def size():
        return self._size

def main():
    # Parse options
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-d", "--dict", help="Sequence Dictionary", type=file, required=True)
    parser.add_argument("-b", "--beds", help="Output BED files", nargs="+", required=True)
    parser.add_argument("-c", "--chunk", help="chunk size in bp (optional ; Default no chunk)", type=int, default=0)
    parser.add_argument("-o", "--overlap", help="chunk overlap in bp (optional ; Default no overlap)", type=int, default=0)
    parser.add_argument("-r", "--remove", help="remove alt contigs (optional ; Default all chromsomes + alts)", type=bool, default=0)
    args = parser.parse_args()

    ordered_dict = parse_dictionary(args.dict)

    total_size = 0
    for seq in ordered_dict:
        total_size += seq.size

    currentSize = 0
    if len(args.beds) > 0:
        targetSize = total_size / len(args.beds)
    else :
        targetSize = total_size 
    currentBED = open(args.beds.pop(0), "w")
    for seq in ordered_dict:
        print targetSize
        if len(args.beds) != 0 and currentSize != 0 and currentSize+seq.size >= targetSize:
            currentBED.close()
            print('File ' + currentBED.name + ' ' + str(currentSize))
            currentBED = open(args.beds.pop(0), "w")
            currentSize = 0
        if args.chunk > 0 :
            start = 0
            end = start + args.chunk
            while end < seq.size:
                currentBED.write(seq.name + '\t'+str(start)+'\t' + str(end) + '\n')
                start += (args.chunk - args.overlap)
                end = start + args.chunk
            currentBED.write(seq.name + '\t'+str(start)+'\t' + str(seq.size) + '\n')
        else :
            currentBED.write(seq.name + '\t0\t' + str(seq.size) + '\n')
        total_size -= seq.size
        currentSize += seq.size
    currentBED.close()
    print('File ' + currentBED.name + ' ' + str(currentSize))

def parse_dictionary(dict):
    ordered_dict=[]
    #with open(dict) as dictFile:
    with dict as dictFile:
        startedSQ = False
        for line in dictFile:
            if line.startswith('@SQ'):
                match = re.search('^\@SQ.+SN:([^\t]+)\t.*LN:(\d+)\t.*', line)
                if args.remove:
                    if "_" in match.group(1) or "." in match.group(1):
                        continue
                    else:
                        seq = Sequence(match.group(1), int(match.group(2)))
                        ordered_dict.append(seq)
                else:
                    seq = Sequence(match.group(1), int(match.group(2)))
                    ordered_dict.append(seq)
            elif startedSQ:
                break
    #ordered_dict.sort(key=lambda x: x.size, reverse=True)
    return ordered_dict

if __name__ == '__main__':
    main()

