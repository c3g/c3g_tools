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
    args = parser.parse_args()

    ordered_dict = parse_dictionary(args.dict)

    total_size = 0
    for seq in ordered_dict:
        total_size += seq.size

    currentSize = 0
    targetSize = total_size / len(args.beds)
    currentBED = open(args.beds.pop(0), "w")
    for seq in ordered_dict:
        if len(args.beds) != 0 and currentSize != 0 and currentSize+seq.size >= targetSize:
            currentBED.close()
            print('File ' + currentBED.name + ' ' + str(currentSize))
            currentBED = open(args.beds.pop(0), "w")
            if len(args.beds) > 0:
                targetSize = total_size / len(args.beds)
            currentSize = 0
        currentBED.write(seq.name + '\t0\t' + str(seq.size) + '\n')
        total_size -= seq.size
        currentSize += seq.size
    currentBED.close()
    print('File2 ' + currentBED.name + ' ' + str(currentSize))

def parse_dictionary(dict):
    ordered_dict=[]
    #with open(dict) as dictFile:
    with dict as dictFile:
        startedSQ = False
        for line in dictFile:
            if line.startswith('@SQ'):
                match = re.search('^\@SQ.+SN:([^\t]+)\t.*LN:(\d+)\t.*', line)
                seq = Sequence(match.group(1), int(match.group(2)))
                ordered_dict.append(seq)
            elif startedSQ:
                break
    #ordered_dict.sort(key=lambda x: x.size, reverse=True)
    return ordered_dict

if __name__ == '__main__':
    main()

