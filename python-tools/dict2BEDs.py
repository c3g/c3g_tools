#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
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
if sys.version_info[0] < 3:
    import httplib
else:
    import http.client

log = logging.getLogger(__name__)

class Sequence:
    def __init__(self, name, size):
        self._name = name
        self._size = size

    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, value):
        self._name = value

    @property
    def size(self):
        return self._size
    @size.setter
    def size(self, value):
        self._size = value

class BedInterval:
    def __init__(self, name, start, end):
        self._name = name
        self._start = start
        self._end = end

    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, value):
        self._name = value

    @property
    def start(self):
        return self._start
    @start.setter
    def start(self, value):
        self._size = value

    @property
    def end(self):
        return self._end
    @end.setter
    def end(self, value):
        self._end = value

def main():
    # Parse options
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-d", "--dict", help="Sequence Dictionary", type=argparse.FileType('r'), required=True)
    parser.add_argument("-i", "--interval", help="Interval File (optional ; parse interval instead)", type=bool, default=0)
    parser.add_argument("-b", "--beds", help="Output BED files", nargs="+", required=True)
    parser.add_argument("-c", "--chunk", help="chunk size in bp (optional ; Default no chunk)", type=int, default=0)
    parser.add_argument("-o", "--overlap", help="chunk overlap in bp (optional ; Default no overlap)", type=int, default=0)
    parser.add_argument("-r", "--remove", help="remove alt contigs (optional ; Default all chromsomes + alts)", type=bool, default=0)
    args = parser.parse_args()

    ordered_dict = parse_dictionary(args.dict, args.remove, args.interval)

    if not args.interval:
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
    
    if args.interval:
        file = args.beds.pop(0)
     
        if os.path.exists(file):
            os.remove(file)
        
        for seq in ordered_dict:
            regions = find_overlapping_regions(seq, args.chunk, args.overlap)

            with open(file, "a") as bed:
                for tuple in regions:
                    print("{}\t{}\t{}".format(*tuple), file=bed)
            

def parse_dictionary(dict, remove, interval):
    ordered_dict=[]
    #with open(dict) as dictFile:
    with dict as dictFile:
        startedSQ = False
        for line in dictFile:
            if not interval:
                if line.startswith('@SQ'):
                    match = re.search('^\@SQ.+SN:([^\t]+)\t.*LN:(\d+)\t.*', line)
                    if remove:
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
            else:
                if not line.startswith('@'):
                    match = line.split('\t')
                    if remove:
                        if "_" in match[0] or "." in match[0]:
                            continue
                        else:
                            seq = BedInterval(match[0], int(match[1]), int(match[2]))
                            ordered_dict.append(seq)
                    else:
                        seq = BedInterval(match[0], int(match[1]), int(match[2]))
                        ordered_dict.append(seq)
                elif startedSQ:
                    break
    #ordered_dict.sort(key=lambda x: x.size, reverse=True)
    return ordered_dict

def find_overlapping_regions(bedInterval, size, overlap):
    """
    Finds all overlapping regions within a given interval.
    :param interval: tuple representing the interval, e.g. (chr1, 0, 100)
    :param size: size of the regions
    :param overlap: overlap between adjacent regions
    :return: list of tuples representing the overlapping regions
    """
    start = bedInterval.start
    regions = []
    if not start == bedInterval.end:
        while start < bedInterval.end - size:
            regions.append((bedInterval.name, start, start + size))
            start += size - overlap
            
        regions.append((bedInterval.name, start, bedInterval.end))
    return regions

if __name__ == '__main__':
    main()

