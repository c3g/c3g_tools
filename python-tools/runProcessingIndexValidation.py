#!/usr/bin/env python

### Edouard Henrion (2020/03/18) - edouard.henrion@computationalgenomics.ca

#from __future__ import unicode_literals
#from __future__ import print_function, division, unicode_literals, absolute_import
import os
import sys
import getopt
import re
import json
import csv
import ast
from itertools import imap
from operator import ne
#from numarray import sum, arraoy
from bs4 import BeautifulSoup

def getarg(argument):
        index_json=""
        index_file=""
        html_report=""
        lane=""
        mismatches=""
        optli,arg = getopt.getopt(argument[1:],"j:f:r:l:m:h",['index_json','index_file','html_report','lane','mismatches','help'])
        if len(optli) == 0 :
            usage()
            sys.exit("Error : No argument given")
        for option, value in optli:
            if option in ("-j","--index_json"):
                index_json=str(value)
            if option in ("-f","--index_file"):
                index_file=str(value)
            if option in ("-r","--html_report"):
                html_report=str(value)
            if option in ("-l","--lane"):
                lane=int(value)
            if option in ("-m","--mismatches"):
                mismatches=int(value)
            if option in ("-h","--help"):
                usage()
                sys.exit()

        if not os.path.exists(index_json):
            sys.exit("Error - JSON index file not found:\n" + str(index_json))
        if not os.path.exists(index_file):
            sys.exit("Error - index file not found:\n" + str(index_file))
        if not os.path.exists(html_report):
            sys.exit("Error - BCL2fastq html report not found:\n" + str(html_report))
        if lane == "":
            sys.exit("Error - missing lane argument -l ...\n" + usage())
        if mismatches == "":
            sys.exit("Error - missing argument -m...\n" + usage())

        return index_json, index_file, html_report, lane, mismatches

#def numeric_hamming_distance(num1, num2):
#    assert len(num1) == len(num2)
#    return sum(num1 != num2)

def distance(
    str1,
    str2
    ):
    """
    Returns the hamming distance. http://code.activestate.com/recipes/499304-hamming-distance/#c2
    """
    return sum(imap(ne, str1, str2))

def parse_bcl2fastq_html_report(html_file):

    with open(html_file, 'r') as file:
        html_string = file.read()

    soup = BeautifulSoup(html_string, 'lxml') # Parse the HTML as a string
    table = soup.find_all('table')[1] # Grab the second table

    rows = table.find_all('tr')
    header = rows[0].find_all('th') # Parse the header row
    values = rows[1].find_all('td') # Parse the value row

    # For now we only parse the # raw reads value
    return int(values[0].get_text().replace(",", ""))

def validate(
    index_json,
    index_file,
    bcl2fastq_html_report,
    lane,
    mismatches,
    ):

    with open(index_json) as json_file:
        index_per_readset = json.load(json_file)

    # First identify all the expected indexes as well as all the unfound indexes
    print "index_per_readset.items() " + str(len(index_per_readset.items()))
    readset_count = 0
    expected = 0
    indexes_to_validate = []
    for readset_name, indexes in index_per_readset.items():
        readset_count += 1
        print readset_count
        # For the current readset, retrieve all the possible indexes (from LIMS)
        # and store them as tuples along their index_name in the indexes_to_test list
        if len(indexes) > 1 :
            for index in indexes:
                expected += 1
                if index['INDEX1'] != "" and index['INDEX2'] != "":
                    indexes_to_validate.append((index['INDEX1']+index['INDEX2'], index['INDEX_NAME'], index['LIBRARY']))
                else:
                    if index['INDEX1'] != "": indexes_to_validate.append((index['INDEX1'], index['INDEX_NAME'], index['LIBRARY']))
                    if index['INDEX2'] != "": indexes_to_validate.append((index['INDEX2'], index['INDEX_NAME'], index['LIBRARY']))
        else:
            index = indexes[0]
            expected += 1
            if index['INDEX1'] != "" and index['INDEX2'] != "":
                indexes_to_validate.append((index['INDEX1']+index['INDEX2'], index['INDEX_NAME'], index['BCL2FASTQ_NAME']))
            else:
                if index['INDEX1'] != "": indexes_to_validate.append((index['INDEX1'], index['INDEX_NAME'], index['LIBRARY']))
                if index['INDEX2'] != "": indexes_to_validate.append((index['INDEX2'], index['INDEX_NAME'], index['LIBRARY']))

    pf_read_count = 0
    index_cluster_count = 0
    found_index = {}
    found_barcode = {}
    unfound_index = {}
    unexpected = {}
    # Parse the index file (i.e. output from `index` job`), skipping the comment lines
    index_csv = csv.DictReader(filter(lambda row: row != '\n' and row[0]!='#', open(index_file, 'rb')), delimiter='\t', quotechar='"')
    for row in index_csv:

        # Skip row if no passed filter reads
        if int(row['PF_READS']) == 0:
            continue

        # Increment the count of PF (Passed Filter) reads
        pf_read_count += int(row['PF_READS'])

        for (idx_seq, idx_name, library) in indexes_to_validate:

            mism_count = distance(row['BARCODE'], idx_seq)

            # If current sequences match with at most the number of mismatches
            if mism_count < mismatches + 1:

                # Consider current index as found and set or update its # PF_READS
                if idx_name in found_index:
                   found_index[idx_name]['pf_reads'] += int(row['PF_READS'])
                   if mism_count == 0 : found_index[idx_name]['0_mis'] += int(row['PF_READS'])
                   if mism_count == 1 : found_index[idx_name]['1_mis'] += int(row['PF_READS'])
                else:
                   found_index[idx_name] = {
                       'pf_reads': int(row['PF_READS']),
                       '0_mis': int(row['PF_READS']) if mism_count == 0 else 0,
                       '1_mis': int(row['PF_READS']) if mism_count == 1 else 0,
                       'idx_seq': idx_seq,
                       'library': library
                   }

                index_cluster_count += int(row['PF_READS'])

                # remove current row from unexpected list if it was already put in
                if (row['BARCODE'] in unexpected):
                    del unexpected[row['BARCODE']]
                # set the current 'BARCODE' as found
                found_barcode[row['BARCODE']] = True

            # If there is no match but # passed filter reads > 0 AND current 'BARCODE' have not been found yet for any other readset
            elif row['PF_READS'] > 0 and row['BARCODE'] not in found_barcode:
                # push current row to the unexpected list
                unexpected[row['BARCODE']] = {
                   'pf_reads': int(row['PF_READS']),
                   'idx_name': row['BARCODE_NAMES'],
                   'idx_seq': row['BARCODE']
                }

    # Once whole index file has been read,
    # check if all indexes of current readset were actually found
    for (idx_seq, idx_name, bcl2fastq) in indexes_to_validate:
        # if not, push index to unfound list
        if idx_name not in found_index:
            unfound_index[idx_name] = idx_seq

    # Compute spread
    pf_read_count_list = [x['pf_reads'] for x in found_index.values()]
    max_pf_reads = max(pf_read_count_list)
    print max_pf_reads
    min_pf_reads = min(pf_read_count_list)
    print min_pf_reads
    spread = max_pf_reads / float(min_pf_reads)
    print spread

    # Count the number of undetermined index reads
    unexpected_read_count = 0
    unexpected_index = 0
    for idx_name, idx_dict in unexpected.items():
        unexpected_read_count += int(idx_dict['pf_reads'])
        if idx_dict['idx_name'] != "": unexpected_index += 1
    print unexpected_read_count

    # % of undetermined index reads
    unexpected_ratio = float(unexpected_read_count) / pf_read_count
    print pf_read_count
    print unexpected_ratio

    # select the 100 most represented unexpeced, to show then in the report
    unexpected_to_show = sorted(unexpected.items(), key=lambda x: x[1]['pf_reads'], reverse=True)[:100]

    # % expected indices
    expected_ratio = float(len(found_index.keys())) / expected

    # Count the number of raw reads
    raw_read_count = parse_bcl2fastq_html_report(bcl2fastq_html_report)

    # Now build the hash table which will be printed into a JSON report file
    index_validation_hash = {
        'raw_reads': raw_read_count,
        'pf_reads': pf_read_count,
        'index_clusters': index_cluster_count,
        '%_pf_reads': float(pf_read_count) / raw_read_count,
        'spread' : spread,
        '%_undetermined_idx_reads' : unexpected_ratio * 100,
        '%_expected_indices' : expected_ratio * 100,
        'expected' : [{
            "pf_reads" : idx_dict['pf_reads'],
            "%total" : float(idx_dict['pf_reads'])*100 / pf_read_count,
            "%index_clusters_in_lane" : float(idx_dict['pf_reads'])*100 / index_cluster_count,
            "0_mismatch" : float(idx_dict['0_mis'])*100 / idx_dict['pf_reads'],
            "1_mismatch" : float(idx_dict['1_mis'])*100 / idx_dict['pf_reads'],
            "index" : idx_dict['idx_seq'],
            "name" : idx_name,
            "library_name" : idx_dict['library']
        } for idx_name, idx_dict in found_index.items()],
        'not_found' : [{
            "index" : idx_seq,
            "name" : idx_name
        } for idx_name, idx_seq in unfound_index.items()],
        'unexpected' : [{
            "pf_reads" : idx_dict['pf_reads'],
            "%total" : float(idx_dict['pf_reads'])*100 / pf_read_count,
            "index" : idx_dict['idx_seq'],
            "name" : idx_dict['idx_name']
        } for idx_name, idx_dict in unexpected_to_show],
    }

    index_validation_str = json.dumps(index_validation_hash, indent=4)

    # Print to file
    filepath = os.path.join(os.path.dirname(index_file), "index_validation_report.json")
    with open(filepath, 'w') as out_json:
        out_json.write(index_validation_str)

def usage():
        print "USAGE : runProcessingIndexValidation.py"
        print "       -j :        lims_index_json_file"
        print "       -f :        index_file"
        print "       -r :        BCL2fastq lane html report file (usually in <BCL2fastq_OUTPUT_DIR>/Reports/html/<FLOWCELL_ID>/all/all/all/lane.html)"
        print "       -l :        lane"
        print "       -m :        mismatches"
        print "       -h :        this help\n"


def main():
        print "\n-------------------------------------------------------------------------------------"
        print "runProcessingIndexValidation.py uses indexes from Clarity LIMS, passed in hash string,"
        print "to validate the indexes produced by CountIlluminaBarcodes and"
        print "output a report in JSON format."
        print "This program was written by Edouard Henrion"
        print "For more information, contact: edouard.henrion@computationalgenomics.ca"
        print "----------------------------------------------------------------------------------\n"
        index_hash, index_file, html_report, lane, mismatches = getarg(sys.argv)
        validate(index_hash, index_file, html_report, lane, mismatches)

main()

