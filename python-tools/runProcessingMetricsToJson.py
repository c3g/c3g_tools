#!/usr/bin/env python

### Edouard Henrion (2020/03/18) - edouard.henrion@computationalgenomics.ca

import errno
import os
import sys
import getopt
import re
import json
import csv
import subprocess
import time
import random
import shutil
import signal
from collections import namedtuple
import xml.etree.ElementTree as ET

def getarg(argument):
        json_file=""
        step=""
        readset=None
        platform=""
        inputs=[]
        optli,arg = getopt.getopt(argument[1:],"j:s:r:p:i:h",['json_file','step','readset','platform','input','help'])
        if len(optli) == 0 :
            usage()
            sys.exit("Error : No argument given")
        for option, value in optli:
            if option in ("-j","--json_file"):
                json_file=str(value)
            if option in ("-s","--step"):
                step=str(value)
            if option in ("-r","--readset"):
                readset=str(value)
            if option in ("-p","--platform"):
                platform=str(value)
            if option in ("-i","--input"):
                inputs.append(str(value))
            if option in ("-h","--help"):
                usage()
                sys.exit()

        if not json_file:
            sys.exit("Error - json_file parameter is required")
        elif json_file and not os.path.exists(json_file):
            sys.exit("Error - JSON index stats file not found:\n" + str(json_file))
        if not step:
            sys.exit("Error - step parameter is required")
        if not readset and step not in ['basecall', 'index', 'fastq']:
            sys.exit("Error - readset parameter is requeired for step " + step)
        if not platform:
            sys.exit("Error - platform parameter is required")
        if len(inputs) == 0:
            sys.exit("Error - input parameter is required")
        else:
            not_found = []
            for in_file in inputs:
                if not os.path.exists(json_file):
                    not_found.append(in_file)
            if not_found:
                sys.exit("Error - input file(s) not found :\n  " + "\n  ".join(not_found))

        return json_file, step, platform, inputs, readset

def lock(filepath):
    unlocked = True
    while unlocked:
        try:
            os.makedirs(filepath + '.lock')
        except OSError as exception:
            if exception.errno == errno.EEXIST and os.path.isdir(filepath + '.lock'):
                # The lock folder already exists, we need to wait for it to be deleted
                sleep_time = random.randint(1, 100)
                time.sleep(sleep_time)
                pass
            else:
                # An unexpected error has occured : let's stop the program and raise the error"
                raise exception
        else:
            # The lock folder was successfully created !"
            unlocked = False

def unlock(filepath):
    shutil.rmtree(filepath + '.lock', ignore_errors=True)

def parseMetricsFile(
    metrics_file
    ):
    reader = open(metrics_file, 'r').readlines()
    metrics_tsv = []
    header = []
    for row in reader:
        if row!="\n" and not str(row).startswith('##') and not str(row).startswith('# '):
            if not header:
                header = row.strip().split("\t")
                header = [head.strip() for head in header]
            else:
                metrics_tsv.append(dict(zip(header, row.strip().split("\t"))))
    return metrics_tsv

def getIndexHash_from_splitBarcode(
    input_file,
    barcode_name,
    barcode_sequences,
    dict_to_update
    ):

    sequence_stat_file = ""
    if "SequenceStat.txt" in input_file:
        sequence_stat_file = input_file
    else:
        sys.exit("Error - Unexpected input file found for splitBarcode index metrics : " + input_file)

    sequence_stat_tsv = parseMetricsFile(sequence_stat_file)
    #Sequence        Barcode         Count   Percentage(%)
    total_pf = sum([int(row['Count']) for row in sequence_stat_tsv])
    total_pf_onindex_inlane = sum([int(row['Count']) for row in sequence_stat_tsv if row['Barcode'] != "undecoded"])

    pf_clusters = sum([int(row['Count']) for row in sequence_stat_tsv if barcode_name in row['Barcode']])
    pf_perfect_matches = sum([int(row['Count']) for row in sequence_stat_tsv if barcode_name in row['Barcode'] and row['#Sequence'] in barcode_sequences])
    pf_one_mismatch_matches = sum([int(row['Count']) for row in sequence_stat_tsv if barcode_name in row['Barcode'] and not row['#Sequence'] in barcode_sequences])

    dict_to_update['pf_clusters'] = pf_clusters
    dict_to_update['pct_of_the_lane'] = 100*pf_clusters/float(total_pf)
    dict_to_update['pct_on_index_in_lane'] = 100*pf_clusters/float(total_pf_onindex_inlane)
    dict_to_update['pct_perfect_barcode'] = 100*pf_perfect_matches/float(pf_clusters)
    dict_to_update['pct_one_mismatch_barcode'] = 100*pf_one_mismatch_matches/float(pf_clusters)

    return dict_to_update

def getIndexHash_from_BCL2fastq(
    input_file,
    readset,
    dict_to_update
    ):

    stats_json_file = ""
    if ".index_stats.json" in input_file:
        stats_json_file = input_file
    else:
        sys.exit("Error - Unexpected input file found for bcl2fasfq index metrics : " + input_file)

    with open(stats_json_file, 'r') as json_file:
        stats_json = json.load(json_file)

    total_raw = stats_json['ConversionResults'][0]['TotalClustersRaw']
    total_pf = stats_json['ConversionResults'][0]['TotalClustersPF']
    total_pf_onindex_inlane = sum([sample['NumberReads'] for sample in stats_json['ConversionResults'][0]['DemuxResults']])

    for i, sample in enumerate(stats_json['ConversionResults'][0]['DemuxResults']):
        if sample['SampleName'] == readset:
            dict_to_update['pf_clusters'] = sample['NumberReads'],
            dict_to_update['pct_of_the_lane'] = 100*sample['NumberReads']/float(total_pf),
            dict_to_update['pct_on_index_in_lane'] = 100*sample['NumberReads']/float(total_pf_onindex_inlane),
            dict_to_update['pct_perfect_barcode'] = 100*sample['IndexMetrics'][0]['MismatchCounts']['0']/float(sample['NumberReads']),
            dict_to_update['pct_one_mismatch_barcode'] = 100*sample['IndexMetrics'][0]['MismatchCounts']['1']/float(sample['NumberReads']),
            dict_to_update['yield'] = sample['Yield'],
            dict_to_update['pct_q30_bases'] = 100*sum([readMetrics['YieldQ30'] for readMetrics in sample['ReadMetrics']])/float(sample['Yield']),
            dict_to_update['mean_quality_score'] = sum([readMetrics['QualityScoreSum'] for readMetrics in sample['ReadMetrics']])/float(sample['Yield'])

    return dict_to_update

def getIndexHash_from_DemuxFastqs(
    input_file,
    readset,
    dict_to_update
    ):

    if ".DemuxFastqs.metrics.txt" in input_file:
        demux_metrics_file = input_file
    else:
        sys.exit("Error - Unexpected input file found for demuxfastq index metrics : " + input_file)

    demux_metrics_tsv = parseMetricsFile(demux_metrics_file)

    total_pf = sum([int(row['pf_templates']) for row in demux_metrics_tsv])
    total_pf_onindex_inlane = sum([int(row['pf_templates']) for row in demux_metrics_tsv if row['library_name'] != "unmatched"])

    # retrieve library from the readset
    readset_libraries = list(set([row['library_name'] for row in demux_metrics_tsv if readset in row['barcode_name']]))
    if not len(readset_libraries) == 1:
        sys.exit("Error - More than one library identified for readset" + readset + " in metrics file " + demux_metrics_file + " :\n  " + "\n  ".join(readset_libraries))
    else:
        library = readset_libraries[0]

    raw_barcode_keys = [row['barcode_name'] for row in demux_metrics_tsv if row['library_name'] == library]
    barcode_keys = []
    barcode_name = ""
    for name in raw_barcode_keys:
        m = re.search(readset+"_(\w_)?(?P<barcode_name>.+)", name)
        if m:
            if barcode_name and barcode_name != m.group('barcode_name'):
                print("Error !! At least more than one barcode (" + barcode_name + ", " + m.group('barcode_name') + ") was found to be associated to library(" + library + ")")
                sys.exit(32)
            barcode_name = m.group('barcode_name')
            barcode_keys.append(name)
        else:
            print(name + " could not be parsed for barcode...")

    pf_clusters = sum([int(row['pf_templates']) for row in demux_metrics_tsv if row['library_name'] == library and row['barcode_name'] in barcode_keys])
    pf_perfect_matches = sum([int(row['pf_perfect_matches']) for row in demux_metrics_tsv if row['library_name'] == library and row['barcode_name'] in barcode_keys])
    pf_one_mismatch_matches = sum([int(row['pf_one_mismatch_matches']) for row in demux_metrics_tsv if row['library_name'] == library and row['barcode_name'] in barcode_keys])

    dict_to_update['pf_clusters'] = pf_clusters
    dict_to_update['pct_of_the_lane'] = 100*pf_clusters/float(total_pf)
    dict_to_update['pct_on_index_in_lane'] = 100*pf_clusters/float(total_pf_onindex_inlane)
    dict_to_update['pct_perfect_barcode'] = 100*pf_perfect_matches/float(pf_clusters)
    dict_to_update['pct_one_mismatch_barcode'] = 100*pf_one_mismatch_matches/float(pf_clusters)

    return dict_to_update

def getIndexHash_from_FastQC(
    input_file,
    dict_to_update
    ):

    if "fastqc_data.txt" in input_file:
        fastqc_file = input_file
    else:
        sys.exit("Error - Unexpected input file found for fastqc index metrics : " + input_file)

    q30_reads = int(subprocess.check_output("sed -n '/#Quality/,/END_MODULE/p' %s | awk '{if($1!=\">>END_MODULE\" && $1!=\"#Quality\" && $1>=30){sum+=$2}} END {print sum}'" % fastqc_file, shell=True).strip())
    seq_length = int(subprocess.check_output("grep 'Sequence length' %s | cut -f 2" % fastqc_file, shell=True).strip())
    total_readset_seq = int(subprocess.check_output("grep 'Total Sequences' %s | cut -f 2" % fastqc_file, shell=True).strip())
    readset_yield = seq_length * float(total_readset_seq)
    q30_bases = seq_length * (q30_reads)
    pct_q30_bases = 100 * q30_bases / float(readset_yield)
    mean_quality_score = float(subprocess.check_output("sed -n '/#Base\tMean/,/END_MODULE/p' %s | awk '{if($1!=\">>END_MODULE\" && $1!=\"#Mean\"){records+=1;sum+=$2}} END {print sum/records}'" % fastqc_file, shell=True).strip())

    # At this point, 'PF Clusters' should be available in the JSON
    if dict_to_update['pf_clusters']:
        pf_clusters = dict_to_update['pf_clusters']
    # If not...
    else:
        pf_clusters = total_readset_seq

    dict_to_update['yield'] = seq_length*float(pf_clusters)
    dict_to_update['pct_q30_bases'] = pct_q30_bases
    dict_to_update['mean_quality_score'] = mean_quality_score

    return dict_to_update

def getQcHash(
    input_file,
    readset,
    dict_to_update
    ):

    pattern = re.compile("mpsQC_" + readset + "_[0-9]+_L00[0-9]_stats.xml$")
    if pattern.search(input_file):
        qc_graph_xml = input_file
    else:
        sys.exit("Error - Unexpected input file found for qc_graphs metrics : " + input_file)

    root = ET.parse(qc_graph_xml).getroot()

    dict_to_update['avg_qual'] = root.attrib['avgQual']
    dict_to_update['duplicate_rate'] = root.attrib['duplicateRate']
    dict_to_update['nb_reads'] = root.attrib['nbReads']
    dict_to_update['nb_bases'] = root.attrib['nbBases']

    return dict_to_update

def getBlastHash(
    input_file,
    readset,
    dict_to_update
    ):

    pattern = re.compile(readset + "_[0-9]+_L00[0-9].R1.RDP.blastHit_20MF_species.txt$")
    if pattern.search(input_file):
        blast_output = input_file
    else:
        sys.exit("Error - Unexpected input file found for blast metrics : " + input_file)

    BlastResult = namedtuple('BlastResult', 'hits species')

    with open(blast_output, 'r') as f:
        blast_result = [line.rstrip() for line in f]

    # best match
    m = re.search("(?P<hits>\w+) (?P<match>[\w\.]+)", blast_result[0].lstrip())
    best_match = None
    if m:
        best_match = BlastResult(m.group('hits'), m.group('match'))

    # 2nd hit
    nd_hit = None
    if len(blast_result) > 1:
        m = re.search("(?P<hits>\w+) (?P<match>[\w\.]+)", blast_result[1].lstrip())
        if m:
            nd_hit = BlastResult(m.group('hits'), m.group('match'))

    # 3rd hit
    rd_hit = None
    if len(blast_result) > 2:
        m = re.search("(?P<hits>\w+) (?P<match>[\w\.]+)", blast_result[2].lstrip())
        if m:
            rd_hit = BlastResult(m.group('hits'), m.group('match'))

    dict_to_update['1st_hit'] = best_match.species + " (" + best_match.hits + ")" if best_match else None
    dict_to_update['2nd_hit'] = nd_hit.species + " (" + nd_hit.hits + ")" if nd_hit else None
    dict_to_update['3rd_hit'] = rd_hit.species + " (" + rd_hit.hits + ")" if rd_hit else None

    return dict_to_update

def getSampleTagHash():
    # Parsing of the kapa tag result file still remains to be done
    return None

def getAlignmentHash_from_picardMarkDup(
    input_file,
    readset,
    dict_to_update
    ):

    dup_metrics_file = ""
    if ".sorted.dup.metrics" in input_file:
        dup_metrics_file = input_file
    else:
        sys.exit("Error - Unexpected input file found for picard_mark_duplicates metrics : " + readset + ".sorted.dup.metrics")

    dup_tsv = parseMetricsFile(dup_metrics_file)

    dict_to_update['aligned_dup_rate'] = dup_tsv[0]['PERCENT_DUPLICATION']

    return dict_to_update

def getAlignmentHash(
    inputs,
    readset,
    dict_to_update
    ):

    alignment_summary_metrics_file = ""
    insert_size_metrics_file = ""
    verify_bam_id_file = ""
    target_coverage_file = ""
    for in_file in inputs:
        if ".sorted.metrics.alignment_summary_metrics" in in_file:
            alignment_summary_metrics_file = in_file
        elif ".sorted.metrics.insert_size_metrics" in in_file:
            insert_size_metrics_file = in_file
        elif ".sorted.metrics.verifyBamId.tsv" in in_file:
            verify_bam_id_file = in_file
        elif ".sorted.metrics.verifyBamId.selfSM" in in_file:
            verify_bam_id_file = in_file
        elif ".sorted.metrics.targetCoverage.txt" in in_file:
            target_coverage_file = in_file
        elif ".metrics.rRNA.tsv" in in_file:
            # nothing iRNA specific has been implemented for now...
            continue
        elif ".rnaseqc.sorted.dup.metrics.tsv" in in_file:
            # nothing iRNA specific has been implemented for now...
            continue
        else:
            sys.exit("Error - Unexpected input file found for alignment metrics : " + in_file)
    not_found = []
    if not alignment_summary_metrics_file:
        not_found.append(readset + ".sorted.metrics.alignment_summary_metrics")
    if not insert_size_metrics_file:
         not_found.append(readset + ".sorted.metrics.insert_size_metrics")
    # if not verify_bam_id_file:
    #     not_found.append(readset + ".sorted.metrics.verifyBamId.tsv")
#    if not target_coverage_file:
        # don't do anythign as it could just be an RNA case...
        #not_found.append(readset + ".sorted.metrics.targetCoverage.txt")
    if not len(not_found) == 0:
        sys.exit("Error - alignement metrics file(s) not found :\n  " + "\n  ".join(not_found))

    align_tsv = parseMetricsFile(alignment_summary_metrics_file)
    insert_tsv = parseMetricsFile(insert_size_metrics_file)
    if os.path.isfile(verify_bam_id_file):
        verifyBamID_tsv = parseMetricsFile(verify_bam_id_file)
    else:
        verifyBamID_tsv = None

    if os.path.isfile(target_coverage_file):
        sex_match_reader = csv.DictReader(open(target_coverage_file, 'r'), delimiter='\t')

        chrX_cov, chrX_covered_bases = 0, 0
        chrY_cov, chrY_covered_bases = 0, 0
        total_cov = 0
        for row in sex_match_reader:
            if row['IntervalName'] == "chrX":
                chrX_cov += int(row['TotalCoverage'])
                chrX_covered_bases += int(row['TotalNbCoveredBases'])
            if row['IntervalName'] == "chrY":
                chrY_cov += int(row['TotalCoverage'])
                chrY_covered_bases += int(row['TotalNbCoveredBases'])
            if row['IntervalName'] == "Total":
                total_cov = row['MeanCoverage']

        chrX_cov = chrX_cov / float(chrX_covered_bases)
        chrY_cov = chrY_cov / float(chrY_covered_bases)

        if chrX_cov > 0.8:
            sex_det = "F"
        elif chrY_cov > 0.25:
            sex_det = "M"
        else:
            sex_det = "?"
        gender = dict_to_update['reported_sex'] # At this point, 'reported_sex' should be available in the JSON
        if sex_det == gender:
            sex_match = True
        else:
            sex_match = False
    else:
        total_cov = "N/A"
        sex_det = "N/A"
        sex_match = "N/A"

    freemix = verifyBamID_tsv[0]['FREEMIX'] if verifyBamID_tsv else 'N/A'

    dict_to_update['pf_read_alignment_rate'] = align_tsv[2]['PCT_PF_READS_ALIGNED']
    dict_to_update['mean_coverage'] = total_cov
    dict_to_update['chimeras'] = align_tsv[2]['PCT_CHIMERAS']
    dict_to_update['adapter_dimers'] = align_tsv[2]['PCT_ADAPTER']
    dict_to_update['average_aligned_insert_size'] = insert_tsv[0]['MEAN_INSERT_SIZE']
    dict_to_update['freemix'] = freemix
    dict_to_update['inferred_sex'] = sex_det
    dict_to_update['sex_concordance'] = sex_match

    return dict_to_update

def report(
    json_file,
    step,
    platform,
    inputs,
    readset=None
    ):

    with open(json_file, 'r') as json_fh:
        run_report_json = json.load(json_fh)

    report_version = run_report_json['version']
    readsets = [readset] if readset else [record['sample'] for record in run_report_json['run_validation']]
    for readset in readsets:
        for record in run_report_json['run_validation']:
            if record['sample'] == readset:

                if step in ['basecall', 'index', 'fastq']:
                        section = 'index'
                        if platform == 'mgit7':
                            # retrieve barcode name for the readset
                            if report_version == "1.0":
                                readset_barcodes = list(set([record['INDEX_NAME'] for record in run_report_json['barcodes'][readset]]))
                                barcode_sequences = [record['BARCODE_SEQUENCE'] for record in run_report_json['barcodes'][readset]]
                            else:
                                readset_barcodes = list(set([record['INDEX_NAME'] for record in run_report_json['readsets'][readset]['barcodes']]))
                                barcode_sequences = [record['BARCODE_SEQUENCE'] for record in run_report_json['readsets'][readset]['barcodes']]
                            if not len(readset_barcodes) == 1:
                                sys.exit("Error - More than one barcode identified for readset" + readset + " in json file " + json_file + " :\n  " + "\n  ".join(readset_barcodes))
                            else:
                                barcode_name = readset_barcodes[0]
                            new_dict = getIndexHash_from_splitBarcode(inputs[0], barcode_name, barcode_sequences, record[section])
                        if platform == 'illumina':
                            new_dict = getIndexHash_from_BCL2fastq(inputs[0], readset, record[section])
                        if platform == 'mgig400':
                            new_dict = getIndexHash_from_DemuxFastqs(inputs[0], readset, record[section])

                elif step == 'fastqc':
                    section = 'index'
                    new_dict = getIndexHash_from_FastQC(inputs[0], record[section])

                elif step == 'qc_graphs':
                    section = 'qc'
                    new_dict = getQcHash(inputs[0], readset, record[section])

                elif step == 'blast':
                    section = step
                    new_dict = getBlastHash(inputs[0], readset, record[section])

                elif step == 'picard_mark_duplicates':
                    section = 'alignment'
                    new_dict = getAlignmentHash_from_picardMarkDup(inputs[0], readset, record[section])

                elif step == 'metrics':
                    section = 'alignment'
                    new_dict = getAlignmentHash(inputs, readset, record[section])

                else:
                    new_dict = None

                if new_dict:
                    record[section] = new_dict
                break
        else:
            sys.exit("Error - no such sample " + readset + " in the provided JSON (" + json_file + ")")

    run_validation_str = json.dumps(run_report_json, indent=4)

    # Print to file
    with open(json_file, 'w') as out_json:
        out_json.write(run_validation_str)

def usage():
    print("\n-------------------------------------------------------------------------------------")
    print("runProcessingMetricsToJson.py" )
    print("This program was written by Edouard Henrion")
    print("For more information, contact: edouard.henrion@computationalgenomics.ca")
    print("----------------------------------------------------------------------------------\n")
    print("USAGE : runProcessingMetricsToJson.py")
    print("    -j    JSON file to be updated")
    print("    -s    step of the pipeline calling for update (basecall, index, fastq, fastqc, qc_graphs, blast, picard_mark_duplicates, metrics)")
    print("    -r    readset for which to fetch metrics and update the JSON")
    print("    -p    sequencing platform (illumina, mgig400, mgit7) : this guides the tool in determining which inputs it has to search for")
    print("    -i    input metrics file being parsed to update the JSON (could be multiple input metrics files)")
    print("    -h    this help\n")

def main():
    json_file, step, platform, inputs, readset = getarg(sys.argv)

    # finally (unlock) will execute even if exceptions occur
    try:

        # Make sure the json_file is unlock if process receives SIGTERM too (not python exception)
        def sigterm_handler(_signo, _stack_frame):
            unlock(json_file)
            sys.exit(0)
        signal.signal(signal.SIGTERM, sigterm_handler)

        # First lock the file to avoid multiple and synchronous writing attemps
        lock(json_file)

        if readset:
            print("Updating " + platform + "run processing report (" + json_file + ") for readset " + readset + " on step " + step)
        else:
            print("Updating " + platform + "run processing report (" + json_file + ") for step " + step)

        report(
            json_file,
            step,
            platform,
            inputs,
            readset
        )

    finally:
        # Finally unlock the file
        unlock(json_file)

if __name__ == '__main__':
    main()
