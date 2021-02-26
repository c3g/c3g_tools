#!/usr/bin/env python

### Edouard Henrion (2020/03/18) - edouard.henrion@computationalgenomics.ca

import os
import sys
import getopt
import re
import json
import csv
import subprocess
from collections import namedtuple
import xml.etree.ElementTree as ET
from bs4 import BeautifulSoup

def getarg(argument):
        project=""
        sample=""
        lane=""
        sample_number=""
        library=""
        stats_json_file=""
        main_info_json=""
        qc_graph_xml=""
        blast_output=""
        kapa_tag_output=""
        alignment_directory=""
        gender=""
        output_file=""
        optli,arg = getopt.getopt(argument[1:],"p:s:n:l:i:m:q:b:k:a:g:o:h",['project','sample','lane_sample_number','library','index_stat','main_json','qc_graph','blast_output','kapa_tag','alignment_dir','gender','output','help'])
        if len(optli) == 0 :
            usage()
            sys.exit("Error : No argument given")
        for option, value in optli:
            if option in ("-p","--project"):
                project=str(value)
            if option in ("-s","--sample"):
                sample=str(value)
            if option in ("-n","--lane_sample_number"):
                try:
                    [lane, sample_number] = value.split("_")
                    lane = int(lane)
                    sample_number = int(sample_number)
                except:
                    sys.exit("Error - Lane Sample Number has to be 2 integers separated by '_':\n" + str(value) + " not allowed !!!")
            if option in ("-l","--library"):
                library=str(value)
            if option in ("-i","--index_stat"):
                stats_json_file=str(value)
            if option in ("-m","--main_json"):
                main_info_json=str(value)
            if option in ("-q","--qc_graph"):
                qc_graph_xml=str(value)
            if option in ("-b","--blast_output"):
                blast_output=str(value)
            if option in ("-k","--kapa_tag"):
                kapa_tag_output=str(value)
            if option in ("-a","--alignment_dir"):
                alignment_directory=str(value)
            if option in ("-g","--gender"):
                gender=str(value)
            if option in ("-o","--output"):
                output_file=str(value)
            if option in ("-h","--help"):
                usage()
                sys.exit()

        if stats_json_file and not os.path.exists(stats_json_file):
            sys.exit("Error - JSON index stats file not found:\n" + str(stats_json_file))
        if main_info_json and not os.path.exists(main_info_json):
            sys.exit("Error - general information JSON file not found:\n" + str(main_info_json))
        if not os.path.exists(qc_graph_xml):
            sys.exit("Error - qc XML file not found:\n" + str(qc_graph_xml))
        if not os.path.exists(blast_output):
            sys.exit("Error - Blast output file not found:\n" + str(blast_output))
        if kapa_tag_output and not os.path.exists(kapa_tag_output):
            sys.exit("Error - Kapa tag output file not found:\n" + str(kapa_tag_output))
        if alignment_directory and not os.path.exists(alignment_directory):
            sys.exit("Error - alignment directory not found:\n" + str(alignment_directory))

        return project, sample, sample_number, lane, library, stats_json_file, main_info_json, qc_graph_xml, blast_output, kapa_tag_output, alignment_directory, gender, output_file

def parseMetricsFile(
    metrics_file
    ):
    reader = open(metrics_file, 'rb').readlines()
    metrics_tsv = []
    header = []
    for row in reader:
        if row!="\n" and not row[0].startswith('#'):
            if not header:
                header = row.strip().split("\t")
            else:
                metrics_tsv.append(dict(zip(header, row.strip().split("\t"))))
    return metrics_tsv

def getIndexHash_from_BCL2fastq(
    stats_json_file,
    sample_name
    ):

    with open(stats_json_file, 'r') as json_file:
        stats_json = json.load(json_file)

    total_raw = stats_json['ConversionResults'][0]['TotalClustersRaw']
    total_pf = stats_json['ConversionResults'][0]['TotalClustersPF']
    total_pf_onindex_inlane = sum([sample['NumberReads'] for sample in stats_json['ConversionResults'][0]['DemuxResults']])

    for i, sample in enumerate(stats_json['ConversionResults'][0]['DemuxResults']):
        if sample['SampleName'] == sample_name: 
            return {   
                'Barcode sequence': sample['IndexMetrics'][0]['IndexSequence'],
                'Barcode': sample['IndexMetrics'][0]['IndexName'] if sample['IndexMetrics'][0]['IndexName'] else "",
                'PF Clusters': sample['NumberReads'],
                '% of the lane': 100*sample['NumberReads']/float(total_pf),
                '% on index in lane': 100*sample['NumberReads']/float(total_pf_onindex_inlane),
                '% Perfect barcode': 100*sample['IndexMetrics'][0]['MismatchCounts']['0']/float(sample['NumberReads']),
                '% One mismatch barcode': 100*sample['IndexMetrics'][0]['MismatchCounts']['1']/float(sample['NumberReads']),
                'Yield (bases)': sample['Yield'],
                '% >= Q30 bases': 100*sum([readMetrics['YieldQ30'] for readMetrics in sample['ReadMetrics']])/float(sample['Yield']),
                'Mean Quality Score': sum([readMetrics['QualityScoreSum'] for readsMetrics in sample['ReadMetrics']])/float(sample['Yield'])
            }

def getIndexHash_from_DemuxFastqs(
    demux_metrics_file,
    project,
    sample,
    sample_number,
    lane,
    library
    ):

    demux_metrics_tsv = parseMetricsFile(demux_metrics_file)

    total_pf = sum([int(row['pf_templates']) for row in demux_metrics_tsv])
    total_pf_onindex_inlane = sum([int(row['pf_templates']) for row in demux_metrics_tsv if row['library_name'] != "unmatched"])

    raw_barcode_keys = [row['barcode_name'] for row in demux_metrics_tsv if row['library_name'] == library]
    barcode_keys = []
    barcode_name = ""
    for name in raw_barcode_keys:
        m = re.search(sample+"_"+library+"_\w?_?(?P<barcode_name>.+)", name)
        if m:
            if barcode_name and barcode_name != m.group('barcode_name'):
                print "Error !! At least more than one barcode (" + barcode_name + ", " + m.group('barcode_name') + ") was found to be associated to library(" + sample + "_" + library + ")"
                sys.exit(32)
            barcode_name = m.group('barcode_name')
            barcode_keys.append(name)
        else:
            print name + " could not be parsed for barcode..."

    pf_clusters = sum([int(row['pf_templates']) for row in demux_metrics_tsv if row['library_name'] == library and row['barcode_name'] in barcode_keys])
    pf_perfect_matches = sum([int(row['pf_perfect_matches']) for row in demux_metrics_tsv if row['library_name'] == library and row['barcode_name'] in barcode_keys])
    pf_one_mismatch_matches = sum([int(row['pf_one_mismatch_matches']) for row in demux_metrics_tsv if row['library_name'] == library and row['barcode_name'] in barcode_keys])

    fastqc_file = os.path.join(os.path.dirname(demux_metrics_file), "Project_" + project, "Sample_" + sample + "_" + library, "fastqc.R1", sample + "_" + library + "_S" + str(sample_number) + "_L00" + str(lane) + "_R1_001_fastqc", "fastqc_data.txt")
    print fastqc_file
    q30_reads = int(subprocess.check_output("sed -n '/#Quality/,/END_MODULE/p' %s | awk '{if($1!=\">>END_MODULE\" && $1!=\"#Quality\" && $1>=30){sum+=$2}} END {print sum}'" % fastqc_file, shell=True).strip())
    seq_length = int(subprocess.check_output("grep 'Sequence length' %s | cut -f 2" % fastqc_file, shell=True).strip())
    total_readset_seq = int(subprocess.check_output("grep 'Total Sequences' %s | cut -f 2" % fastqc_file, shell=True).strip())
    readset_yield = seq_length * float(total_readset_seq)
    q30_bases = seq_length * (q30_reads)
    pct_q30_bases = 100 * q30_bases / float(readset_yield)
    mean_quality_score = float(subprocess.check_output("sed -n '/#Base\tMean/,/END_MODULE/p' %s | awk '{if($1!=\">>END_MODULE\" && $1!=\"#Mean\"){records+=1;sum+=$2}} END {print sum/records}'" % fastqc_file, shell=True).strip())

    return {
        'Barcode sequence': ','.join([row['barcode'] for row in demux_metrics_tsv if row['library_name'] == library and row['barcode_name'] in barcode_keys]),
        'Barcode': barcode_name if barcode_name else "",
        'PF Clusters': pf_clusters,
        '% of the lane': 100*pf_clusters/float(total_pf),
        '% on index in lane': 100*pf_clusters/float(total_pf_onindex_inlane),
        '% Perfect barcode': 100*pf_perfect_matches/float(pf_clusters),
        '% One mismatch barcode': 100*pf_one_mismatch_matches/float(pf_clusters),
        'Yield (bases)': seq_length*float(pf_clusters),
        '% >= Q30 bases': pct_q30_bases,
        'Mean Quality Score': mean_quality_score
    }

def getIndexHash_from_AlignmentMetrics(
    alignment_directory,
    main_info_json,
    sample,
    library
    ):

    with open(main_info_json, 'r') as json_file:
        main_json = json.load(json_file)

    alignment_summary_metrics_file = os.path.join(alignment_directory, sample + "." + library + ".sorted.metrics.alignment_summary_metrics")
    align_tsv = parseMetricsFile(alignment_summary_metrics_file)

    for i, readset_index in enumerate(main_json["barcodes"][sample + "_" + library]):
        index_name = readset_index['INDEX_NAME']
        if readset_index['INDEX1'] and readset_index['INDEX2']:
            index_seq = readset_index['INDEX1'] + "-" + readset_index['INDEX2'] 
        elif readset_index['INDEX1']:
            index_seq = readset_index['INDEX1']
        elif readset_index['INDEX2']:
            index_seq = readset_index['INDEX2']
        else:
            index_seq = ""
        if index_name:
            break

    return {
        'Barcode sequence': index_seq,
        'Barcode': index_name,
        'PF Clusters': align_tsv[0]['PF_READS'],
        '% of the lane': 0,
        '% on index in lane': 0,
        '% Perfect barcode': 0,
        '% One mismatch barcode': 0,
        'Yield (bases)': 0,
        '% >= Q30 bases': 0,
        'Mean Quality Score': 0
    }


def getQcHash(
    qc_graph_xml
    ):

    root = ET.parse(qc_graph_xml).getroot()
#    nb_reads = root.attrib['nbReads']
#    nb_bases = root.attrib['nbBases']
    return {
        'avgQual' : root.attrib['avgQual'],
        'duplicateRate' : root.attrib['duplicateRate'],
    }

def getBlastHash(blast_output):
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

    return {
        '1st_hit' : best_match.species + " (" + best_match.hits + ")" if best_match else None,
        '2nd_hit' : nd_hit.species + " (" + nd_hit.hits + ")" if nd_hit else None,
        '3rd_hit' : rd_hit.species + " (" + rd_hit.hits + ")" if rd_hit else None
    }
def getSampleTagHash():
    # Parsing of the kapa tag result file still remains to be done
    return None

def getAlignmentHash(
    alignment_directory,
    sample,
    library,
    gender
    ):

    alignment_summary_metrics_file = os.path.join(alignment_directory, sample + "." + library + ".sorted.metrics.alignment_summary_metrics")
    align_tsv = parseMetricsFile(alignment_summary_metrics_file)

    insert_size_metrics_file = os.path.join(alignment_directory, sample + "." + library + ".sorted.metrics.insert_size_metrics")
    insert_tsv = parseMetricsFile(insert_size_metrics_file)

    dup_metrics_file = os.path.join(alignment_directory, sample + "." + library + ".sorted.dup.metrics")
    dup_tsv = parseMetricsFile(dup_metrics_file)

    verify_bam_id_file = os.path.join(alignment_directory, sample + "." + library + ".sorted.metrics.verifyBamId.tsv")
    verifyBamID_tsv = parseMetricsFile(verify_bam_id_file)

    target_coverage_file = os.path.join(alignment_directory, sample + "." + library + ".sorted.metrics.targetCoverage.txt")
    if os.path.isfile(target_coverage_file):
        sex_match_reader = csv.DictReader(open(target_coverage_file, 'rb'), delimiter='\t')

        for row in sex_match_reader:
            #print row
            if row['IntervalName'] == "chrX":
                chrX_cov = row['MeanCoverage']
            if row['IntervalName'] == "chrY":
                chrY_cov = row['MeanCoverage']
            if row['IntervalName'] == "Total":
                total_cov = row['MeanCoverage']

        #print chrX_cov
        #print chrY_cov
        if float(chrX_cov) > 0.8:
            sex_det = "F"
        elif float(chrY_cov) > 0.25:
            sex_det = "M"
        else:
            sex_det = "?"
        if sex_det == gender:
            sex_match = True
        else:
            sex_match = False
    else:
        total_cov = "N/A"
        sex_det = "N/A"
        sex_match = "N/A"
    
    return {
        'pf_read_alignment_rate': align_tsv[2]['PCT_PF_READS_ALIGNED'],
        'mean_coverage': total_cov,
        'chimeras': align_tsv[2]['PCT_CHIMERAS'],
        'adapter_dimers': align_tsv[2]['PCT_ADAPTER'],
        'average_aligned_insert_size': insert_tsv[0]['MEAN_INSERT_SIZE'],
        'aligned_dup_rate': dup_tsv[0]['PERCENT_DUPLICATION'],
        'Freemix': verifyBamID_tsv[0]['FREEMIX'],
        'reported_sex': gender,
        'inferred_sex': sex_det,
        'sex_concordance': sex_match
    }

def report(
    project,
    sample,
    sample_number,
    lane,
    library,
    demultiplex_stats_file,
    main_info_json,
    qc_graph_xml,
    blast_output,
    kapa_tag_output,
    alignment_directory,
    gender,
    output_file
    ):

    # Now start building run validation hash table (list of dict for all the samples)
    sample_report_hash = {
        'project':project,
        'sample': sample+"_"+library,
        'qc': getQcHash(qc_graph_xml),
        'blast': getBlastHash(blast_output)
    }

    if demultiplex_stats_file:
        if demultiplex_stats_file.endswith('.json'):
            sample_report_hash['index'] = getIndexHash_from_BCL2fastq(demultiplex_stats_file, sample+"_"+library)
        elif demultiplex_stats_file.endswith('.DemuxFastqs.metrics.txt'):
            sample_report_hash['index'] = getIndexHash_from_DemuxFastqs(demultiplex_stats_file, project, sample, sample_number, lane, library)
        else:
            sys.exit("Error - unrecognized barcodes stats file :\n" + str(demultiplex_stats_file))
    elif alignment_directory:
        sample_report_hash['index'] = getIndexHash_from_AlignmentMetrics(alignment_directory, main_info_json, sample, library)
    else:
        sample_report_hash['index'] = None

    sample_report_hash['sample_tag'] = getSampleTagHash() if kapa_tag_output else None
    sample_report_hash['alignment'] = getAlignmentHash(alignment_directory, sample, library, gender) if alignment_directory else None

    run_validation_str = json.dumps(sample_report_hash, indent=4)

    # Print to file
    with open(output_file, 'w') as out_json:
        out_json.write(run_validation_str)

def usage():
        print "\n-------------------------------------------------------------------------------------"
        print "runProcessingSampleReport.py uses indexes from Clarity LIMS, passed in hash string,"
        print "to validate the indexes produced by CountIlluminaBarcodes and"
        print "output a report in JSON format."
        print "This program was written by Edouard Henrion"
        print "For more information, contact: edouard.henrion@computationalgenomics.ca"
        print "----------------------------------------------------------------------------------\n"
        print "USAGE : runProcessingSampleReport.py"
        print "    -p    project ID of the sample for which the report will be created"
        print "    -s    sample name for which the report will be created"
        print "    -n    sample number on lane, asually assigned by the pipeline"
        print "    -l    library name for which the report will be created"
        print "    -i    index stats file : could be a BCL2fastq JSON stats file or a DemuxFastqs metrics output file (will be determined upon file extension)"
        print "    -m    main_info_json"
        print "    -q    qc graph output xml file"
        print "    -b    blast output file"
        print "    -k    kapa_tag_output"
        print "    -a    folder containing all the alignemnt file metrics"
        print "    -g    gender"
        print "    -o    output_file"
        print "    -h    this help\n"

def main():
        project, sample, sample_number, lane, library, index_stats_file, main_info_json, qc_graph_xml, blast_output, kapa_tag_output, alignment_directory, gender, output_file = getarg(sys.argv)
        print main_info_json
        print "Building the run processing report for sample " + sample + " (library " + library + ") from project " + project

        report(project, sample, sample_number, lane, library, index_stats_file, main_info_json, qc_graph_xml, blast_output, kapa_tag_output, alignment_directory, gender, output_file)

main()

