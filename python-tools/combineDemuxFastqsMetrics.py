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
#
# Author : Edouard Henrion
# Contact : edouard.henrion@computationalgenomics.ca
#

import argparse
import os
import csv

def parse_metrics_file(
    metrics_file
    ):
    metrics_dict = {}
    #with open(filename, 'r') as metrics_file:
    metrics_tsv = csv.reader(metrics_file, delimiter='\t')
    headers = next(metrics_tsv)[1:]
    for row in metrics_tsv:
        metrics_dict[row[0]] = {key: str(value) for key, value in zip(headers, row[1:])}
    return metrics_dict

def sum_metrics(
    chunck_metrics,
    global_metrics
    ):
    for sample in chunck_metrics:
        if sample in metrics:
            # Check if library matches
            if not metrics[sample]['library_name'] == chunck_metrics[sample]['library_name']:
                raise Exception("Error : library should be the same for " + sample + ".\n  CONFLICT: " +  metrics[sample]['library_name'] + " vs. " + chunck_metrics[sample]['library_name'])
            # Check if library matches
            if not metrics[sample]['barcode'] == chunck_metrics[sample]['barcode']:
                raise Exception("Error : barcode should be the same for " + sample + ".\n  CONFLICT: " +  metrics[sample]['barcode'] + " vs. " + chunck_metrics[sample]['barcode'])
            global_metrics[sample]['templates'] += int(chunck_metrics[sample]['templates'])
            global_metrics[sample]['pf_templates'] += int(chunck_metrics[sample]['pf_templates'])
            global_metrics[sample]['perfect_matches'] += int(chunck_metrics[sample]['perfect_matches'])
            global_metrics[sample]['pf_perfect_matches'] += int(chunck_metrics[sample]['pf_perfect_matches'])
            global_metrics[sample]['one_mismatch_matches'] += int(chunck_metrics[sample]['one_mismatch_matches'])
            global_metrics[sample]['pf_one_mismatch_matches'] += int(chunck_metrics[sample]['pf_one_mismatch_matches'])
        else:
            global_metrics[sample] = chunck_metrics[sample]
            global_metrics[sample]['templates'] = int(global_metrics[sample]['templates'])
            global_metrics[sample]['pf_templates'] = int(global_metrics[sample]['pf_templates'])
            global_metrics[sample]['perfect_matches'] = int(global_metrics[sample]['perfect_matches'])
            global_metrics[sample]['pf_perfect_matches'] = int(global_metrics[sample]['pf_perfect_matches'])
            global_metrics[sample]['one_mismatch_matches'] = int(global_metrics[sample]['one_mismatch_matches'])
            global_metrics[sample]['pf_one_mismatch_matches'] = int(global_metrics[sample]['pf_one_mismatch_matches'])
    return global_metrics

def finalize_metrics(
    metrics_dict
    ):
    # Finish combining metrics (mean, ratio, etc...)
    sum_templates = sum([metrics_dict[sample]['templates'] for sample in metrics_dict])
    max_templates = max([metrics_dict[sample]['templates'] for sample in metrics_dict])
    sum_pf_templates = sum([metrics_dict[sample]['pf_templates'] for sample in metrics_dict])
    max_pf_templates = sum([metrics_dict[sample]['pf_templates'] for sample in metrics_dict])
    avg_pf_templates = float(sum_templates) / len(metrics_dict)
    for sample in metrics:
        metrics_dict[sample]['fraction_matches'] = metrics_dict[sample]['templates'] / float(sum_templates)
        metrics_dict[sample]['ratio_this_barcode_to_best_barcode'] = metrics_dict[sample]['templates'] / float(max_templates)
        metrics_dict[sample]['pf_fraction_matches'] = metrics_dict[sample]['pf_templates'] / float(sum_pf_templates)
        metrics_dict[sample]['pf_ratio_this_barcode_to_best_barcode'] = metrics_dict[sample]['pf_templates'] / float(max_pf_templates)
        metrics_dict[sample]['pf_normalized_matches'] = metrics_dict[sample]['pf_templates'] / avg_pf_templates
    return metrics

def print_metrics(
    metrics_dict,
    output_file
    ):
    with open(output_file, 'w') as metrics_file:
        header = []
        for i, sample in enumerate([sample for sample in metrics_dict]):
            # Print header if not already done
            if i == 0:
                header = list(metrics_dict[sample].keys())
                metrics_file.write("\t".join(['barcode_name'] + header) + "\n")
            # Print current sample metrics
            metrics_file.write("\t".join([sample] + [str(metrics_dict[sample][head]) for head in header]) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Combine fgbio DemuxFasts metrics files. First, for each demux barcode sum all the counts metrics from all the input files. Then rebuild all the ratio metrics.""")
    parser.add_argument("-i", "--input_files", help="List of input metrics files to combine ; format of the input is the fgbio DemuxFastqs metric format.", nargs="+", type=argparse.FileType('r'), required=True)
    parser.add_argument("-o", "--output", help="Output metric file, result of the combination of all the input metrics files.", type=str, required=True)
    args = parser.parse_args()

    # Loop over all the input files and build a pre-combined table
    metrics = {}
    for filename in args.input_files:
        chunck_metrics = parse_metrics_file(
            filename
        )
        # make the sum of all the <counts> columns 
        metrics = sum_metrics(
            chunck_metrics,
            metrics
        )

    # Finish combining metrics (mean, ratio, etc...)
    metrics = finalize_metrics(metrics)

    # Print combined metrics to output file
    print_metrics(
        metrics,
        args.output
    )

