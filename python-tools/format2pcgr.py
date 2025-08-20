#! python3
##!/usr/bin/env python

import argparse
import sys

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from cyvcf2 import VCF, Writer
import numpy as np

def process_variant(variant, callers, type, normal, tumor):
    if callers[0] == "clairS":
        variant.INFO['TAL'] = ",".join(callers)
        variant.INFO['TDP'] = variant.format('DP')[tumor].item()
        variant.INFO['TVAF'] = variant.format('AF')[tumor].item()
        variant.INFO['NDP'] = variant.format('NDP')[tumor].item()
        variant.INFO['NVAF'] = variant.format('NAF')[tumor].item()
        
    elif callers[0] == "varscan2":
        TDP4 = variant.format('DP4').item(tumor).split(",")
        NDP4 = variant.format('DP4').item(normal).split(",")
        variant.INFO['TAL'] = ",".join(callers)
        variant.INFO['TDP'] = int(TDP4[0]) + int(TDP4[1]) + int(TDP4[2]) + int(TDP4[3])
        variant.INFO['TVAF'] = float((int(TDP4[2]) + int(TDP4[3]))/ \
                                    (int(TDP4[0])+ int(TDP4[1]) + int(TDP4[2]) + int(TDP4[3])))
        variant.INFO['NDP'] = int(NDP4[0]) + int(NDP4[1]) + int(NDP4[2]) + int(NDP4[3])
        variant.INFO['NVAF'] = float((int(NDP4[2]) + int(NDP4[3])) / \
                                    (int(NDP4[0]) + int(NDP4[1]) + int(NDP4[2]) + int(NDP4[3])))
        #print(callers[0],DP4)

    elif callers[0] == "strelka2" and variant.is_indel and type == "somatic":
        TAR = variant.format('TAR')[tumor]
        TIR = variant.format('TIR')[tumor]
        NAR = variant.format('TAR')[normal]
        NIR = variant.format('TIR')[normal]
        variant.INFO['TAL'] = ",".join(callers)
        variant.INFO['TDP'] = int(TAR.item(0) + int(TIR.item(1)))

        if variant.INFO['TDP'] != 0:
            variant.INFO['TVAF'] = float(int(TAR.item(0))/ \
                                    (int(TAR.item(0)) + int(TIR.item(1))))
        else:
            variant.INFO['TVAF'] = 0

        variant.INFO['NDP'] = int(NAR.item(0)) + int(NIR.item(1))

        if variant.INFO['NDP'] != 0:
            variant.INFO['NVAF'] = float(int(NAR.item(0))/ \
                                    (int(NAR.item(0)) + int(NIR.item(1))))
        else:
            variant.INFO['NVAF'] = 0

    elif callers[0] == "GATK":
        TAD = variant.format('AD')[tumor]
        variant.INFO['TDP'] = int(TAD.item(0))+int(TAD.item(1))

        if variant.INFO['TDP'] != 0:
            variant.INFO['TVAF'] = float(int(TAD.item(1))/ \
                                    (int(TAD.item(0))+int(TAD.item(1))))
        else:
            variant.INFO['TVAF'] = 0

    #elif variant.INFO['PURPLE_AF'] :
    #    variant.INFO['TVAF'] = variant.INFO['PURPLE_AF']
    #    TAD = variant.format('AD')[0]
    #    NAD = variant.format('AD')[1]
    #    variant.INFO['TDP'] = int(TAD.item(0)) + int(TAD.item(1))
    #    variant.INFO['NDP'] = int(NAD.item(0)) + int(NAD.item(1))

    #    if int(NAD.item(0))+int(NAD.item(1)) != 0:
    #        variant.INFO['NVAF'] = float(int(NAD.item(1))/ \
    #                                (int(NAD.item(0))+int(NAD.item(1))))
    #    else:
    #        variant.INFO['NVAF'] = 0

    else:
        if variant.format('AD')[tumor].item(0) < 0:
            TDP =  variant.format('DP')[tumor]
            TAD = np.array([int(TDP),0])

        else:
            TAD = variant.format('AD')[tumor]

        if variant.format('AD')[normal].item(0) < 0:
            NDP =  variant.format('DP')[normal]
            NAD = np.array([int(NDP),0])

        else:
            NAD = variant.format('AD')[normal]

        variant.INFO['TAL'] = ",".join(callers)
        variant.INFO['TDP'] = int(TAD.item(0))+int(TAD.item(1))
        variant.INFO['NDP'] = int(NAD.item(0))+int(NAD.item(1))

        if variant.INFO['TDP'] != 0:
            variant.INFO['TVAF'] = float(int(TAD.item(1))/ \
                                    (int(TAD.item(0))+int(TAD.item(1))))
        else:
            variant.INFO['TVAF'] = 0

        if variant.INFO['NDP'] != 0:
                variant.INFO['NVAF'] = float(int(NAD.item(1))/ \
                                    (int(NAD.item(0))+int(NAD.item(1))))
        else:
                variant.INFO['NVAF'] = 0

    return variant

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', help="input vcf file", required=True)
    parser.add_argument('-o', '--output_file', help="output vcf", default="stdout")
    parser.add_argument('-v', '--variant_type', help="type of variant assessed (default:somatic)", default="somatic")
    parser.add_argument('-t', '--tumor_name', help="ordinal location of tumor in vcf", required=True)
    parser.add_argument('-f', '--filter', help="minimum number of caller support (default:1)", default=1)
    parser.add_argument('-fo', '--filter_only', help="filter based on number of callers only", action="store_true", default=False)
    parser.add_argument('-to', '--tumor_only', help="add TDP and TVAF to tumor only vcf", action="store_true", default=False)
    parser.add_argument('-tag', '--caller_tag', help="add caller name in cases on calling on one caller", required=False, default="GATK")
    args = parser.parse_args()
    args.logLevel = "INFO"

    sys.stderr.write("Process input VCF\n")

    input = VCF(args.input_file, threads=3)

    sample_list = input.samples

    if args.tumor_only:
        tumor_ordinal = sample_list.index(args.tumor_name)
        normal_ordinal = None
    else:
        if sample_list.index(args.tumor_name) == 0:
            normal_ordinal = 1
            tumor_ordinal = 0
        else:
            normal_ordinal = 0
            tumor_ordinal = 1

    sys.stderr.write(
        "Sample positions - Normal:" + str(normal_ordinal) +
        " Tumor: " + str(tumor_ordinal) + "\n"
    )

    if args.tumor_only == True:
        #input.add_info_to_header({'ID': 'TAL', 'Number': '1', 'Type': 'String', 'Description': 'Confidence of call i.e. number of callers supporting the variant'})
        input.add_info_to_header({'ID': 'TDP', 'Number': '1', 'Type': 'Integer', 'Description': 'Tumor depth derived from tumor AD field'})
        input.add_info_to_header({'ID': 'TVAF', 'Number': '1', 'Type': 'Float', 'Description': 'Tumor variant allele frequency derived from tumor AD field'})

    else:
        input.add_info_to_header({'ID': 'TAL', 'Number': '1', 'Type': 'String', 'Description': 'Confidence of call i.e. number of callers supporting the variant'})
        input.add_info_to_header({'ID': 'TDP', 'Number': '1', 'Type': 'Integer', 'Description': 'Tumor depth derived from tumor AD field'})
        input.add_info_to_header({'ID': 'TVAF', 'Number': '1', 'Type': 'Float', 'Description': 'Tumor variant allele frequency derived from tumor AD field'})
        input.add_info_to_header({'ID': 'NDP', 'Number': '1', 'Type': 'Integer', 'Description': 'Normal depth derived from tumor AD field'})
        input.add_info_to_header({'ID': 'NVAF', 'Number': '1', 'Type': 'Float', 'Description': 'Normal variant allele frequency derived from tumor AD field'})

    output = Writer(args.output_file, input)

    for variant in input:
        if not args.tumor_only and args.caller_tag:
            if args.caller_tag != "GATK":
                callers = [args.caller_tag]
            else:
                try:
                    callers = variant.INFO.get('CALLERS').split(",")
                except AttributeError:
                    callers = [args.caller_tag]

            if len(callers) >= int(args.filter) and args.filter_only:

                output.write_record(variant)

            elif len(callers) >= int(args.filter) and not args.filter_only:
                new_variant = process_variant(variant, callers, args.variant_type, normal_ordinal, tumor_ordinal)
                output.write_record(new_variant)

        elif args.tumor_only:
            callers = [args.caller_tag]
            new_variant = process_variant(variant, callers, args.variant_type, normal_ordinal, tumor_ordinal)
            output.write_record(new_variant)
    
    output.close()
    input.close()
