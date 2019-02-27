#!/usr/bin/env python

import argparse

import os
import re
import sys
import pysam
from collections import defaultdict

from pysam import VariantFile

def name_to_gt(val):
    if val.lower() == "het":
        return (0,1)
    elif val.lower() == "hom":
        return (1,1)
    elif val.lower() in ["ref", "confict"]:
        return (0,0)
    else:
    # Non-standard representations, het is our best imperfect representation
    # print(fname, coords, ref, alt, info, val)
        return (0,1)

def alleles_to_gt(val, ref, alt):
    gt_indices = {gt.upper(): i for i, gt in enumerate(str(ref) + str(alt[0]))}
    tumor_gts = [gt_indices[x.upper()] for x in val if x in gt_indices]
    if tumor_gts and val not in known_names:
        if max(tumor_gts) == 0:
            tumor_gt = (0,0)
        elif 0 in tumor_gts:
            tumor_gt = (0,1)
            #tumor_gt = (0,"%s") % min([x for x in tumor_gts if x > 0])
        else:
            tumor_gt = (1,1)
            #tumor_gt = ("%s","%s") % (min(tumor_gts), max(tumor_gts))
    else:
        tumor_gt = name_to_gt(val)
    return tumor_gt


def process_vcfs(input_vcf, output_vcf, normal, tumor):
    input = VariantFile(input_vcf, 'rb')

    input.header.formats.add("GT", ".", "String", "Genotype")

    output = VariantFile(output_vcf, 'w', header=input.header)

    global known_names
    known_names = ["het", "hom", "ref", "conflict"]

    for variant in input:
        """Retrieve standard 0/0, 0/1, 1/1 style genotypes from INFO field.
        Normal -- NT field (ref, het, hom, conflict)
        Tumor -- SGT field
          - for SNPs specified as GG->TT for the normal and tumor diploid alleles. These
            can also represent more complex alleles in which case we set at heterozygotes
            pending longer term inclusion of genotypes in Strelka2 directly
            (https://github.com/Illumina/strelka/issues/16)
          - For indels, uses the ref, het, hom convention
        """
        nt_val = variant.info["NT"].split("=")[-1]
        normal_gt = name_to_gt(nt_val)

        sgt_val = variant.info["SGT"].split("->")[-1]

        if not sgt_val:
            tumor_gt = (0,0)
        else:
            tumor_gt = alleles_to_gt(sgt_val, variant.ref, variant.alts)

        variant.samples[normal]['GT'] = normal_gt
        variant.samples[tumor]['GT'] = tumor_gt

        #print str(normal_gt) + "\t" + str(tumor_gt) + "\t" + str(variant.samples['NORMAL']['GT']) + "\t" + str(variant.samples['TUMOR']['GT'])

        output.write(variant)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file',  default=sys.stdin, help="Input strelka vcf file")
    parser.add_argument('-o', '--output_file', help="Output vcf file")
    parser.add_argument('-n', '--normal', help="normal sample name")
    parser.add_argument('-t', '--tumor', help="normal sample name")
    args = parser.parse_args()
    args.logLevel = "INFO"

    sys.stderr.write("Process strelka VCF\n")
    process_vcfs(args.input_file, args.output_file, args.normal, args.tumor)
    if not os.path.exists(args.output_file + '.tbi'):
            sys.stderr.write('Indexing vcf file.\n')
            pysam.tabix_index(args.output_file, preset="vcf",force=True)
