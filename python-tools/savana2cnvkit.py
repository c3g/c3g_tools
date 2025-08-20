#! python3
##!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse

def convert_savana_to_cnvkit(input_file, output_file):
    # Read the input file
    df = pd.read_csv(input_file, sep='\t')

    # Clean or replace problematic values
    df['copyNumber'] = pd.to_numeric(df['copyNumber'], errors='coerce')
    df['minorAlleleCopyNumber'] = pd.to_numeric(df['minorAlleleCopyNumber'], errors='coerce')

    # Avoid division by zero or log2(0)
    df['copyNumber'].replace(0, np.nan, inplace=True)

    # Compute Segment_Mean safely
    df['Segment_Mean'] = np.log2(df['copyNumber'] / 2)

    # Compute nMinor and nMajor safely
    df['nMinor'] = (df['minorAlleleCopyNumber']).round()
    df['nMajor'] = (df['copyNumber'] - df['minorAlleleCopyNumber']).round()

    # Drop rows with non-finite values
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=['Segment_Mean', 'nMinor', 'nMajor'])

    # Convert to integer types safely
    df['nMinor'] = df['nMinor'].astype(int)
    df['nMajor'] = df['nMajor'].astype(int)

    # Select final columns
    df_out = df[['chromosome', 'start', 'end', 'Segment_Mean', 'nMajor', 'nMinor']]
    df_out.columns = ['Chromosome', 'Start', 'End', 'Segment_Mean', 'nMajor', 'nMinor']

    # Save to file
    df_out.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert Savana CNV file to CNVkit-like format.')
    parser.add_argument('--input', required=True, help='Input Savana CNV file (TSV)')
    parser.add_argument('--output', required=True, help='Output CNVkit-like file (TSV)')
    args = parser.parse_args()

    convert_savana_to_cnvkit(args.input, args.output)

