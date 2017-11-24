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

# Python Standard Modules
import os
import csv
import numpy as np
import pandas as pd
import sys

# MUGQIC Modules

## extract genome assembly:
genome = sys.argv[1]

## open metrics/rnaseqRep/metrics.tsv file and load its content into an array

metrics_file1 = "metrics/rnaseqRep/metrics.tsv"


with open(metrics_file1) as f:
    reader = csv.reader(f, delimiter="\t")
    d = list(reader)

d = np.array(d)

## columns to keep:

keep1 = ["Sample", "Duplication Rate of Mapped", "Intergenic Rate", "Intragenic Rate", "Intronic Rate", "End 2 % Sense"]

header = ["Sample", "IntragenicRate", "IntergenicRate",  "IntronicRate", "DuplicationRate", "StrandSpecificity"]


columns_to_keep = []
for i, name in enumerate(d[:1][0]):
    if name in keep1:
        columns_to_keep.append(i)

metrics1 = pd.DataFrame(d[1:,columns_to_keep], columns=header)

#### extract data from report/trimAlignmentTable.tsv

header = ['Sample', 'RawReads', 'SurvivingReads', 'SurvivingReads_perc', 'AlignedReads',
        'AlignedReads_perc', 'AlternativeAlignments', 'AlternativeAlignments_perc', 'rRNAReads', 'rRNAReads_perc', 'Coverage', 'ExonicRate', 'GenesDetected']

metrics_file2 = "report/trimAlignmentTable.tsv"

with open(metrics_file2) as f:
    reader = csv.reader(f, delimiter="\t")
    d = list(reader)

d = np.array(d)
metrics2= pd.DataFrame(d[1:, :], columns=header)



## merge metrics1 and metrics 2 by Sample ID

All = pd.merge(metrics1, metrics2, on='Sample', how='outer')
All['genomeAssembly'] = genome

## calculate mitochondrial reads:
samples = list(All['Sample'])

MTreads={}

for sample in samples:
    cmd = "samtools view -c alignment/{sample}/{sample}.sorted.mdup.bam MT".format(sample =sample)
    MTreads[sample] = int(os.popen(cmd).read().strip())


MTreads = pd.DataFrame.from_dict(MTreads, 
                            orient='index')

MTreads.index.name = 'Sample'
MTreads.reset_index(inplace=True)

MTreads.columns = ['Sample','MitoReads']

All = pd.merge(All, MTreads, on='Sample', how='outer')


## reorder columns:

header = ['genomeAssembly', 'Sample', 'RawReads', 'SurvivingReads', 'SurvivingReads_perc', 'AlignedReads', 'AlignedReads_perc', 'AlternativeAlignments', 'AlternativeAlignments_perc', 'rRNAReads', 'rRNAReads_perc', 'DuplicationRate', 'Coverage', 'IntergenicRate', 'MitoReads', 'IntragenicRate', 'ExonicRate', 'IntronicRate', 'GenesDetected', 'StrandSpecificity']

All = All.loc[:,header]

file_name = 'report/IHEC_metrics_rnaseq_All.txt'
All.to_csv(file_name, sep='\t', index=False)


### write one file per sample:

for sample in samples:   
    file_name = 'report/IHEC_metrics_rnaseq_' + sample +'.txt'
    sub = All.loc[All['Sample'] == sample] 
    sub.to_csv(file_name, sep='\t', index=False)






