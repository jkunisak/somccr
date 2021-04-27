#!/usr/bin/env python

import pandas as pd
import sys

## Read in the file
df = pd.read_csv("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz", compression="gzip", header=None, sep="\t")
df = df.iloc[:, [1,2,3, 5]]
df.columns = ['chr', 'start', 'stop', 'score'] ## Score based on the raw BLAST alignment score. Set to 0 and not used in later versions.

## Define chromosomes to keep
chrom = list(range(1,22+1))
chrom.extend(['X', 'Y'])
chrom = list(map(str, chrom))
chr_to_keep = list("chr" + c for c in chrom)

## Subset segmental duplications file for chr of interest
segdup = df.loc[df['chr'].isin(chr_to_keep)]

lines = []
for index, row in segdup.iterrows():
    a = [row['chr'], row['start'], row['stop']]
    lines.append(a)
    
lines.sort()
for l in lines: 
    l[1] = str(l[1])
    l[2] = str(l[2])
    print("\t".join(l))
