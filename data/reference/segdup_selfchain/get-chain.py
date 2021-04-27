#!/usr/bin/env python

import pandas as pd

df = pd.read_csv("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/chainSelf.txt.gz", compression="gzip",
                header=None, sep="\t")
df.columns = ["bin", "score", "tName", "tSize", "tStart", "tEnd", "qName", "qSize", "qStrand", "qStart", "qEnd", "id", "normScore"]

## Select for regions with an alignment score >= 90.0
df = df[df["normScore"] >= 90.0]

## Define chromosomes of interest
chrom = list(range(1,22+1))
chrom.extend(['X', 'Y'])
chrom = list(map(str, chrom))
chr_to_keep = list("chr" + c for c in chrom)

## Subset self-chain file for chr of interest            
selfchain = df.loc[df['tName'].isin(chr_to_keep) & df['qName'].isin(chr_to_keep)]

lines = []
for index, row in selfchain.iterrows():
    if not row['tName'].startswith('chr'): continue
    if not row['qName'].startswith('chr'): continue

    a = [row['tName'], row['tStart'], row['tEnd']]
    b = [row['qName'], row['qStart'], row['qEnd']]

    lines.append(a)
    lines.append(b)

lines.sort()
for l in lines:
    l[1] = str(l[1])
    l[2] = str(l[2])
    print("\t".join(l))

