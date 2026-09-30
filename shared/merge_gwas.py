#!/usr/bin/env python3
"""Merge chunked userGWAS outputs into one file per factor (main OSA analysis, gout null control), checking that every chunk is present,
SNPs are unique, and (when the sumstats RDS row count is given) the total equals the prepared SNP count.
Chunk files are deleted after a verified merge (only the latest outputs are kept).
Usage: merge_gwas.py GWAS_DIR PREFIX NCHUNK FACTOR [FACTOR ...] [--expected N]
  e.g. merge_gwas.py derived/gwas userGWAS_ 8 F_BMI F_NonBMI
       merge_gwas.py derived/gwas null_gout_userGWAS_ 8 F_BMI F_NonBMI"""
import os, sys
import pandas as pd

a = sys.argv[1:]
expected = None
if '--expected' in a:
    i = a.index('--expected'); expected = int(a[i + 1]); a = a[:i] + a[i + 2:]
D, prefix, n, facs = a[0], a[1], int(a[2]), a[3:]
for fac in facs:
    files = [f'{D}/{prefix}{fac}_all_chunk{c:02d}of{n:02d}.tsv.gz' for c in range(1, n + 1)]
    missing = [f for f in files if not os.path.exists(f)]
    assert not missing, f'missing chunks: {missing}'
    d = pd.concat([pd.read_csv(f, sep='\t', low_memory=False) for f in files], ignore_index=True)
    assert d.SNP.is_unique, 'duplicate SNPs across chunks'
    if expected: assert len(d) == expected, f'{len(d)} rows, expected {expected}'
    nerr = int((d['error'].astype(str) != '0').sum())
    out = f'{D}/{prefix}{fac}_all.tsv.gz'
    d.to_csv(out, sep='\t', index=False, compression='gzip')
    chk = pd.read_csv(out, sep='\t', usecols=['SNP'])
    assert len(chk) == len(d)
    for f in files: os.remove(f)
    print(fac, 'rows', len(d), 'chunks', n, 'rows with lavaan error', nerr, 'GWS', int((d.Pval_Estimate < 5e-8).sum()),
          'min P', d.Pval_Estimate.min(), '| chunk files removed')
