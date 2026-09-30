#!/usr/bin/env python3
"""Supplementary Tables 1 to 3: LDSC heritability (S1), full rg matrix (S2), Genomic SEM estimates (S3).
Each written as a formatted tsv in figures/ and a raw-data Excel in figures/raw_data/.
LDSC intercepts and ratios are parsed from the GenomicSEM ldsc log of the full-matrix run.
Usage (from the repository root): python3 shared/build_stables_ldsc_sem.py"""
import glob, os, re, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figkit as fk

topic = 'C9'
T = '.'
R = T + '/results'
pre = '05'
h2 = pd.read_csv(f'{R}/{pre}_h2.tsv', sep='\t'); rg = pd.read_csv(f'{R}/{pre}_rg.tsv', sep='\t')
log = open(f'{R}/{pre}_ldsc_matrix_ldsc.log').read()
blocks = re.split(r'Heritability Results for trait: ', log)[1:]
qc = []
for b in blocks:
    name = b.split('.sumstats.gz')[0]
    g = lambda pat: (re.search(pat, b).group(1) if re.search(pat, b) else None)
    qc.append(dict(trait=name, mean_chi2=g(r'Mean Chi\^2 across remaining SNPs: ([\d.]+)'), lambda_GC=g(r'Lambda GC: ([\d.]+)'),
                   intercept=g(r'Intercept: ([\d.]+ \([\d.]+\))'), ratio=g(r'Ratio: ([-\d.]+ \([\d.]+\)|.*?)\n'),
                   observed_h2=g(r'Total Observed Scale h2: ([-\d.]+ \([\d.]+\))'), h2_z=g(r'h2 Z: ([-\d.]+)')))
S1 = h2.merge(pd.DataFrame(qc), on='trait', how='left')
S1 = S1.rename(columns={'liability_or_observed_h2': 'h2_used_in_SEM', 'intercept_x': 'intercept_matrix', 'intercept_y': 'intercept_log'})
fk.write_raw_xlsx(f'{T}/figures/raw_data/STable1_raw.xlsx', S1, [f'{R}/{pre}_h2.tsv', f'{R}/{pre}_ldsc_matrix_ldsc.log'], os.path.abspath(__file__),
                  {'h2_used_in_SEM': 'liability scale for binary traits (population prevalence = sample prevalence), observed for continuous',
                   'observed_h2': 'observed-scale h2 (SE) from the LDSC log', 'intercept_log': 'LDSC intercept (SE); near 1 indicates little confounding',
                   'ratio': '(intercept - 1) / (mean chi2 - 1); share of inflation not due to polygenicity'},
                  [(f'{r.trait}: h2 / SE', round(r.h2_used_in_SEM / r.h2_SE, 2)) for r in S1.itertuples()])
S1.to_csv(f'{T}/figures/STable1_ldsc_h2.tsv', sep='\t', index=False)

fk.write_raw_xlsx(f'{T}/figures/raw_data/STable2_raw.xlsx', rg, [f'{R}/{pre}_rg.tsv'], os.path.abspath(__file__),
                  {'rg': 'LDSC genetic correlation', 'SE': 'block-jackknife SE from GenomicSEM ldsc (standardized)', 'P': 'two-sided normal P = 2 x Phi(-|rg/SE|)'},
                  [(f'{r.trait1} vs {r.trait2}: z = rg / SE', round(r.rg / r.SE, 3)) for r in rg.head(5).itertuples()])
rg.to_csv(f'{T}/figures/STable2_rg.tsv', sep='\t', index=False)

est_files = sorted(glob.glob(f'{R}/*model_estimates*.tsv'))
S3 = pd.concat([pd.read_csv(f, sep='\t').assign(source_file=os.path.basename(f)) for f in est_files], ignore_index=True)
fk.write_raw_xlsx(f'{T}/figures/raw_data/STable3_raw.xlsx', S3, est_files, os.path.abspath(__file__),
                  {'Unstand_Est': 'unstandardized estimate (genetic covariance scale)', 'STD_Genotype': 'standardized on the genetic scale (loadings; correlations)',
                   'STD_Genotype_SE': 'SE of the standardized estimate', 'source_file': 'which model run (primary, sensitivity or null) produced the row'},
                  [('number of model runs', len(est_files))])
S3.to_csv(f'{T}/figures/STable3_sem_estimates.tsv', sep='\t', index=False)
print(topic, 'S1 rows', len(S1), 'S2 rows', len(rg), 'S3 rows', len(S3), 'from', len(est_files), 'runs')
