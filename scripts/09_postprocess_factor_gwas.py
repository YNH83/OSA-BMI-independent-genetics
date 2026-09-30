#!/usr/bin/env python3
"""C9 Stage 2 addendum 03, A9: post-processing of the multivariable GWAS (F_NonBMI and F_BMI).
QC per deviation record 02: drop SNPs with error != 0, with any warning other than the fixed lavaan 'slotOptions: se' notice,
or with Z_smooth > 1.96. Effective sample size N_hat = mean of 1 / (2p(1-p) SE^2) over SNPs with MAF 0.1 to 0.4.
Writes derived/gwas/<factor>_QC.tsv.gz in the shared standardized format (SNP A1 A2 BETA SE P N MAF), with N = N_hat,
and results/09_postprocess_summary.tsv. Usage: 09_postprocess_factor_gwas.py [tag=all]"""
import sys
import numpy as np
import pandas as pd

# Usage: 09_postprocess_factor_gwas.py [TAG=all] [PREFIX=userGWAS_] [TOPIC_DIR=.] [FACTORS=F_NonBMI,F_BMI]
TAG = sys.argv[1] if len(sys.argv) > 1 else 'all'
PREFIX = sys.argv[2] if len(sys.argv) > 2 else 'userGWAS_'
T = sys.argv[3] if len(sys.argv) > 3 else '.'
FACS = sys.argv[4].split(',') if len(sys.argv) > 4 else ['F_NonBMI', 'F_BMI']
OUTP = '' if PREFIX == 'userGWAS_' else PREFIX.replace('userGWAS_', '')
BENIGN = 'override(s) the options in slotOptions: se'
rows = []
for fac in FACS:
    d = pd.read_csv(f'{T}/derived/gwas/{PREFIX}{fac}_{TAG}.tsv.gz', sep='\t', low_memory=False)
    n0 = len(d)
    err = pd.to_numeric(d['error'], errors='coerce').fillna(1) != 0
    w = d['warning'].astype(str).str.replace('\n', ' ')
    benign_only = w.str.strip().isin(['0', '', 'nan']) | (w.str.contains(BENIGN, regex=False) &
                  (w.str.replace(r'lavaan->lav_lavaan_step02_options\(\):\s*the following argument\(s\) override\(s\) the options in slotOptions: se', '', regex=True)
                   .str.strip(' "') == ''))
    bad_warn = ~benign_only
    zs = pd.to_numeric(d['Z_smooth'], errors='coerce').fillna(0) > 1.96
    keep = ~(err | bad_warn | zs)
    q = d[keep].copy()
    m = q.MAF.between(0.1, 0.4)
    n_hat = float(np.mean(1 / (2 * q.MAF[m] * (1 - q.MAF[m]) * q.SE[m] ** 2)))
    out = pd.DataFrame(dict(SNP=q.SNP, A1=q.A1, A2=q.A2, BETA=q.est, SE=q.SE, P=q.Pval_Estimate, N=round(n_hat, 1), MAF=q.MAF))
    out.to_csv(f'{T}/derived/gwas/{OUTP}{fac}_QC.tsv.gz', sep='\t', index=False, compression='gzip', float_format='%.6g')
    chi = (q.Z_Estimate.astype(float)) ** 2
    rows += [(fac, 'SNPs in', n0), (fac, 'dropped: error', int(err.sum())), (fac, 'dropped: warning other than the fixed se notice', int((bad_warn & ~err).sum())),
             (fac, 'dropped: Z_smooth > 1.96', int((zs & ~err & ~bad_warn).sum())), (fac, 'SNPs after QC', len(q)),
             (fac, 'share dropped', round(1 - len(q) / n0, 5)), (fac, 'N_hat (MAF 0.1 to 0.4)', round(n_hat, 1)),
             (fac, 'SNPs used for N_hat', int(m.sum())), (fac, 'lambda GC', round(float(np.median(chi) / 0.4549364), 4)),
             (fac, 'GWS SNPs P < 5e-8', int((q.Pval_Estimate < 5e-8).sum())), (fac, 'GWS SNPs P < 2.5e-8', int((q.Pval_Estimate < 2.5e-8).sum()))]
S = pd.DataFrame(rows, columns=['factor', 'item', 'value'])
S.to_csv(f'{T}/results/09_postprocess_summary_{OUTP}{TAG}.tsv', sep='\t', index=False)
print(S.to_string(index=False))
