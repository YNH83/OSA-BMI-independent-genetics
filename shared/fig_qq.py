#!/usr/bin/env python3
"""Supplementary Figure 1 (QQ plots), one panel per GWAS (OSA, BMI, F_NonBMI), with lambda GC and LDSC intercept.
Points are thinned for plotting (all P < 1e-3 kept, a 5% random sample of the rest, fixed seed); lambda uses all SNPs.
Raw data: figures/raw_data/SFig1_raw.xlsx (500 expected/observed quantiles per GWAS + summary).
Usage (from the repository root): python3 shared/fig_qq.py"""
import os, re, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figkit as fk

STD = 'data/hm3_std'
T = '.'
G = [('OSA (FinnGen R12)', f'{STD}/OSA_FG.tsv.gz', 'P', 'OSA_FG', fk.INK2), ('BMI (GIANT)', f'{STD}/BMI_GIANT.tsv.gz', 'P', 'BMI_GIANT', fk.ORANGE_LIGHT),
     ('F$_{NonBMI}$ (after QC)', f'{T}/derived/gwas/F_NonBMI_QC.tsv.gz', 'P', 'F_NonBMI', fk.BLUE)]
ldlog = f'{T}/results/05_ldsc_matrix_ldsc.log'
log = open(ldlog).read() if os.path.exists(ldlog) else ''


def ldsc_int(trait):
    if not trait: return None
    fl = f'{T}/results/15_ldsc_{trait}_ldsc.log'   # residual components: single-trait LDSC on the post-QC GWAS (ldsc_factor_h2.R)
    if os.path.exists(fl):
        m = re.search(r'Intercept: ([\d.]+ \([\d.]+\))', open(fl).read()); return m.group(1) if m else None
    b = log.split(f'Heritability Results for trait: {trait}.sumstats.gz')
    if len(b) < 2: return None
    m = re.search(r'Intercept: ([\d.]+ \([\d.]+\))', b[1]); return m.group(1) if m else None


G = [g for g in G if os.path.exists(g[1])]
rng = np.random.default_rng(20260926)
fig, axes = plt.subplots(1, len(G), figsize=(2.5 * len(G), 2.6), squeeze=False)
qrows, srows, inputs = [], [], [ldlog]
for ax, (lab, f, pc, tr, col) in zip(axes[0], G):
    p = pd.read_csv(f, sep='\t', usecols=[pc])[pc].to_numpy(float); inputs.append(f)
    p = np.sort(p[(p > 0) & (p <= 1)]); n = len(p)
    lam = np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n); obs = -np.log10(p)
    keep = (p < 1e-3) | (rng.random(n) < 0.05)
    ax.scatter(exp[keep], obs[keep], s=2, color=col, linewidths=0, rasterized=True)
    mx = exp.max(); ax.plot([0, mx], [0, mx], color=fk.INK2, lw=0.7)
    li = ldsc_int(tr)
    ax.text(0.04, 0.97, f'λGC = {lam:.2f}' + (f'\nLDSC intercept\n{li}' if li else ''), transform=ax.transAxes, va='top', ha='left', fontsize=6.5, linespacing=1.3)
    if li: inputs.append(f'{T}/results/15_ldsc_{tr}_ldsc.log') if os.path.exists(f'{T}/results/15_ldsc_{tr}_ldsc.log') else None
    ax.set_title(lab, fontsize=7.5, loc='left', color=col, fontweight='bold'); ax.set_xlabel('Expected −log$_{10}$ P', fontsize=7); ax.set_ylabel('Observed −log$_{10}$ P', fontsize=7)
    qi = np.unique(np.linspace(0, n - 1, 500).astype(int))
    plain = lab.replace('$_{', '_').replace('}$', '')
    qrows.append(pd.DataFrame(dict(gwas=plain, rank=qi + 1, expected_mlog10P=exp[qi], observed_mlog10P=obs[qi])))
    srows.append(dict(gwas=plain, n_snps=n, lambda_GC=lam, LDSC_intercept=li, min_P=p[0], n_P_lt_5e8=int((p < 5e-8).sum())))
fig.tight_layout()
out = fk.save_fig(fig, T + '/figures', 'SupplementaryFigure1_QQ')
S = pd.DataFrame(srows)
fk.write_raw_xlsx(f'{T}/figures/raw_data/SFig1_raw.xlsx', {'summary': S, 'quantiles': pd.concat(qrows)}, inputs, os.path.abspath(__file__),
                  {'lambda_GC': 'median(chi2 from P) / 0.4549, all SNPs', 'LDSC_intercept': 'input GWAS: LDSC matrix log; residual component: single-trait LDSC on the post-QC GWAS with N = N_hat (SE)',
                   'expected_mlog10P': '-log10((rank - 0.5)/n)', 'observed_mlog10P': '-log10(P) at that rank'},
                  [(f'{r.gwas}: SNPs', r.n_snps) for r in S.itertuples()],
                  notes=['Plot thinning: all P < 1e-3 plus a 5% random sample of the rest (seed 20260926); lambda and quantiles use all SNPs.'])
print(S.to_string(index=False)); print('saved', out)
