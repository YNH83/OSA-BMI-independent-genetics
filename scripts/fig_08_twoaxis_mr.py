#!/usr/bin/env python3
"""C9 Stage 2 figure: (a) two-axis test, correlation of F_NonBMI with T2D versus a second axis (A1, A1b, A2, post hoc same registry);
(b) bidirectional MR between BMI and OSA (A3). Raw data: figures/raw_data/Fig_twoaxis_mr_raw.xlsx.
Also writes Supplementary Table 10 (MR full results + instruments) and Supplementary Figure 5 (MR scatter and funnel)."""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, 'shared')
import figkit as fk

T = '.'
R = T + '/results'
TESTS = [('A1_HTN', 'Hypertension\n(FinnGen)', 'A1, primary'), ('A1b_CAD', 'Coronary artery\ndisease', 'A1b'),
         ('A2_SBP', 'Systolic BP\n(adj. BMI)', 'A2'), ('A2_DBP', 'Diastolic BP\n(adj. BMI)', 'A2'),
         ('posthoc_HTN_vs_T2DFG', 'Hypertension vs\nT2D, both FinnGen', 'post hoc')]
rows, inputs = [], []
for lab, name, tag in TESTS:
    p = f'{R}/07_twoaxis_{lab}.tsv'; inputs.append(p)
    r = pd.read_csv(p, sep='\t').iloc[0]
    rows.append(dict(test=lab, label=name.replace('\n', ' '), prereg=tag, metab=r.metab, axis2=r.axis2,
                     cor_metab=r.cor_NonBMI_metab, se_metab=r.se_metab, cor_axis2=r.cor_NonBMI_axis2, se_axis2=r.se_axis2,
                     chisq_1df=r.constrained_chisq, P=r.p_difference, grade=r.grade))
A = pd.DataFrame(rows)
mr = pd.read_csv(f'{R}/08_mr_A3.tsv', sep='\t'); inputs.append(f'{R}/08_mr_A3.tsv')
ins = pd.read_csv(f'{R}/08_mr_A3_instruments.tsv', sep='\t'); inputs.append(f'{R}/08_mr_A3_instruments.tsv')

fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.4), gridspec_kw=dict(width_ratios=[1.25, 1]))
y = np.arange(len(A))[::-1]
for yi, r in zip(y, A.itertuples()):
    a.errorbar(r.cor_metab, yi + 0.14, xerr=1.96 * r.se_metab, fmt='o', ms=4.5, color=fk.ORANGE, elinewidth=0.9, capsize=0)
    a.errorbar(r.cor_axis2, yi - 0.14, xerr=1.96 * r.se_axis2, fmt='s', ms=4.5, color=fk.BLUE, elinewidth=0.9, capsize=0)
    a.text(0.47, yi, f'P = {r.P:.2g}', va='center', fontsize=6.5, color=fk.INK, fontweight='bold' if r.P < 0.05 else 'normal')
a.axvline(0, color=fk.GREY, lw=0.5)
a.set_yticks(y); a.set_yticklabels([f'{n}\n[{t}]' for _, n, t in TESTS], fontsize=6.5)
a.set_xlim(-0.05, 0.56); a.set_xlabel('Genetic correlation with F$_{NonBMI}$ (95% CI)')
a.plot([], [], 'o', color=fk.ORANGE, label='type 2 diabetes (metabolic axis)'); a.plot([], [], 's', color=fk.BLUE, label='second axis')
a.legend(loc='upper center', bbox_to_anchor=(0.45, -0.17), ncol=2, fontsize=6.3)
a.set_title('a  Equality test: metabolic versus second axis', loc='left', fontsize=8, fontweight='bold')

meth = ['Inverse variance weighted', 'Weighted median', 'MR Egger', 'Weighted mode']
short = {'Inverse variance weighted': 'IVW', 'Weighted median': 'Weighted median', 'MR Egger': 'MR-Egger', 'Weighted mode': 'Weighted mode'}
B = []
for d, lab in [('BMI_to_OSA', 'BMI → OSA'), ('OSA_to_BMI', 'OSA → BMI')]:
    x = mr[mr.direction == d].set_index('method')
    for m_ in meth:
        B.append(dict(direction=d, label=lab, method=short[m_], nsnp=x.loc[m_, 'nsnp'], b=x.loc[m_, 'b'], lo=x.loc[m_, 'lo'], hi=x.loc[m_, 'hi'], P=x.loc[m_, 'pval']))
    f = x.iloc[0]
    B.append(dict(direction=d, label=lab, method='IVW, Steiger-concordant', nsnp=f.ivw_steiger_filtered_n, b=f.ivw_steiger_filtered_b,
                  lo=f.ivw_steiger_filtered_b - 1.96 * f.ivw_steiger_filtered_se, hi=f.ivw_steiger_filtered_b + 1.96 * f.ivw_steiger_filtered_se,
                  P=np.nan))
B = pd.DataFrame(B)
yy = []
pos = 0
for d in ['BMI_to_OSA', 'OSA_to_BMI']:
    for _ in range((B.direction == d).sum()):
        yy.append(-pos); pos += 1
    pos += 0.8
B['y'] = yy
for r in B.itertuples():
    c = fk.ORANGE if r.direction == 'BMI_to_OSA' else fk.BLUE
    b.errorbar(r.b, r.y, xerr=[[r.b - r.lo], [r.hi - r.b]], fmt='D' if 'IVW' in r.method else 'o', ms=4, color=c, elinewidth=0.9,
               mfc='white' if 'Steiger' in r.method else c)
    b.text(-0.22, r.y, f'{r.method} ({int(r.nsnp)})', va='center', ha='right', fontsize=6.3, color=fk.INK)
for d, lab in [('BMI_to_OSA', 'BMI → OSA (log odds per SD BMI)'), ('OSA_to_BMI', 'OSA → BMI (SD BMI per log odds OSA)')]:
    b.text(-0.95, B[B.direction == d].y.max() + 0.75, lab, fontsize=6.6, fontweight='bold', color=fk.INK,
           bbox=dict(facecolor='white', edgecolor='none', pad=0.5), zorder=5)
b.axvline(0, color=fk.GREY, lw=0.5)
b.set_yticks([]); b.spines['left'].set_visible(False)
b.set_xlim(-0.2, 1.05); b.set_ylim(B.y.min() - 0.8, B.y.max() + 1.4)
b.set_xlabel('MR estimate (95% CI)')
b.set_title('b  Bidirectional MR, BMI and OSA', loc='left', fontsize=8, fontweight='bold', x=-0.62)
fig.subplots_adjust(left=0.2, right=0.98, wspace=0.95, bottom=0.25, top=0.9)
# The main figure is superseded by Figure3v2_two_axes (fig_03_axes_v2.py, 2026-09-27); it is no longer saved.
# This script now only writes Supplementary Table 10 and Supplementary Figure 5.
plt.close(fig); out = None
if False:  # kept for reference only
  fk.write_raw_xlsx(f'{T}/figures/raw_data/Fig_twoaxis_mr_raw.xlsx', {'a_twoaxis': A, 'b_mr': B.drop(columns='y')}, inputs, os.path.abspath(__file__),
                  {'cor_metab': 'Genomic SEM correlation of F_NonBMI with T2D (single-indicator factor), free model',
                   'cor_axis2': 'correlation of F_NonBMI with the second-axis trait, free model',
                   'chisq_1df': 'chi-square of the model constraining the two unstandardized covariances equal (free model saturated, df 0)',
                   'b': 'MR estimate; BMI to OSA in log odds per SD BMI, OSA to BMI in SD BMI per unit log odds of OSA liability',
                   'nsnp': 'number of instruments after clumping, palindromic removal and harmonisation'},
                  [(f'{r.test}: 1-df chi-square P', r.P) for r in A.itertuples()] +
                  [('BMI to OSA IVW b', float(B[(B.direction == 'BMI_to_OSA') & (B.method == 'IVW')].b.iloc[0]))],
                  notes=['A2 blood-pressure GWAS (Evangelou 2018) were adjusted for BMI, which is collider-prone for this contrast.',
                         'The post hoc row compares two FinnGen registry diagnoses; it is not preregistered.'])

# Supplementary Table 10: MR full results and instruments
stab = mr.copy()
fk.write_raw_xlsx(f'{T}/figures/raw_data/STable10_raw.xlsx', {'mr_results': stab, 'instruments': ins}, inputs[-2:], os.path.abspath(__file__),
                  {'Q_ivw': 'Cochran Q for IVW', 'egger_intercept': 'MR-Egger intercept (directional pleiotropy)',
                   'n_steiger_correct': 'instruments whose exposure R2 exceeds outcome R2 (Steiger)', 'steiger_dir': 'TRUE if the instrument passes Steiger'},
                  [('instruments BMI_to_OSA', int((ins.direction == 'BMI_to_OSA').sum())), ('instruments OSA_to_BMI', int((ins.direction == 'OSA_to_BMI').sum()))])
stab.to_csv(f'{T}/figures/STable10_mr.tsv', sep='\t', index=False)

# Supplementary Figure 5: MR scatter and funnel plots for both directions
fig, ax = plt.subplots(2, 2, figsize=(7.0, 6.0))
for j, (d, xl, yl) in enumerate([('BMI_to_OSA', 'SNP effect on BMI (SD)', 'SNP effect on OSA (log odds)'),
                                  ('OSA_to_BMI', 'SNP effect on OSA (log odds)', 'SNP effect on BMI (SD)')]):
    h = ins[ins.direction == d].copy()
    s = np.sign(h['beta.exposure']); h['bx'] = h['beta.exposure'] * s; h['by'] = h['beta.outcome'] * s   # orient to exposure-increasing allele
    x = mr[mr.direction == d].set_index('method')
    s0 = ax[0, j]
    s0.errorbar(h.bx, h.by, xerr=h['se.exposure'], yerr=h['se.outcome'], fmt='none', ecolor=fk.GREY, elinewidth=0.4, alpha=0.6)
    s0.plot(h.bx, h.by, 'o', ms=2.8, color=fk.INK2 if d == 'OSA_to_BMI' else fk.ORANGE, mec='white', mew=0.3)
    xs = np.linspace(0, h.bx.max() * 1.05, 10)
    s0.plot(xs, x.loc['Inverse variance weighted', 'b'] * xs, color=fk.INK, lw=1, label='IVW')
    s0.plot(xs, x.loc['Weighted median', 'b'] * xs, color=fk.INK2, lw=1, ls=(0, (3, 2)), label='weighted median')
    s0.plot(xs, x['egger_intercept'].iloc[0] + x.loc['MR Egger', 'b'] * xs, color=fk.GREY, lw=1, ls=(0, (1, 1.5)), label='MR-Egger')
    s0.axhline(0, color=fk.GREY, lw=0.4); s0.set_xlabel(xl); s0.set_ylabel(yl); s0.legend(fontsize=6, loc='upper left')
    s0.set_title(f'{"ab"[j]}  {d.replace("_to_", " to ")}: scatter ({len(h)} SNPs)', loc='left', fontsize=8, fontweight='bold')
    s1 = ax[1, j]
    ratio = h.by / h.bx; prec = h.bx / h['se.outcome']
    s1.plot(ratio, prec, 'o', ms=2.8, color=fk.INK2 if d == 'OSA_to_BMI' else fk.ORANGE, mec='white', mew=0.3)
    s1.axvline(x.loc['Inverse variance weighted', 'b'], color=fk.INK, lw=1)
    s1.set_xlabel('Wald ratio estimate'); s1.set_ylabel('Precision (|beta exposure| / SE outcome)')
    s1.set_title(f'{"cd"[j]}  {d.replace("_to_", " to ")}: funnel', loc='left', fontsize=8, fontweight='bold')
fig.tight_layout()
fk.save_fig(fig, T + '/figures', 'SupplementaryFigure5_MR')
fk.write_raw_xlsx(f'{T}/figures/raw_data/SFig5_raw.xlsx', ins, inputs[-2:], os.path.abspath(__file__),
                  {'beta.exposure': 'instrument effect on the exposure (plotted after orienting to the exposure-increasing allele)',
                   'beta.outcome': 'instrument effect on the outcome, same allele'},
                  [(f'{d} IVW b (line slope)', float(mr[(mr.direction == d) & (mr.method == 'Inverse variance weighted')].b.iloc[0])) for d in ('BMI_to_OSA', 'OSA_to_BMI')])
print(A[['test', 'cor_metab', 'cor_axis2', 'P', 'grade']].round(4).to_string(index=False)); print(out)
