#!/usr/bin/env python3
"""Sensitivity table + forest plot.
Supplementary Table 6 and Supplementary Figure 3 (P1 and P2 across versions).
Usage (from the repository root): python3 shared/build_sensitivity.py"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figkit as fk

rows, inputs = [], []
T = '.'; R = T + '/results'
V = [('02_numeric_T2D_XUE.tsv', 'Primary: FinnGen OSA, GIANT BMI, T2D Xue'),
     ('02_numeric_T2D_XUE__OSA_FG_INCLAVO__BMI_GIANT.tsv', 'OSA definition including primary-care records'),
     ('02_numeric_T2D_FG.tsv', 'T2D from FinnGen'),
     ('02_numeric_T2D_XUE__OSA_CAMPOS__BMI_GIANT.tsv', 'OSA from Campos 2023 (overlaps FinnGen)')]
for f, lab in V:
    p = f'{R}/{f}'
    if not os.path.exists(p): continue
    d = pd.read_csv(p, sep='\t').iloc[0]; inputs.append(p)
    rows.append(dict(version=lab, metric='P1', est=d.P1, se=d.P1_se))
    rows.append(dict(version=lab, metric='P2', est=d.P2, se=d.P2_se))
thr = {'P1': (0.60, 0.40), 'P2': (0.50, 0.25)}
tlab = {'P1': ('pass ≥ 0.60', 'partial ≥ 0.40'), 'P2': ('pass ≥ 0.50', 'partial ≥ 0.25')}   # 2026-09-28: lowercase threshold labels
names = {'P1': 'P1  BMI-independent share of OSA genetic variance', 'P2': 'P2  attenuation of rg with type 2 diabetes'}
stab, sfig = 'STable6', 'SupplementaryFigure3_sensitivity'
D = pd.DataFrame(rows); D['lo'] = D.est - 1.96 * D.se; D['hi'] = D.est + 1.96 * D.se

mets = list(names)
fig, axes = plt.subplots(1, len(mets), figsize=(7.2, 0.45 * D.version.nunique() + 1.6), sharey=True, gridspec_kw=dict(wspace=0.15))
vers = list(dict.fromkeys(D.version))
y = np.arange(len(vers))[::-1]
for ax, m in zip(axes, mets):
    s = D[D.metric == m].set_index('version').reindex(vers)
    # analysis role is not a component, so no component colour (design spec rule 1): primary = black filled,
    # sensitivity = grey filled, null control = grey open; grey/white zebra rows (spec, supplementary figures)
    for i, yi in enumerate(y):
        if i % 2 == 0: ax.axhspan(yi - 0.5, yi + 0.5, color='#f3f2ef', zorder=0, lw=0)
    for i, (yi, (v, r)) in enumerate(zip(y, s.iterrows())):
        if pd.isna(r.est): continue
        c = fk.INK if i == 0 else fk.GREY
        ax.errorbar(r.est, yi, xerr=1.96 * r.se, fmt='o', color=c, ms=4.5, lw=1, capsize=0, zorder=3,
                    mfc='white' if 'Null' in v else c)
        _lines = [v for v in thr[m] if v is not None] + [0.0]   # value label kept clear of threshold and zero lines
        _span = (D.hi.max() - D.lo.min()) or 1.0; _near = [v for v in _lines if abs(r.est - v) < 0.04 * _span]
        _ha, _dx = ('center', 0.0) if not _near else (('left', 0.012 * _span) if r.est >= _near[0] else ('right', -0.012 * _span))
        ax.text(r.est + _dx, yi + 0.2, f'{r.est:.2f}'.replace('-', '\u2212'), ha=_ha, va='bottom', fontsize=7, color=fk.INK, zorder=4)
    a, b = thr[m]
    ax.axvline(a, color=fk.INK2, lw=0.7, ls=(0, (3, 2)))
    if b is not None: ax.axvline(b, color=fk.GREY, lw=0.5, ls=(0, (1, 2)))
    ax.axvline(0, color=fk.INK2, lw=0.5)
    # label the threshold lines so the dashed/dotted meaning is readable without the legend text
    top = len(vers) - 0.25
    # 3 pt gap between each caption and its line (layout pass 2026-09-27)
    FS_CAP = 6.5   # 2026-09-28: captions raised to >= 6.5 pt
    ax.annotate(tlab[m][0], (a, top), xytext=(3, 0), textcoords='offset points', fontsize=FS_CAP, color=fk.INK2, ha='left', va='center')
    # PARTIAL caption on two lines so it stays clear of the zero line
    if b is not None: ax.annotate(tlab[m][1].replace(' ≥', '\n≥'), (b, top), xytext=(-3, 0), textcoords='offset points', fontsize=FS_CAP, color=fk.GREY,
                                  ha='right', va='center', multialignment='right')
    ax.set_title(names[m], loc='left', fontsize=7.5, fontweight='bold')
    ax.set_ylim(-0.5, len(vers) + 0.35)
    ax.set_xlabel('Estimate (95% CI)')
axes[0].set_yticks(y); axes[0].set_yticklabels(vers, fontsize=7)
from matplotlib.lines import Line2D
_h = [Line2D([], [], marker='o', ls='', color=fk.INK, ms=4.5, label='primary analysis'), Line2D([], [], marker='o', ls='', color=fk.GREY, ms=4.5, label='sensitivity analysis')]
if D.version.str.contains('Null').any(): _h.append(Line2D([], [], marker='o', ls='', mfc='white', mec=fk.GREY, ms=4.5, label='null control'))
fig.legend(handles=_h, loc='upper center', ncol=len(_h), fontsize=7, frameon=False, bbox_to_anchor=(0.55, -0.01))   # below the axis titles; tight bbox keeps it in the saved file
out = fk.save_fig(fig, T + '/figures', sfig)
fk.write_raw_xlsx(f'{T}/figures/raw_data/{stab}_raw.xlsx', D, inputs, os.path.abspath(__file__),
                  {'est': 'point estimate', 'se': 'standard error (delta method for P1 and P2; Genomic SEM SE for rg)',
                   'lo': 'est - 1.96 se', 'hi': 'est + 1.96 se'},
                  [('dashed line = PASS threshold; dotted = PARTIAL threshold', str(thr))])
fk.write_raw_xlsx(f'{T}/figures/raw_data/SFig{sfig.split("_")[0][-1]}_raw.xlsx', D, inputs, os.path.abspath(__file__),
                  {'est': 'plotted point', 'lo': 'left whisker', 'hi': 'right whisker'})
D.to_csv(f'{T}/figures/{stab}_sensitivity.tsv', sep='\t', index=False)
print(D.round(3).to_string(index=False)); print('saved', out)
