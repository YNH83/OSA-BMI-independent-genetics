#!/usr/bin/env python3
"""Tissue enrichment figure (design spec: Foote 2025 Fig 4a small multiples), MAGMA gene-property on GTEx v8 (54 tissues).
One column per component, shared tissue order (sorted by the BMI P value), bars = -log10 P,
filled = Bonferroni-significant (0.05/54), pale otherwise. Also writes the gene-level and tissue-level supplementary tables.
Usage (from the repository root): python3 shared/fig_enrichment_v2.py
  OSA, BMI, F_NonBMI, gout residual (null)  -> Figure6_enrichment, STable8_magma, STable9_tissue, SupplementaryFigure4_null"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figkit as fk

GL = 'tools/genetics_env/magma/NCBI37.3/NCBI37.3.gene.loc'
T = '.'
COMP = [('OSA_FG', 'OSA', fk.INK2), ('BMI_GIANT', 'BMI (GIANT)', fk.ORANGE_LIGHT), ('F_NonBMI', 'F$_{NonBMI}$', fk.BLUE), ('GOUT_NonBMI', 'Gout minus BMI (null)', fk.GREY)]
KEY, FIG, NULL, SORTBY = 'F_NonBMI', 'Figure6_enrichment', 'GOUT_NonBMI', 'BMI_GIANT'
M = T + '/results/magma'
tis, genes, inputs = {}, {}, []
for k, lab, c in COMP:
    f = f'{M}/{k}.gtex.gsa.out'
    if not os.path.exists(f): sys.exit(f'missing MAGMA output {f}')
    inputs += [f, f'{M}/{k}.genes.out']
    d = pd.read_csv(f, sep=r'\s+', comment='#'); tis[k] = d.set_index('FULL_NAME')   # VARIABLE is truncated at 30 characters by MAGMA
    genes[k] = pd.read_csv(f'{M}/{k}.genes.out', sep=r'\s+')
NT = len(tis[KEY]); BONF = 0.05 / NT
order = tis[SORTBY].sort_values('P').index.tolist()   # ordered by the subtracted trait so its signal pattern is the visual reference
nice = lambda v: v.replace('___', ': ').replace('_', ' ')
# 2026-09-28 layout pass for journal artwork limits (final width <= 174 mm, all text >= 6.5 pt at final size)
S = dict(fs=(6.8, 7.7), wspace=0.14, title=7.5, pad=16, sub=6.5, xlab=7, tick=6.5, ytick=6.5, foot=6.5, verb='pass',
         adj=dict(left=0.285, right=0.965, top=0.935, bottom=0.085))
fig, axes = plt.subplots(1, len(COMP), figsize=S['fs'], sharey=True, gridspec_kw=dict(wspace=S['wspace']))
y = np.arange(len(order))[::-1]
for ax, (k, lab, c) in zip(axes, COMP):
    p = tis[k].loc[order, 'P'].values; v = -np.log10(p); sig = p < BONF
    ax.barh(y[sig], v[sig], color=c, height=0.72)                  # Bonferroni-significant: full colour
    ax.barh(y[~sig], v[~sig], color=c, alpha=0.28, height=0.72)    # not significant: pale
    ax.axvline(-np.log10(BONF), color=fk.INK, lw=0.6, ls=(0, (3, 2)))
    ax.set_title(lab, fontsize=S['title'], fontweight='bold', color=c, loc='left', pad=S['pad'])
    ax.text(0, 1.012, f'{int(sig.sum())} of {NT} {S["verb"]} Bonferroni', transform=ax.transAxes, ha='left', va='bottom', fontsize=S['sub'], color=fk.INK2)
    ax.set_xlabel('−log$_{10}$ P', fontsize=S['xlab']); ax.tick_params(labelsize=S['tick'])
    ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
axes[0].set_ylim(-0.7, len(order) - 0.3)   # no empty rows above the first or below the last tissue
axes[0].set_yticks(y); axes[0].set_yticklabels([nice(o) for o in order], fontsize=S['ytick'])
_foot = (f'MAGMA gene-property analysis, GTEx v8 (one-sided, conditional on average expression). Dashed line: Bonferroni 0.05/{NT}. '
         f'Tissues ordered by the {dict(BMI_GIANT="BMI")[SORTBY]} P value.')
_foot = _foot.replace('Dashed line', '\nDashed line')
fig.text(0.01, 0.005, _foot, fontsize=S['foot'], color=fk.INK2)
fig.subplots_adjust(**S['adj'])
fk.save_fig(fig, T + '/figures', FIG)
long = pd.concat([tis[k].reset_index().assign(component=k) for k, _, _ in COMP], ignore_index=True)
raw_name = 'Fig6_raw.xlsx'
fk.write_raw_xlsx(f'{T}/figures/raw_data/{raw_name}', long, inputs, os.path.abspath(__file__),
                  {'P': 'MAGMA gene-property P (one-sided)', 'BETA': 'regression coefficient of gene Z on tissue expression'},
                  [(f'Bonferroni threshold 0.05/{NT}', BONF)] + [(f'{k}: tissues below Bonferroni', int((tis[k].P < BONF).sum())) for k, _, _ in COMP])
# gene-level table
loc = pd.read_csv(GL, sep=r'\s+', header=None, names=['GENE', 'c', 's', 'e', 'st', 'SYMBOL'])
G = pd.concat([genes[k].merge(loc[['GENE', 'SYMBOL']], on='GENE', how='left').assign(component=k) for k, _, _ in COMP], ignore_index=True)
G['bonferroni'] = G.groupby('component').P.transform(lambda p: p < 0.05 / len(p))
Gs = G[G.bonferroni].sort_values(['component', 'P'])
Gs.to_csv(f'{T}/figures/STable8_magma.tsv', sep='\t', index=False)
fk.write_raw_xlsx(f'{T}/figures/raw_data/STable8_raw.xlsx', {'gene_wide_significant': Gs}, inputs, os.path.abspath(__file__),
                  {'P': 'MAGMA gene-based P', 'bonferroni': 'P < 0.05 / genes tested for that component'}, [(k, int((Gs.component == k).sum())) for k, _, _ in COMP])
long.to_csv(f'{T}/figures/STable9_tissue.tsv', sep='\t', index=False)
fk.write_raw_xlsx(f'{T}/figures/raw_data/STable9_raw.xlsx', long, inputs, os.path.abspath(__file__), {'P': 'MAGMA gene-property P'}, [('tissues', NT)])
# Supplementary Figure 4: null comparison, F_NonBMI vs gout residual (scatter of -log10 P per tissue)
fig, ax = plt.subplots(figsize=(4.3, 4.0))
a_, b_ = -np.log10(tis['F_NonBMI'].loc[order, 'P']), -np.log10(tis[NULL].loc[order, 'P'])
ax.scatter(a_, b_, s=12, color=fk.BLUE, edgecolor='white', linewidth=0.4)
lim = max(a_.max(), b_.max(), -np.log10(BONF)) * 1.1
ax.plot([0, lim], [0, lim], color=fk.INK2, lw=0.6, ls=(0, (3, 2)))
ax.axvline(-np.log10(BONF), color=fk.GREY, lw=0.5); ax.axhline(-np.log10(BONF), color=fk.GREY, lw=0.5)
ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_xlabel('F$_{NonBMI}$ −log$_{10}$ P', fontsize=7); ax.set_ylabel('Gout minus BMI −log$_{10}$ P', fontsize=7)
lab = list(dict.fromkeys(list(a_.sort_values().index[-3:]) + list(b_.sort_values().index[-3:])))
ax.text(-np.log10(BONF), lim * 0.35, 'Bonferroni', fontsize=6.5, color=fk.INK2, va='center', ha='right', rotation=90,
        transform=ax.transData + __import__('matplotlib').transforms.ScaledTranslation(-2.5 / 72, 0, fig.dpi_scale_trans))
ax.text(lim * 0.03, -np.log10(BONF), 'Bonferroni', fontsize=6.5, color=fk.INK2, va='bottom', ha='left',
        transform=ax.transData + __import__('matplotlib').transforms.ScaledTranslation(0, 2 / 72, fig.dpi_scale_trans))
fig.tight_layout()
# label placement (2026-09-27 layout pass): try offsets in turn and keep the first whose text box stays inside the axes and clears
# earlier labels, every point (with its marker radius), the y = x line, the Bonferroni lines and the two Bonferroni captions
from matplotlib.text import Text as _T
fig.canvas.draw(); rend = fig.canvas.get_renderer()
boxes = [_T.get_window_extent(t, rend).expanded(1.1, 1.3) for t in ax.texts]
pts = ax.transData.transform(np.c_[a_.values, b_.values]); rad = np.sqrt(12) / 2 * fig.dpi / 72 + 1.5
axbb = ax.get_window_extent(rend)
diag = ax.transData.transform(np.c_[np.linspace(0, lim, 400), np.linspace(0, lim, 400)])
bl = ax.transData.transform([(-np.log10(BONF), -np.log10(BONF))])[0]
def clear(bb, own):
    if bb.x0 < axbb.x0 + 2 or bb.x1 > axbb.x1 - 2 or bb.y0 < axbb.y0 + 2 or bb.y1 > axbb.y1 - 2: return 'axes'
    if any(bb.overlaps(o) for o in boxes): return 'label'
    if any(bb.x0 - rad < px < bb.x1 + rad and bb.y0 - rad < py < bb.y1 + rad for px, py in pts if abs(px - own[0]) + abs(py - own[1]) > 1): return 'point'
    if any(bb.contains(dx_, dy_) for dx_, dy_ in diag): return 'diag'
    if bb.x0 < bl[0] < bb.x1 or bb.y0 < bl[1] < bb.y1: return 'bonf'
    return ''
def unambiguous(bb, own):   # an unleadered label must sit nearer its own point than any other point
    dist = lambda px, py: np.hypot(max(bb.x0 - px, 0, px - bb.x1), max(bb.y0 - py, 0, py - bb.y1))
    return all(dist(px, py) > dist(*own) + 2 for px, py in pts if abs(px - own[0]) + abs(py - own[1]) > 1)
def leader_ok(bb, own):   # the leader (own point to the text box centre) must not cross an earlier label or pass over another point
    seg = np.linspace(own, ((bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2), 60)
    if any(o.contains(sx, sy) for o in boxes for sx, sy in seg): return False
    return not any(np.hypot(sx - px, sy - py) < rad for px, py in pts if abs(px - own[0]) + abs(py - own[1]) > 1 for sx, sy in seg[5:])
CAND = sorted([(dx, dy, 'left' if dx > 0 else ('right' if dx < 0 else 'center')) for dx in (0, 4, -4, 8, -8, 12, -12, 18, -18, 26, -26, 34, -34)
               for dy in (2, -2, 6, -6, 10, -10, 16, -16, 22, 28)], key=lambda c: np.hypot(c[0], c[1] * 1.3))
# two neighbouring brain tissues near (0.5 to 0.6, 1.2) get hand-set label positions (data coords) with leaders, set after visual check
MANUAL = {'Brain___Amygdala': (0.38, 1.62, 'center'), 'Brain___Cerebellum': (0.63, 1.45, 'left')}
for t in [t for t in lab if t in MANUAL]:
    mx, my, mha = MANUAL[t]
    tx = ax.annotate(nice(t), (a_[t], b_[t]), xytext=(mx, my), textcoords='data', fontsize=6.5, ha=mha, va='center',
                     arrowprops=dict(arrowstyle='-', lw=0.3, color=fk.GREY, shrinkA=0, shrinkB=2.5, relpos=(0, 0.5) if mha == 'left' else (0.5, 0)))
    tx.update_positions(rend); boxes.append(_T.get_window_extent(tx, rend).expanded(1.05, 1.15))
for t in sorted([t for t in lab if t not in MANUAL], key=lambda t: -a_[t]):   # right to left
    own = ax.transData.transform((a_[t], b_[t]))
    for dx, dy, ha in CAND:
        far = abs(dx) > 6 or abs(dy) >= 10
        tx = ax.annotate(nice(t), (a_[t], b_[t]), xytext=(dx, dy), textcoords='offset points', fontsize=6.5, ha=ha, va='bottom' if dy > 0 else 'top',
                         arrowprops=dict(arrowstyle='-', lw=0.3, color=fk.GREY, shrinkA=0, shrinkB=2.5) if far else None)
        tx.update_positions(rend); bb = _T.get_window_extent(tx, rend).expanded(1.05, 1.15)
        why = clear(bb, own) or ('' if not far or leader_ok(bb, own) else 'leader') or ('' if far or unambiguous(bb, own) else 'ambiguous')
        if not why: boxes.append(bb); break
        tx.remove()
    else:
        print('no clear label position for', t)
        ax.annotate(nice(t), (a_[t], b_[t]), xytext=(5, 3), textcoords='offset points', fontsize=6.5)
ax.set_title('Tissue enrichment: F$_{NonBMI}$ vs null residual', loc='left', fontsize=7.5, fontweight='bold')
fk.save_fig(fig, T + '/figures', 'SupplementaryFigure4_null')
fk.write_raw_xlsx(f'{T}/figures/raw_data/SFig4_raw.xlsx', pd.DataFrame({'tissue': order, 'mlog10P_F_NonBMI': a_.values, 'mlog10P_gout_residual': b_.values}),
                  inputs, os.path.abspath(__file__), {'mlog10P_gout_residual': 'same subtraction pipeline with gout in place of OSA'}, [('Bonferroni', BONF)])
print('C9', 'tissues', NT, {k: int((tis[k].P < BONF).sum()) for k, _, _ in COMP}, 'gene-wide sig', Gs.groupby('component').size().to_dict())
