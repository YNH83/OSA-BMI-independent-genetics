#!/usr/bin/env python3
"""C9 Figure 5 and Supplementary Table 5: locus-level triangulation across three BMI-subtraction methods (design: Xu 2026 Fig 6A bubble matrix;
Zheng 2025 Fig 4 gene-name weight by method agreement).
Rows: the 32 FinnGen OSA loci (500 kb clumping, post-fix). Per locus the TRIANGULATION SNP is the SNP with the smallest OSA P within
500 kb of the OSA lead among SNPs present in all five GWAS (the OSA lead itself when available), so no locus is dropped for coverage. Columns: OSA (reference), F_NonBMI (Genomic SEM Cholesky), Campos 2023 mtCOJO
conditioned on BMI, Strausz 2025 OSA in BMI < 30, and BMI (GIANT, reference). Value: the OSA lead SNP in each GWAS, allele-aligned to the
OSA effect allele. Colour = same (wine) or opposite (green) direction as OSA; size = -log10 P (capped at 15); black outline = P < 5e-8;
pale = P >= 0.05; x = SNP absent from that GWAS. Gene name black bold = P < 0.05 and same direction in all three subtraction methods (grey otherwise).
Note: Campos 2023 includes FinnGen cases and Strausz 2025 is FinnGen, so the three methods are not independent samples.
Right-hand columns (2026-09-27): colocalisation (coloc, results/18_coloc.tsv) of F_NonBMI with BMI (GIANT) and with snoring (Campos) at the
5 loci retained in F_NonBMI; filled diamond = shared causal variant (PP.H4 >= 0.80), open diamond = distinct causal variants (PP.H3 >= 0.80),
grey bar = inconclusive; number = PP.H4. Rows are matched by gene name and checked against the coloc lead SNP position (within 500 kb)."""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, 'shared')
import figkit as fk

T = '.'
STD = 'data/hm3_std'
GL = 'tools/genetics_env/magma/NCBI37.3/NCBI37.3.gene.loc'
POS, NEG = '#9e2a55', '#3f7f55'
SRC = [('OSA', f'{STD}/OSA_FG.tsv.gz', 'OSA (FinnGen)'), ('F_NonBMI', f'{T}/derived/gwas/F_NonBMI_QC.tsv.gz', 'F$_{NonBMI}$ (Cholesky)'),
       ('mtCOJO', f'{STD}/OSA_CAMPOS_ADJBMI.tsv.gz', 'BMI-adjusted OSA (Campos)'), ('BMIlt30', f'{STD}/OSA_STRAUSZ_BMILT30.tsv.gz', 'OSA in BMI < 30 (Strausz)'),
       ('BMI', f'{STD}/BMI_GIANT.tsv.gz', 'BMI (GIANT)')]
METHODS = ['F_NonBMI', 'mtCOJO', 'BMIlt30']
L = pd.read_csv(f'{T}/results/04_OSA_FG_loci_all_fixP.tsv', sep='\t').sort_values(['CHR', 'BP']).reset_index(drop=True)
genes = pd.read_csv(GL, sep=r'\s+', header=None, names=['id', 'CHR', 'start', 'end', 'strand', 'gene'], dtype={'CHR': str})
genes = genes[~genes.gene.str.match(r'^(LOC|LINC|MIR|SNOR|FLJ|KIAA)|orf|-AS\d|^RNU')]
genes = genes[genes.CHR.str.isdigit()].assign(CHR=lambda g: g.CHR.astype(int))


def nearest(c, bp):
    g = genes[genes.CHR == c]
    dist = np.where((g.start <= bp) & (g.end >= bp), 0, np.minimum((g.start - bp).abs(), (g.end - bp).abs()))
    return g.gene.iloc[int(np.argmin(dist))]


L['gene'] = [nearest(r.CHR, r.BP) for r in L.itertuples()]
data = {k: pd.read_csv(f, sep='\t', usecols=['SNP', 'A1', 'A2', 'BETA', 'P']).set_index('SNP') for k, f, _ in SRC}
ref = data['OSA']
posd = pd.read_csv('tools/genetics_env/ref/reference.1000G.EUR.maf.txt', sep='\t', usecols=['SNP', 'CHR', 'BP']).set_index('SNP')
common = set.intersection(*[set(d.index) for d in data.values()])
tri = []
for r in L.itertuples():
    if r.lead_SNP in common: tri.append(r.lead_SNP); continue
    cand = posd[(posd.CHR == r.CHR) & ((posd.BP - r.BP).abs() <= 500_000)].index
    cand = [c for c in cand if c in common]
    tri.append(ref.loc[cand, 'P'].idxmin())
L['osa_lead_SNP'] = L.lead_SNP; L['lead_SNP'] = tri; L['triangulation_SNP_is_lead'] = L.lead_SNP == L.osa_lead_SNP
rows = []
for r in L.itertuples():
    a1 = ref.loc[r.lead_SNP, 'A1']
    for k, _, lab in SRC:
        d = data[k]
        if r.lead_SNP not in d.index:
            rows.append(dict(lead_SNP=r.lead_SNP, CHR=r.CHR, BP=r.BP, gene=r.gene, gwas=k, beta_aligned=np.nan, P=np.nan)); continue
        x = d.loc[r.lead_SNP]; s = 1 if x.A1 == a1 else (-1 if x.A2 == a1 else np.nan)
        rows.append(dict(lead_SNP=r.lead_SNP, CHR=r.CHR, BP=r.BP, gene=r.gene, gwas=k, beta_aligned=x.BETA * s, P=x.P))
D = pd.DataFrame(rows)
D['same_direction_as_OSA'] = np.sign(D.beta_aligned) == np.sign(D.merge(D[D.gwas == 'OSA'][['lead_SNP', 'beta_aligned']], on='lead_SNP', suffixes=('', '_osa')).beta_aligned_osa)
W = D.pivot(index='lead_SNP', columns='gwas', values=['P', 'same_direction_as_OSA'])
agree = pd.Series({s: all((W.loc[s, ('P', m)] < 0.05) and bool(W.loc[s, ('same_direction_as_OSA', m)]) for m in METHODS) for s in L.lead_SNP})
L['all3_nominal_same_direction'] = L.lead_SNP.map(agree)
L['n_methods_gws'] = [int(sum(W.loc[s, ('P', m)] < 5e-8 for m in METHODS)) for s in L.lead_SNP]

# ---------- coloc evidence (5 retained F_NonBMI loci) ----------
COLOC = f'{T}/results/18_coloc.tsv'
CO = pd.read_csv(COLOC, sep='\t')
CO = CO[(CO.trait1 == 'F_NonBMI') & CO.trait2.isin(['BMI_GIANT', 'SNORING_CAMPOS'])].copy()
CO = CO.merge(posd[['CHR', 'BP']].rename(columns={'CHR': 'coloc_lead_CHR', 'BP': 'coloc_lead_BP'}), left_on='lead_SNP', right_index=True, how='left')
g2row = L.set_index('gene')
CO['figure_row_osa_lead_SNP'] = CO.locus.map(g2row.osa_lead_SNP); CO['figure_row_triangulation_SNP'] = CO.locus.map(g2row.lead_SNP)
CO['figure_row_CHR'] = CO.locus.map(g2row.CHR); CO['figure_row_BP'] = CO.locus.map(g2row.BP)
CO['distance_coloc_lead_to_osa_lead_bp'] = (CO.coloc_lead_BP - CO.figure_row_BP).abs()
CO['row_match_ok'] = (CO.coloc_lead_CHR == CO.figure_row_CHR) & (CO.distance_coloc_lead_to_osa_lead_bp <= 500_000)
assert CO.locus.isin(L.gene).all() and CO.row_match_ok.all(), CO[['locus', 'lead_SNP', 'row_match_ok']]
assert (L.gene.value_counts()[CO.locus.unique()] == 1).all()
CO['PP_H4'] = CO['PP.H4']
CO_COLS = [('BMI_GIANT', 'coloc with BMI'), ('SNORING_CAMPOS', 'coloc with snoring')]

# 2026-09-28 layout pass for journal artwork limits (final height <= 234 mm, width <= 174 mm, all text >= 6.5 pt):
# wider and shorter canvas, column labels horizontal and wrapped above the matrix, legends and notes below in figure coordinates.
FW, FH = 6.4, 8.55
fig = plt.figure(figsize=(FW, FH))
AX_B, AX_T = 1.55, FH - 0.9   # axes bottom and top (inches)
ax = fig.add_axes([0.95 / FW, AX_B / FH, (FW - 0.95 - 0.08) / FW, (AX_T - AX_B) / FH])   # row pitch about 0.2 in so the largest bubbles do not touch
ymap = {s: len(L) - 1 - i for i, s in enumerate(L.lead_SNP)}
for j, (k, _, lab) in enumerate(SRC):
    for r in D[D.gwas == k].itertuples():
        y = ymap[r.lead_SNP]
        if np.isnan(r.P):
            ax.plot(j, y, 'x', color=fk.GREY, ms=4, mew=0.8); continue
        size = 6 + 9 * min(-np.log10(r.P), 15) ** 1.15
        c = POS if r.same_direction_as_OSA else NEG
        ax.scatter(j, y, s=size, color=c, alpha=1 if r.P < 0.05 else 0.25, edgecolor=fk.INK if r.P < 5e-8 else 'none', linewidth=0.9, zorder=3)
ax.axvspan(0.5, 3.5, color='#f3f2ef', zorder=0)
XC = [len(SRC) + 0.35 + 0.95 * i for i in range(len(CO_COLS))]   # coloc columns, right of the matrix
ax.axvline(len(SRC) - 0.4 + 0.2, color='#c9c8c3', lw=0.5, ymin=0.01, ymax=0.99)
for x, (t2, _) in zip(XC, CO_COLS):
    for r in CO[CO.trait2 == t2].itertuples():
        y = ymap[g2row.loc[r.locus, 'lead_SNP']]
        if r.call == 'shared causal variant':
            ax.plot(x - 0.2, y, 'D', ms=4.2, color=fk.INK, mec=fk.INK)
        elif r.call == 'distinct causal variants':
            ax.plot(x - 0.2, y, 'D', ms=4.2, mfc='white', mec=fk.INK, mew=0.8)
        else:
            ax.plot(x - 0.2, y, '_', ms=5, color=fk.GREY, mew=1.2)
        h4 = r.PP_H4
        ax.text(x + 0.0, y, '<0.01' if h4 < 0.005 else f'{h4:.2f}', ha='left', va='center', fontsize=7,
                color=fk.INK if r.call == 'shared causal variant' else fk.INK2, fontweight='bold' if r.call == 'shared causal variant' else 'normal')
COLLAB = {'OSA': 'OSA\n(FinnGen)', 'F_NonBMI': 'F$_{NonBMI}$\n(Cholesky)', 'mtCOJO': 'BMI-adjusted\nOSA\n(Campos)',
          'BMIlt30': 'OSA in\nBMI < 30\n(Strausz)', 'BMI': 'BMI\n(GIANT)'}
ax.xaxis.tick_top()
ax.set_xticks(list(range(len(SRC))) + XC)
ax.set_xticklabels([COLLAB[k] for k, _, _ in SRC] + ['coloc\nwith BMI', 'coloc\nwith snoring'], rotation=0, ha='center', va='bottom', fontsize=7,
                   linespacing=1.1)
ax.set_yticks([ymap[s] for s in L.lead_SNP])
ax.set_yticklabels([f'{g}  ' for g in L.gene], fontsize=7, style='italic')
for t, r in zip(ax.get_yticklabels(), L.itertuples()):
    # black bold = concordant in all three methods, grey = not (spec: weight grey to black by agreement); not blue, which marks F_NonBMI retention in Figure 4
    t.set_color(fk.INK if r.all3_nominal_same_direction else fk.GREY)
    if r.all3_nominal_same_direction: t.set_fontweight('bold')
ax.set_xlim(-0.6, XC[-1] + 0.62); ax.set_ylim(-0.7, len(L) - 0.3); ax.tick_params(length=0, pad=3)
for s_ in ax.spines.values(): s_.set_visible(False)
ax.annotate('three BMI-subtraction methods', (2, 1), xycoords=('data', 'axes fraction'), xytext=(0, 36), textcoords='offset points',
            ha='center', va='bottom', fontsize=7, color=fk.INK2, fontweight='bold')
for v, lab in ((5, 'P = 10$^{-5}$'), (15, 'P ≤ 10$^{-15}$')):
    ax.scatter([], [], s=6 + 9 * v ** 1.15, color=fk.GREY, label=lab)
ax.scatter([], [], s=60, color=POS, label='same direction as OSA'); ax.scatter([], [], s=60, color=NEG, label='opposite direction')
ax.scatter([], [], s=60, facecolor='white', edgecolor=fk.INK, label='P < 5 × 10$^{-8}$')
_h, _l = ax.get_legend_handles_labels()
leg1 = fig.legend(_h, _l, loc='upper center', bbox_to_anchor=(0.5, (AX_B - 0.08) / FH), ncol=3, fontsize=6.8, frameon=False, handletextpad=0.6,
                  columnspacing=1.6, labelspacing=1.3, borderpad=0.2)
from matplotlib.lines import Line2D
hco = [Line2D([], [], ls='none', marker='D', ms=4.2, color=fk.INK, label='shared causal variant (PP.H4 ≥ 0.80)'),
       Line2D([], [], ls='none', marker='D', ms=4.2, mfc='white', mec=fk.INK, mew=0.8, label='distinct causal variants (PP.H3 ≥ 0.80)'),
       Line2D([], [], ls='none', marker='_', ms=5, mew=1.2, color=fk.GREY, label='inconclusive')]
fig.legend(handles=hco, loc='upper center', bbox_to_anchor=(0.5, (AX_B - 0.62) / FH), ncol=3, fontsize=6.8, frameon=False, handletextpad=0.4,
           columnspacing=1.4, title='coloc of F$_{NonBMI}$ with BMI or snoring', title_fontsize=6.8, borderpad=0.2)
fig.text(0.04 / FW, 0.03 / FH, 'Black bold gene: P < 0.05 and same direction in all three methods. Pale: P ≥ 0.05. Per locus, the SNP with the smallest\n'
         'OSA P present in all five GWAS. The BMI-adjusted OSA GWAS (Campos 2023) includes FinnGen cases; the BMI < 30 GWAS is\n'
         'FinnGen: methods are not independent samples. coloc columns: F$_{NonBMI}$ vs BMI (GIANT) and vs snoring (Campos), 5 loci\n'
         'retained in F$_{NonBMI}$ only (blank = not tested); number = PP.H4.', fontsize=6.5, color=fk.INK2, va='bottom', linespacing=1.25)
fk.save_fig(fig, T + '/figures', 'Figure5_triangulation')
S5 = L[['osa_lead_SNP', 'lead_SNP', 'triangulation_SNP_is_lead', 'CHR', 'BP', 'gene', 'P', 'retained_in_NonBMI', 'all3_nominal_same_direction', 'n_methods_gws']] \
      .rename(columns={'P': 'P_OSA_lead', 'lead_SNP': 'triangulation_SNP'})
S5['lead_SNP'] = S5.triangulation_SNP
for k in ['F_NonBMI', 'mtCOJO', 'BMIlt30', 'BMI']:
    S5[f'beta_{k}'] = S5.lead_SNP.map(D[D.gwas == k].set_index('lead_SNP').beta_aligned)
    S5[f'P_{k}'] = S5.lead_SNP.map(D[D.gwas == k].set_index('lead_SNP').P)
S5.to_csv(f'{T}/figures/STable5_triangulation.tsv', sep='\t', index=False)
fk.write_raw_xlsx(f'{T}/figures/raw_data/Fig5_raw.xlsx', {'long': D, 'loci': S5, 'coloc': CO.drop(columns=['env'])}, [f for _, f, _ in SRC] + [f'{T}/results/04_OSA_FG_loci_all_fixP.tsv', GL, COLOC],
                  os.path.abspath(__file__), {'beta_aligned': 'effect of the OSA effect allele of the lead SNP (scale of each GWAS)',
                                              'all3_nominal_same_direction': 'P < 0.05 and same direction as OSA in F_NonBMI, mtCOJO and BMI < 30',
                                              'coloc sheet': 'rows of results/18_coloc.tsv with trait1 = F_NonBMI (plotted); call: PP.H4 >= 0.80 shared, PP.H3 >= 0.80 distinct, otherwise inconclusive',
                                              'row_match_ok': 'coloc lead SNP on the same chromosome and within 500 kb of the OSA lead of the figure row with that gene label'},
                  [('loci', len(L)), ('all three methods nominal and concordant', int(L.all3_nominal_same_direction.sum())),
                   ('coloc rows plotted', len(CO)), ('coloc rows matched to figure row', int(CO.row_match_ok.sum()))])
S5 = S5.drop(columns='lead_SNP')
S5.to_csv(f'{T}/figures/STable5_triangulation.tsv', sep='\t', index=False)
fk.write_raw_xlsx(f'{T}/figures/raw_data/STable5_raw.xlsx', S5, [f'{T}/figures/raw_data/Fig5_raw.xlsx'], os.path.abspath(__file__),
                  {'n_methods_gws': 'number of the three subtraction methods with P < 5e-8 at the OSA lead SNP'}, [('loci', len(S5))])
print(S5[['gene', 'retained_in_NonBMI', 'all3_nominal_same_direction', 'n_methods_gws', 'P_F_NonBMI', 'P_mtCOJO', 'P_BMIlt30', 'P_BMI']].to_string(index=False))
