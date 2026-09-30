#!/usr/bin/env python3
"""C9 Figure 4 (v2, following the internal figure design specification): Miami plot of OSA (FinnGen R12, top) and the
BMI-independent component F_NonBMI (bottom), HapMap3 SNPs, MHC removed. Each half has its own y scale so the weaker
multivariable signal is not flattened (Foote 2025 Fig 3; Xu 2026 Fig 3A). Top in grey shades, bottom in F_NonBMI blue shades.
OSA loci (P3 rule, 500 kb): retained = filled blue diamond with italic nearest gene; lost = hollow grey circle (top 6 labelled).
Nearest gene = closest gene body in MAGMA NCBI37.3 gene locations (LOC, LINC, MIR and orf-type names skipped).
Raw data: figures/raw_data/Fig4v2_raw.xlsx (P < 1e-5 SNPs and loci) + Fig4v2_full_points.tsv.gz."""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, 'shared')
import figkit as fk

T = '.'
STD = 'data/hm3_std'
REF = 'tools/genetics_env/ref/reference.1000G.EUR.maf.txt'
GL = 'tools/genetics_env/magma/NCBI37.3/NCBI37.3.gene.loc'
NB = f'{T}/derived/gwas/F_NonBMI_QC.tsv.gz'
C_NON = fk.BLUE

pos = pd.read_csv(REF, sep='\t', usecols=['SNP', 'CHR', 'BP'])
osa = pd.read_csv(f'{STD}/OSA_FG.tsv.gz', sep='\t', usecols=['SNP', 'P']).rename(columns={'P': 'P_OSA'})
nb = pd.read_csv(NB, sep='\t', usecols=['SNP', 'P', 'N']).rename(columns={'P': 'P_NonBMI'})
n_hat = float(nb.N.iloc[0])
d = osa.merge(nb[['SNP', 'P_NonBMI']], on='SNP').merge(pos, on='SNP').sort_values(['CHR', 'BP'])
d = d[~((d.CHR == 6) & d.BP.between(25e6, 35e6))]
off = d.groupby('CHR').BP.max().cumsum().shift(fill_value=0)
d['x'] = d.BP + d.CHR.map(off)

genes = pd.read_csv(GL, sep=r'\s+', header=None, names=['id', 'CHR', 'start', 'end', 'strand', 'gene'], dtype={'CHR': str})
genes = genes[~genes.gene.str.match(r'^(LOC|LINC|MIR|SNOR|FLJ|KIAA)|orf|-AS\d|^RNU')]
genes = genes[genes.CHR.str.isdigit()].assign(CHR=lambda g: g.CHR.astype(int))


def nearest(c, bp):
    g = genes[genes.CHR == c]
    dist = np.where((g.start <= bp) & (g.end >= bp), 0, np.minimum((g.start - bp).abs(), (g.end - bp).abs()))
    i = int(np.argmin(dist)); return g.gene.iloc[i], int(dist[i])


L = pd.read_csv(f'{T}/results/04_OSA_FG_loci_all_fixP.tsv', sep='\t')
Ln = pd.read_csv(f'{T}/results/04_NonBMI_loci_all_fixP.tsv', sep='\t')
for tab in (L, Ln):
    ng = [nearest(r.CHR, r.BP) for r in tab.itertuples()]
    tab['nearest_gene'] = [g for g, _ in ng]; tab['distance_bp'] = [k for _, k in ng]
    tab['x'] = tab.BP + tab.CHR.map(off)

# 2026-09-28 layout pass for journal artwork limits: 6.55 in wide canvas (final width <= 174 mm), all text >= 6.5 pt
fig = plt.figure(figsize=(6.55, 4.6))
top = fig.add_axes([0.1, 0.53, 0.885, 0.40]); bot = fig.add_axes([0.1, 0.10, 0.885, 0.40], sharex=top)
g = -np.log10(5e-8)
ev = d.CHR % 2 == 0
# y-axis break for the OSA panel (design spec: avoid one peak flattening the rest): above CAP the scale is compressed 8-fold,
# marked by a break on the axis; tick labels show true -log10 P.
CAP, SQ = 16.0, 8.0
tr = lambda v: np.where(v > CAP, CAP + (v - CAP) / SQ, v)
top.scatter(d.x, tr(-np.log10(d.P_OSA)), s=1.1, c=np.where(ev, '#52514e', '#a9a8a4'), rasterized=True, linewidths=0)
bot.scatter(d.x, -np.log10(d.P_NonBMI), s=1.1, c=np.where(ev, C_NON, '#8fbbee'), rasterized=True, linewidths=0)
for ax in (top, bot):
    ax.axhline(g, color=fk.INK2, lw=0.6, ls=(0, (3, 2)))
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ymax_t = float(tr(-np.log10(d.P_OSA.min()))) * 1.22; ymax_b = max(-np.log10(d.P_NonBMI.min()) * 1.35, 10)
top.set_ylim(0, ymax_t); bot.set_ylim(ymax_b, 0)
hi_ticks = [t for t in (20, 30, 40, 50) if t <= -np.log10(d.P_OSA.min()) + 2]
top.set_yticks([0, 5, 10, 15] + [float(tr(t)) for t in hi_ticks]); top.set_yticklabels(['0', '5', '10', '15'] + [str(t) for t in hi_ticks], fontsize=7)
for yy in (CAP + 0.25, CAP + 0.75):   # break marks on the y axis
    top.plot([-0.012, 0.012], [yy - 0.35, yy + 0.35], transform=top.get_yaxis_transform(), color=fk.INK2, lw=0.8, clip_on=False)
top.axhspan(CAP + 0.2, CAP + 0.8, color='white', zorder=3, lw=0)
bot.set_yticks([t for t in (0, 5, 10, 15) if t <= ymax_b]); bot.tick_params(axis='y', labelsize=7)
top.set_ylabel('OSA\n−log$_{10}$ P', fontsize=7.5); bot.set_ylabel('F$_{NonBMI}$\n−log$_{10}$ P', fontsize=7.5)
plt.setp(top.get_xticklabels(), visible=False); top.tick_params(axis='x', length=0)
# OSA loci: retained vs lost
lost_lab = L[~L.retained_in_NonBMI].nsmallest(6, 'P').lead_SNP.tolist()
for r in L.itertuples():
    y = float(tr(-np.log10(r.P)))   # marker drawn at the lead SNP itself
    if r.retained_in_NonBMI:
        top.plot(r.x, y, 'D', ms=5, color=C_NON, mec='white', mew=0.5, zorder=4)
        # CDKN2B label right-aligned (2026-09-28) so it clears the lost-locus circle just right of it on chr9
        top.annotate(r.nearest_gene, (r.x, y), xytext=(3, 6) if r.nearest_gene == 'CDKN2B' else (0, 6), textcoords='offset points',
                     ha='right' if r.nearest_gene == 'CDKN2B' else 'center', fontsize=7, style='italic', color=C_NON, fontweight='bold')
    else:
        top.plot(r.x, y, 'o', ms=3.8, mfc='white', mec=fk.GREY, mew=0.8, zorder=4)
        if r.lead_SNP in lost_lab:
            # label offsets set after visual check (neighbouring peaks on chr11 and chr19/20)
            # 2026-09-28: BDNF label lifted above the chr9 to chr10 gap with a thin leader so it no longer touches the CDKN2B label
            dx, dy, ha = {'BDNF': (-22, 13, 'right'), 'GIPR': (-4, 4, 'right'), 'ZFP64': (4, 4, 'left')}.get(r.nearest_gene, (0, 4, 'center'))
            top.annotate(r.nearest_gene, (r.x, y), xytext=(dx, dy), textcoords='offset points', ha=ha, fontsize=6.8, style='italic', color=fk.INK2,
                         arrowprops=dict(arrowstyle='-', lw=0.4, color=fk.GREY, shrinkA=1, shrinkB=3, relpos=(1, 0)) if r.nearest_gene == 'BDNF' else None)
for r in Ln.itertuples():
    y = -np.log10(r.P)
    bot.annotate(r.nearest_gene, (r.x, y), xytext=(0, -8), textcoords='offset points', ha='center', va='top', fontsize=7, style='italic',
                 color=C_NON, fontweight='bold')
mids = d.groupby('CHR').x.median()
bot.set_xticks(mids.values); bot.set_xticklabels([str(c) if c <= 12 or c % 2 == 0 else '' for c in mids.index], fontsize=6.8)
bot.set_xlabel('Chromosome', fontsize=7.5); bot.set_xlim(d.x.min(), d.x.max())
nret = int(L.retained_in_NonBMI.sum())
top.text(0.005, 0.97, f'OSA, FinnGen R12 (N$_{{eff}}$ 201,550): {len(L)} loci', transform=top.transAxes, fontsize=7.2, va='top', color=fk.INK, fontweight='bold')
bot.text(0.005, 0.04, f'F$_{{NonBMI}}$, GWAS-by-subtraction (N$_{{hat}}$ {n_hat:,.0f}): {len(Ln)} loci', transform=bot.transAxes, fontsize=7.2,
         va='bottom', color=C_NON, fontweight='bold')
top.plot([], [], 'D', color=C_NON, ms=4.5, label=f'OSA locus retained in F$_{{NonBMI}}$ ({nret})')
top.plot([], [], 'o', mfc='white', mec=fk.GREY, ms=4, label=f'OSA locus lost ({len(L) - nret})')
top.legend(loc='upper left', bbox_to_anchor=(0.0, 0.9), fontsize=6.8, frameon=False, borderaxespad=0.2)
out = fk.save_fig(fig, T + '/figures', 'Figure4v2_miami')

full = f'{T}/figures/raw_data/Fig4v2_full_points.tsv.gz'
d[['SNP', 'CHR', 'BP', 'P_OSA', 'P_NonBMI']].to_csv(full, sep='\t', index=False, compression='gzip')
topsnp = d[(d.P_OSA < 1e-5) | (d.P_NonBMI < 1e-5)][['SNP', 'CHR', 'BP', 'P_OSA', 'P_NonBMI']]
fk.write_raw_xlsx(f'{T}/figures/raw_data/Fig4v2_raw.xlsx', {'snps_P_lt_1e-5': topsnp, 'OSA_loci_P3': L.drop(columns='x'), 'NonBMI_loci': Ln.drop(columns='x')},
                  [f'{STD}/OSA_FG.tsv.gz', NB, full, GL], os.path.abspath(__file__),
                  {'P_OSA': 'FinnGen R12 G6_SLEEPAPNO P', 'P_NonBMI': 'Genomic SEM F_NonBMI ~ SNP P (after A9 QC)',
                   'retained_in_NonBMI': 'any F_NonBMI SNP with P < 5e-8 within 500 kb of the OSA lead (preregistered P3 rule)',
                   'nearest_gene': 'closest gene body (MAGMA NCBI37.3), LOC/LINC/MIR/orf names skipped; 0 bp = lead inside the gene'},
                  [('points plotted', len(d)), ('OSA loci', len(L)), ('retained', nret), ('P3', round(nret / len(L), 4)), ('N_hat', n_hat)],
                  notes=['MHC (hg19 chr6 25 to 35 Mb) removed; positions hg19 from 1000G. Each half has its own y scale.'])
print(L[['lead_SNP', 'CHR', 'BP', 'P', 'retained_in_NonBMI', 'nearest_gene', 'distance_bp']].head(8).to_string(index=False))
print(Ln[['lead_SNP', 'CHR', 'BP', 'P', 'GWS_in_BMI_input', 'nearest_gene']].to_string(index=False)); print(out)
