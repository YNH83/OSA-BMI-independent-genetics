#!/usr/bin/env python3
"""C9 Supplementary Figure 2: regional association plots of the 5 F_NonBMI lead loci (Table 4), one column per locus, rows =
OSA (FinnGen R12), BMI (GIANT) and F_NonBMI; points coloured by r2 with the F_NonBMI lead (1000G EUR, plink 1.9), gene track
below (MAGMA NCBI37.3 gene bodies, LOC/LINC/MIR/orf skipped). Window: lead +/- 500 kb (the locus definition).
Raw data: figures/raw_data/SFig2_raw.xlsx"""
import os, subprocess, sys, tempfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, 'shared')
import figkit as fk

T = '.'
STD = 'data/hm3_std'
REF = 'tools/genetics_env/ref/reference.1000G.EUR.maf.txt'
GL = 'tools/genetics_env/magma/NCBI37.3/NCBI37.3.gene.loc'
PLINK = 'tools/genetics_env/bin/plink19/plink'; BFILE = 'tools/genetics_env/magma/g1000_eur/g1000_eur'
W = 500_000
L = pd.read_csv(f'{T}/figures/Table4_loci.tsv', sep='\t').sort_values(['CHR', 'BP_hg19'])
pos = pd.read_csv(REF, sep='\t', usecols=['SNP', 'CHR', 'BP'])
src = {'OSA': (f'{STD}/OSA_FG.tsv.gz', 'P', fk.INK2), 'BMI': (f'{STD}/BMI_GIANT.tsv.gz', 'P', fk.ORANGE),
       'F_NonBMI': (f'{T}/derived/gwas/F_NonBMI_QC.tsv.gz', 'P', fk.BLUE)}
data = {k: pd.read_csv(f, sep='\t', usecols=['SNP', c]).rename(columns={c: 'P'}).merge(pos, on='SNP') for k, (f, c, _) in src.items()}
genes = pd.read_csv(GL, sep=r'\s+', header=None, names=['id', 'CHR', 'start', 'end', 'strand', 'gene'], dtype={'CHR': str})
genes = genes[~genes.gene.str.match(r'^(LOC|LINC|MIR|SNOR|FLJ|KIAA)|orf|-AS\d|^RNU')]
genes = genes[genes.CHR.str.isdigit()].assign(CHR=lambda g: g.CHR.astype(int))
cmap = LinearSegmentedColormap.from_list('r2', ['#d6d5d0', '#b9a6d6', '#7a4fb0', '#3d1a6e'])   # grey to purple: avoids the component colours


def r2_with(lead, chrom, lo, hi):
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, 'ld')
        subprocess.run([PLINK, '--bfile', BFILE, '--r2', '--ld-snp', lead, '--ld-window-kb', str(W // 1000 + 10), '--ld-window', '99999',
                        '--ld-window-r2', '0', '--chr', str(chrom), '--from-bp', str(lo), '--to-bp', str(hi), '--out', out],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        f = out + '.ld'
        if not os.path.exists(f): return pd.DataFrame(columns=['SNP', 'r2'])
        d = pd.read_csv(f, sep=r'\s+'); return d[['SNP_B', 'R2']].rename(columns={'SNP_B': 'SNP', 'R2': 'r2'})


n = len(L)
fig, axes = plt.subplots(4, n, figsize=(7.2, 6.2), gridspec_kw=dict(height_ratios=[1, 1, 1, 0.8], hspace=0.12, wspace=0.28),
                         sharex='col')
raw = []
for j, r in enumerate(L.itertuples()):
    lo, hi = r.BP_hg19 - W, r.BP_hg19 + W
    ld = r2_with(r.lead_SNP, r.CHR, lo, hi)
    for i, (k, (f, c, col)) in enumerate(src.items()):
        ax = axes[i, j]
        d = data[k][(data[k].CHR == r.CHR) & data[k].BP.between(lo, hi)].merge(ld, on='SNP', how='left')
        d['r2'] = d.r2.fillna(0); d = d.sort_values('r2')
        ax.scatter(d.BP / 1e6, -np.log10(d.P), c=d.r2, cmap=cmap, vmin=0, vmax=1, s=5, linewidths=0, rasterized=True)
        ld_lead = d[d.SNP == r.lead_SNP]
        if len(ld_lead): ax.plot(ld_lead.BP / 1e6, -np.log10(ld_lead.P), 'D', ms=5, color=col, mec='white', mew=0.6)
        ax.axhline(-np.log10(5e-8), color=fk.INK2, lw=0.5, ls=(0, (3, 2)))
        ax.tick_params(labelsize=5.8); ax.set_ylim(bottom=0)
        if j == 0: ax.set_ylabel(f"{ {'F_NonBMI': 'F$_{NonBMI}$'}.get(k, k)}\n−log$_{{10}}$ P", fontsize=6.5, color=col, fontweight='bold')
        if i == 0: ax.set_title(f'{r.nearest_gene}\n{r.lead_SNP}', fontsize=6.8, style='italic')
        raw.append(d.assign(locus=r.nearest_gene, gwas=k)[['locus', 'gwas', 'SNP', 'CHR', 'BP', 'P', 'r2']])
    g = genes[(genes.CHR == r.CHR) & (genes.end >= lo) & (genes.start <= hi)].copy()
    g['d'] = np.where((g.start <= r.BP_hg19) & (g.end >= r.BP_hg19), 0, np.minimum((g.start - r.BP_hg19).abs(), (g.end - r.BP_hg19).abs()))
    g = g.nsmallest(4, 'd').sort_values('start')   # the 5 genes closest to the lead, to keep labels legible
    gx = axes[3, j]; gx.set_yticks([]); gx.spines['left'].set_visible(False)
    cw = 0.060 * (2 * W / 1e6)                     # approximate label width per character (Mb) at 5.2 pt in this column width
    ends = []                                      # greedy row packing so that no two labels or bodies overlap
    for gr in g.itertuples():
        a, b = max(gr.start, lo) / 1e6, min(gr.end, hi) / 1e6
        half = len(gr.gene) * cw / 2
        xl = min(max((a + b) / 2, lo / 1e6 + half), hi / 1e6 - half)
        left, right = min(a, xl - half), max(b, xl + half)
        row = next((q for q, e in enumerate(ends) if e < left - 0.01), len(ends))
        if row == len(ends): ends.append(right)
        else: ends[row] = right
        gx.plot([a, b], [-row, -row], color=fk.INK2, lw=1.6, solid_capstyle='butt')
        gx.text(xl, -row + 0.18, gr.gene, fontsize=5.2, style='italic', color=fk.INK, ha='center', va='bottom')
    gx.set_ylim(-3.4, 1.0)
    gx.set_xlim(lo / 1e6, hi / 1e6); gx.set_xlabel(f'chr{r.CHR} (Mb, hg19)', fontsize=6); gx.tick_params(labelsize=5.8)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
fig.canvas.draw(); top, bot = axes[0, -1].get_position(), axes[2, -1].get_position()
cax = fig.add_axes([top.x1 + 0.012, bot.y0, 0.009, top.y1 - bot.y0])   # own axis, so the gene row keeps the same column width
cb = fig.colorbar(sm, cax=cax); cb.set_label('r$^2$ with the F$_{NonBMI}$ lead (1000G EUR)', fontsize=6)
cb.ax.tick_params(labelsize=5.5)
out = fk.save_fig(fig, T + '/figures', 'SupplementaryFigure2_regional')
R = pd.concat(raw, ignore_index=True)
fk.write_raw_xlsx(f'{T}/figures/raw_data/SFig2_raw.xlsx', R, [f for f, _, _ in src.values()] + [f'{T}/figures/Table4_loci.tsv', GL], os.path.abspath(__file__),
                  {'r2': 'LD r2 with the F_NonBMI lead SNP in 1000G EUR (plink 1.9); 0 when not computable', 'gwas': 'OSA, BMI or F_NonBMI'},
                  [(f'{r.nearest_gene} window', f'chr{r.CHR}:{r.BP_hg19 - W}-{r.BP_hg19 + W}') for r in L.itertuples()])
print(out)
