#!/usr/bin/env python3
"""C9 Stage 2 addendum 03, A6: summarise LAVA output with the preregistered reading rules.
Significance: Benjamini-Hochberg FDR q < 0.05 within each trait pair (over the loci tested for that pair); Bonferroni also reported.
 A6a (descriptive): number of loci with significant local rg between F_NonBMI and each trait (HTN vs T2D contrast named in the rule).
 A6b (primary): partial local rg(OSA, X | BMI), X = HTN (FinnGen, Pan-UKB) and T2D: supported for HTN if >= 1 locus with q < 0.05 and positive.
 A6c: cancellation flag for F_NonBMI with FEV1/FVC, sleep duration and chronotype: >= 1 significant positive AND >= 1 significant negative locus.
 Added after the main-axis change (descriptive, not a preregistered rule): local rg of F_NonBMI with snoring.
Writes results/11_lava_A6_summary.tsv, figures/STable11_lava.tsv, figures/raw_data/STable11_raw.xlsx and Supplementary Figure 6."""
import glob, os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, 'shared')
import figkit as fk

T = '.'
D = T + '/derived/lava'
LOC = 'tools/genetics_env/src_pkgs/LAVA-main/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile'


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


rd = lambda pat: pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(f'{D}/{pat}_chunk*of04.tsv'))], ignore_index=True)
U, B, P = rd('univ'), rd('bivar'), rd('pcor')
# LAVA writes the target phenotype as phen2 in run.bivar output; put the target (NONBMI or OSA_FG) first
sw = B.phen2.isin(['NONBMI', 'OSA_FG']) & ~B.phen1.isin(['NONBMI', 'OSA_FG'])
B.loc[sw, ['phen1', 'phen2']] = B.loc[sw, ['phen2', 'phen1']].values
assert B.phen1.isin(['NONBMI', 'OSA_FG']).all(), B.phen1.unique()
loci = pd.read_csv(LOC, sep=r'\s+')
for df in (B, P):
    df['q_fdr'] = np.nan
    for (a, b), g in df.groupby(['phen1', 'phen2']):
        df.loc[g.index, 'q_fdr'] = bh(g.p.values)
        df.loc[g.index, 'n_tested'] = len(g)
    df['p_bonf_sig'] = df.p < 0.05 / df.n_tested
    df['sig'] = df.q_fdr < 0.05
B = B.merge(loci, left_on='locus', right_on='LOC', how='left').drop(columns='LOC')
P = P.merge(loci, left_on='locus', right_on='LOC', how='left').drop(columns='LOC')
nb = B[B.phen1 == 'NONBMI']
rows = []
for x, g in nb.groupby('phen2'):
    s = g[g.sig]
    rows.append(dict(analysis='A6a/A6c bivariate F_NonBMI', trait=x, loci_tested=len(g), sig_fdr=len(s), sig_pos=int((s.rho > 0).sum()),
                     sig_neg=int((s.rho < 0).sum()), sig_bonf=int(g.p_bonf_sig.sum())))
osa = B[B.phen1 == 'OSA_FG']
for x, g in osa.groupby('phen2'):
    s = g[g.sig]
    rows.append(dict(analysis='reference: bivariate OSA', trait=x, loci_tested=len(g), sig_fdr=len(s), sig_pos=int((s.rho > 0).sum()),
                     sig_neg=int((s.rho < 0).sum()), sig_bonf=int(g.p_bonf_sig.sum())))
pc_col = 'pcor' if 'pcor' in P.columns else [c for c in P.columns if 'pcor' in c or c == 'r'][0]
for x, g in P.groupby('phen2'):
    s = g[g.sig]
    rows.append(dict(analysis='A6b partial rg(OSA, X | BMI)', trait=x, loci_tested=len(g), sig_fdr=len(s), sig_pos=int((s[pc_col] > 0).sum()),
                     sig_neg=int((s[pc_col] < 0).sum()), sig_bonf=int(g.p_bonf_sig.sum())))
S = pd.DataFrame(rows)
# preregistered verdicts
def cnt(an, tr, col='sig_pos'):
    r = S[(S.analysis == an) & (S.trait == tr)]; return int(r[col].iloc[0]) if len(r) else 0
verd = []
for h in ('HTN_FG', 'HTN_PANUKB'):
    verd.append(('A6b', f'partial rg(OSA, {h} | BMI) >= 1 positive FDR locus', 'SUPPORTED' if cnt('A6b partial rg(OSA, X | BMI)', h) >= 1 else 'not supported'))
verd.append(('A6b', 'partial rg(OSA, T2D | BMI) positive FDR loci (reported, no rule)', str(cnt('A6b partial rg(OSA, X | BMI)', 'T2D_XUE'))))
verd.append(('A6a', 'F_NonBMI significant loci: HTN (FinnGen) / HTN (Pan-UKB) / T2D',
             f"{cnt('A6a/A6c bivariate F_NonBMI', 'HTN_FG', 'sig_fdr')} / {cnt('A6a/A6c bivariate F_NonBMI', 'HTN_PANUKB', 'sig_fdr')} / {cnt('A6a/A6c bivariate F_NonBMI', 'T2D_XUE', 'sig_fdr')}"))
for x in ('FEV1FVC_SHRINE', 'SLEEPDUR_DASHTI', 'CHRONOTYPE_JONES'):
    flag = cnt('A6a/A6c bivariate F_NonBMI', x, 'sig_pos') >= 1 and cnt('A6a/A6c bivariate F_NonBMI', x, 'sig_neg') >= 1
    verd.append(('A6c', f'cancellation flag F_NonBMI x {x}', 'FLAG' if flag else 'no flag'))
verd.append(('added (descriptive)', 'F_NonBMI x snoring significant loci (pos/neg)',
             f"{cnt('A6a/A6c bivariate F_NonBMI', 'SNORING_CAMPOS')}/{cnt('A6a/A6c bivariate F_NonBMI', 'SNORING_CAMPOS', 'sig_neg')}"))
V = pd.DataFrame(verd, columns=['item', 'rule', 'result'])
S.to_csv(f'{T}/results/11_lava_A6_summary.tsv', sep='\t', index=False)
V.to_csv(f'{T}/results/11_lava_A6_verdicts.tsv', sep='\t', index=False)
sigB = B[B.sig & B.phen1.isin(['NONBMI'])].sort_values(['phen2', 'CHR', 'START'])
sigP = P[P.sig].sort_values(['phen2', 'CHR', 'START'])
st = pd.concat([sigB.assign(analysis='bivariate F_NonBMI'), sigP.assign(analysis='partial rg(OSA, X | BMI)')], ignore_index=True)
st.to_csv(f'{T}/figures/STable11_lava.tsv', sep='\t', index=False)
fk.write_raw_xlsx(f'{T}/figures/raw_data/STable11_raw.xlsx', {'verdicts': V, 'counts': S, 'significant_loci': st, 'all_bivariate': B, 'all_pcor': P,
                  'univariate': U}, sorted(glob.glob(f'{D}/*_chunk*of04.tsv')) + [f'{D}/sample_overlap.txt', LOC], os.path.abspath(__file__),
                  {'rho': 'local genetic correlation (LAVA bivariate)', 'q_fdr': 'BH FDR within the trait pair', 'sig': 'q_fdr < 0.05',
                   'p_bonf_sig': 'p < 0.05 / loci tested for the pair'},
                  [(f'{r.item}: {r.rule}', r.result) for r in V.itertuples()])
# Supplementary Figure 6: (a) number of FDR-significant loci per trait (positive / negative), (b) forest of the loci for snoring and
# the two cancellation-flag traits. All 181 loci are in Supplementary Table 11.
NAMES = {'SNORING_CAMPOS': 'Snoring', 'INSOMNIA_FG': 'Insomnia (FG)', 'SLEEPDUR_DASHTI': 'Sleep duration', 'CHRONOTYPE_JONES': 'Chronotype',
         'FEV1FVC_SHRINE': 'FEV1/FVC', 'DEPRESSION_FG': 'Depression (FG)', 'T2D_XUE': 'Type 2 diabetes', 'HTN_FG': 'Hypertension (FG)',
         'HTN_PANUKB': 'Hypertension (Pan-UKB)'}
ORDER = ['SNORING_CAMPOS', 'INSOMNIA_FG', 'SLEEPDUR_DASHTI', 'CHRONOTYPE_JONES', 'FEV1FVC_SHRINE', 'DEPRESSION_FG', 'T2D_XUE', 'HTN_FG', 'HTN_PANUKB']
POS, NEG = '#9e2a55', '#3f7f55'
fig = plt.figure(figsize=(7.2, 4.6)); gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.05], wspace=0.75, left=0.2, right=0.97, top=0.9, bottom=0.12)
ax = fig.add_subplot(gs[0, 0])
bars = [('bivariate F_NonBMI', 'A6a/A6c bivariate F_NonBMI', t) for t in ORDER] + [('partial rg(OSA, X | BMI)', 'A6b partial rg(OSA, X | BMI)', t) for t in ('T2D_XUE', 'HTN_FG', 'HTN_PANUKB')]
ylab, yy, y = [], [], 0
for k, (sec, an, t) in enumerate(bars):
    if k == 0 or bars[k - 1][0] != sec:
        ax.text(-0.02, y + 0.1, 'Local rg with F$_{NonBMI}$' if sec.startswith('bivariate') else 'Partial local rg(OSA, X | BMI)',
                transform=ax.get_yaxis_transform(), ha='right', fontsize=6.6, fontweight='bold'); y -= 1
    r = S[(S.analysis == an) & (S.trait == t)]
    if not len(r): continue
    r = r.iloc[0]
    ax.barh(y, r.sig_pos, color=POS, height=0.62); ax.barh(y, -r.sig_neg, color=NEG, height=0.62)
    ax.text(max(r.sig_pos, 0) + 1.2, y, f'{int(r.sig_pos)}/{int(r.sig_neg)} of {int(r.loci_tested)}', va='center', fontsize=5.8, color=fk.INK2)
    ylab.append(NAMES[t]); yy.append(y); y -= 1
ax.set_yticks(yy); ax.set_yticklabels(ylab, fontsize=6.4); ax.tick_params(axis='y', length=0); ax.axvline(0, color=fk.INK2, lw=0.6)
ax.set_xlim(-8, 115); ax.set_xlabel('FDR-significant genomic blocks (negative | positive)', fontsize=6.8); ax.spines['left'].set_visible(False)
ax.set_title('a', loc='left', x=-0.68, fontsize=9, fontweight='bold')   # clear of the section header, flush with the left edge
bx = fig.add_subplot(gs[0, 1])
fb = st[(st.analysis == 'bivariate F_NonBMI') & st.phen2.isin(['SNORING_CAMPOS', 'FEV1FVC_SHRINE', 'SLEEPDUR_DASHTI'])].copy()
fb['o'] = fb.phen2.map({'SNORING_CAMPOS': 0, 'FEV1FVC_SHRINE': 1, 'SLEEPDUR_DASHTI': 2}); fb = fb.sort_values(['o', 'CHR', 'START'])
yb = np.arange(len(fb))[::-1]
for yi, rho, lo_, hi_ in zip(yb, fb['rho'], fb['rho.lower'], fb['rho.upper']):
    c = POS if rho > 0 else NEG
    bx.plot([lo_, hi_], [yi, yi], color=c, lw=1)
    bx.plot(rho, yi, 'D', color=c, ms=4)
bx.set_yticks(yb); bx.set_yticklabels([f'{NAMES[r.phen2]}  chr{int(r.CHR)}:{r.START / 1e6:.1f} Mb' for r in fb.itertuples()], fontsize=6)
bx.axvline(0, color=fk.INK2, lw=0.5); bx.set_xlim(-1.05, 1.05); bx.tick_params(axis='y', length=0); bx.spines['left'].set_visible(False)
bx.set_xlabel('Local rg with F$_{NonBMI}$ (95% CI)', fontsize=6.8)
bx.set_title('b', loc='left', x=-0.75, fontsize=9, fontweight='bold')
fig.text(0.02, 0.015, 'FDR q < 0.05 within each trait pair. FG, FinnGen (same registry as OSA). FEV1/FVC and sleep duration show significant blocks in both directions.',
         fontsize=5.6, color=fk.INK2)
fk.save_fig(fig, T + '/figures', 'SupplementaryFigure6_lava')
fk.write_raw_xlsx(f'{T}/figures/raw_data/SFig6_raw.xlsx', st, [f'{T}/figures/STable11_lava.tsv'], os.path.abspath(__file__),
                  {'rho': 'local rg', 'q_fdr': 'BH FDR within the trait pair'}, [('points plotted', len(st))])
print(V.to_string(index=False)); print(S.to_string(index=False))
