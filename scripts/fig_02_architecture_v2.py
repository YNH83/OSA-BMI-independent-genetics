#!/usr/bin/env python3
"""C9 Figure 2 (v2, following the internal figure design specification).
a  Partition of OSA SNP genetic variance (P1) for the primary model, sensitivity versions and the gout null control.
b  Table-style heatmap of genetic correlations of 17 external traits with OSA, F_BMI and F_NonBMI (Foote 2025 Fig 2 layout),
   grouped by domain; text weight encodes significance (bold = Bonferroni, black = P < 0.05, grey = not significant).
c  The same correlations with 95% CI, rows aligned to b (Zhou 2026 Fig 3 logic: original trait and both components side by side).
Raw data: figures/raw_data/Fig2v2_raw.xlsx"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from scipy.stats import norm
sys.path.insert(0, 'shared')
import figkit as fk

T = '.'
R = T + '/results'
C_OSA, C_BMI, C_NON = fk.INK2, fk.ORANGE, fk.BLUE
DOMAINS = [   # 2026-09-27: upper airway first (main axis)
    ('Sleep and upper airway', [('SNORING_CAMPOS', 'Snoring'), ('SNORING_CAMPOS_ADJBMI', 'Snoring (BMI-adj. GWAS)'),
                                 ('INSOMNIA_FG', 'Insomnia, FinnGen'), ('SLEEPDUR_DASHTI', 'Sleep duration'), ('CHRONOTYPE_JONES', 'Chronotype')]),
    ('Lower airway', [('FEV1FVC_SHRINE', 'FEV1/FVC'), ('ASTHMA_FG', 'Asthma, FinnGen')]),
    ('Neuropsychiatric', [('DEPRESSION_FG', 'Depression, FinnGen')]),
    ('Cardiometabolic', [('T2D_XUE', 'Type 2 diabetes'), ('T2D_FG', 'Type 2 diabetes, FinnGen'), ('HTN_FG', 'Hypertension, FinnGen'),
                         ('HTN_PANUKB', 'Hypertension, Pan-UKB'), ('SBP_EVANGELOU', 'Systolic BP (BMI-adj. GWAS)'), ('DBP_EVANGELOU', 'Diastolic BP (BMI-adj. GWAS)'),
                         ('SBP_PANUKB', 'Systolic BP, Pan-UKB'), ('DBP_PANUKB', 'Diastolic BP, Pan-UKB'),
                         ('CAD_ARAGAM', 'Coronary artery disease')]),
    ('Adiposity', [('WHRADJBMI_PULIT', 'WHR (BMI-adj. GWAS)')]),
    ('Socioeconomic and cognitive', [('INCOME_HILL', 'Household income'), ('COGNITION_LEE', 'Cognitive performance')]),
]
inputs, rows = [], []
for dom, items in DOMAINS:
    for tr, lab in items:
        pn, pe = f'{R}/02_numeric_{tr}.tsv', f'{R}/02_model_estimates_{tr}.tsv'
        inputs += [pn, pe]
        n = pd.read_csv(pn, sep='\t').iloc[0]
        e = pd.read_csv(pe, sep='\t'); e = e[(e.model == 'P2_with_T2D') & (e.op == '~~')]
        ext = lambda x: not str(x).startswith('F_')
        g = lambda f: e[((e.lhs == f) & e.rhs.map(ext) & (e.rhs != e.lhs)) | ((e.rhs == f) & e.lhs.map(ext) & (e.rhs != e.lhs))]
        b, nb = g('F_BMI'), g('F_NonBMI'); assert len(b) == 1 and len(nb) == 1, tr
        rows.append(dict(domain=dom.replace('\n', ' ').replace('- ', ''), trait_id=tr, trait=lab,
                         rg_OSA=n.rg_osa_ext, se_OSA=n.rg_osa_ext_se,
                         rg_F_BMI=float(b.STD_Genotype.iloc[0]), se_F_BMI=float(b.STD_Genotype_SE.iloc[0]),
                         rg_F_NonBMI=float(nb.STD_Genotype.iloc[0]), se_F_NonBMI=float(nb.STD_Genotype_SE.iloc[0])))
D = pd.DataFrame(rows)
NT = len(D); BONF = 0.05 / NT
for c in ('OSA', 'F_BMI', 'F_NonBMI'):
    D[f'P_{c}'] = 2 * norm.sf(np.abs(D[f'rg_{c}'] / D[f'se_{c}']))

# ---------- panel a data: variance partition ----------
PART = [('T2D_XUE', 'OSA, FinnGen (primary)'), ('T2D_XUE__OSA_FG_INCLAVO__BMI_GIANT', 'OSA incl. primary care'),
        ('T2D_XUE__OSA_CAMPOS__BMI_GIANT', 'OSA, Campos 2023*'), ('T2D_XUE__GOUT_FG__BMI_GIANT', 'Gout (null control)')]
A = []
for tag, lab in PART:
    p = f'{R}/02_numeric_{tag}.tsv'; inputs.append(p)
    r = pd.read_csv(p, sep='\t').iloc[0]
    A.append(dict(version=lab, P1=r.P1, P1_se=r.P1_se, share_BMI=1 - r.P1))
A = pd.DataFrame(A)

# ---------- figure ----------
# 2026-09-28 layout pass for journal artwork limits: final width <= 174 mm (figure 6.8 in wide, nothing drawn outside it),
# all text >= 6.5 pt at final size; panel letters placed in figure coordinates
FW, FH = 6.8, 9.05
fig = plt.figure(figsize=(FW, FH))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 5.2], width_ratios=[1.0, 1.12], hspace=0.17, wspace=0.07,
                      left=0.235, right=0.975, top=0.95, bottom=0.1)
ax = fig.add_subplot(gs[0, :])
ya = np.arange(len(A))[::-1]
for y, r in zip(ya, A.itertuples()):
    null = 'Gout' in r.version
    ax.barh(y, r.share_BMI, color=C_BMI, alpha=0.45 if null else 1, height=0.5, edgecolor='white', lw=1.2)
    ax.barh(y, r.P1, left=r.share_BMI, color=C_NON, alpha=0.45 if null else 1, height=0.5, edgecolor='white', lw=1.2)
    ax.errorbar(r.share_BMI, y - 0.42, xerr=1.96 * r.P1_se, color=fk.INK, lw=0.7, capsize=1.6, capthick=0.7)   # CI of the boundary = CI of P1
    ax.text(max(r.share_BMI / 2, 0.035), y, f'{r.share_BMI:.2f}', ha='center', va='center', fontsize=7, color='white', fontweight='bold')
    ax.text(r.share_BMI + r.P1 / 2 + 0.06, y, f'{r.P1:.2f}', ha='center', va='center', fontsize=7, color='white', fontweight='bold')
ax.axvline(0.40, color=fk.INK, lw=0.6, ls=(0, (2, 2)))
ax.text(0.412, len(A) - 0.35, 'P1 pass: boundary left of this line (BMI share ≤ 0.40)', fontsize=6.8, color=fk.INK2, va='bottom')
ax.set_yticks(ya); ax.set_yticklabels(A.version, fontsize=7)
ax.set_xlim(-0.32, 1); ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); ax.set_ylim(-0.75, len(A) - 0.1)   # room for whiskers that cross 0 (gout)
ax.spines['left'].set_bounds(-0.75, len(A) - 0.1)
ax.set_xlabel('Share of SNP genetic variance (whisker: 95% CI of the boundary)', fontsize=7)
ax.text(0.0, 1.15, 'via BMI (F$_{BMI}$)', color=C_BMI, fontsize=7.5, fontweight='bold', transform=ax.transAxes)
ax.text(1.0, 1.15, 'BMI-independent (F$_{NonBMI}$)', color=C_NON, fontsize=7.5, fontweight='bold', ha='right', transform=ax.transAxes)
fig.text(0.01, 0.985, 'a', fontsize=9, fontweight='bold', va='top')
ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)

# ---------- panel b: table heatmap ----------
bx = fig.add_subplot(gs[1, 0]); cx = fig.add_subplot(gs[1, 1], sharey=bx)
cmap = LinearSegmentedColormap.from_list('rg', ['#3f7f55', '#f7f6f3', '#9e2a55'])
ypos, yl, blocks, y = [], [], [], 0
for dom, items in DOMAINS:
    y += 0.95                     # room for the domain header row
    start = y
    for _ in items:
        ypos.append(y); y += 1
    blocks.append((dom, start, y - 1))
ypos = np.array(ypos)
cols = [('OSA', 'OSA', C_OSA), ('F_BMI', 'F$_{BMI}$', C_BMI), ('F_NonBMI', 'F$_{NonBMI}$', C_NON)]
for j, (k, lab, c) in enumerate(cols):
    bx.add_patch(Rectangle((j + 0.04, -0.72), 0.92, 0.7, color=c, clip_on=False))
    bx.text(j + 0.5, -0.37, lab, ha='center', va='center', color='white', fontsize=7.5, fontweight='bold')
    for yy, r in zip(ypos, D.itertuples()):
        v, p = getattr(r, f'rg_{k}'), getattr(r, f'P_{k}')
        bx.add_patch(Rectangle((j + 0.04, yy - 0.44), 0.92, 0.88, facecolor=cmap((np.clip(v, -0.7, 0.7) + 0.7) / 1.4), edgecolor='none'))
        txt = f'{v:.2f}'.replace('-', '−')
        col = fk.INK if p < 0.05 else '#9a9994'
        bx.text(j + 0.5, yy, txt, ha='center', va='center', fontsize=7, color='white' if abs(v) > 0.45 else col,
                fontweight='bold' if p < BONF else 'normal')
for dom, s, e in blocks:
    bx.text(-0.06, s - 0.82, dom, ha='right', va='center', fontsize=7, fontweight='bold', color=fk.INK, clip_on=False)
    bx.plot([-2.0, 3.0], [s - 0.55, s - 0.55], color='#c9c8c3', lw=0.5, clip_on=False)
bx.set_xlim(0, 3); bx.set_ylim(ypos.max() + 0.7, -0.75)
bx.set_yticks(ypos); bx.set_yticklabels(D.trait, fontsize=7)
bx.tick_params(axis='y', length=0, pad=4); bx.set_xticks([])
for s in bx.spines.values(): s.set_visible(False)
fig.text(0.01, bx.get_position().y1 + 0.004, 'b', fontsize=9, fontweight='bold', va='bottom')

# ---------- panel c: aligned CI plot ----------
off = {'OSA': -0.28, 'F_BMI': 0.0, 'F_NonBMI': 0.28}
for k, lab, c in cols:
    for yy, r in zip(ypos, D.itertuples()):
        v, s, p = getattr(r, f'rg_{k}'), getattr(r, f'se_{k}'), getattr(r, f'P_{k}')
        # same encoding as Figure 3a: filled = Bonferroni (0.05/NT), open = not; faded = P >= 0.05
        cx.plot([v - 1.96 * s, v + 1.96 * s], [yy + off[k]] * 2, color=c, lw=0.9, alpha=1 if p < 0.05 else 0.45, solid_capstyle='butt')
        cx.plot(v, yy + off[k], 'o', ms=3.6, color=c, mfc=c if p < BONF else 'white', mec=c, mew=0.8, alpha=1 if p < 0.05 else 0.6)
for _, s, e in blocks:
    cx.axhspan(s - 0.5, e + 0.5, color='#f3f2ef', zorder=-2, lw=0)
    cx.plot([-0.35, 0.75], [s - 0.55, s - 0.55], color='#c9c8c3', lw=0.5)
cx.axvline(0, color=fk.INK2, lw=0.6)
cx.set_xlim(-0.35, 0.75); cx.set_xlabel('Genetic correlation (95% CI)', fontsize=7)
cx.tick_params(axis='y', left=False, labelleft=False); cx.spines['left'].set_visible(False)
for k, lab, c in cols:
    cx.plot([], [], 'o', color=c, ms=3.6, label={'OSA': 'OSA', 'F_BMI': 'F$_{BMI}$', 'F_NonBMI': 'F$_{NonBMI}$'}[k])
cx.plot([], [], 'o', mfc='white', mec=fk.INK2, ms=3.6, label=f'P ≥ 0.05/{NT}')
cx.legend(loc='lower right', fontsize=6.8, ncol=2, handletextpad=0.2, columnspacing=0.8, borderaxespad=0.2)
fig.text(cx.get_position().x0 - 0.012, cx.get_position().y1 + 0.004, 'c', fontsize=9, fontweight='bold', va='bottom')
# colour bar for b
cb_ax = fig.add_axes([0.255, 0.058, 0.27, 0.008])
cb = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-0.7, 0.7)), cax=cb_ax, orientation='horizontal')
cb.set_ticks([-0.6, -0.3, 0, 0.3, 0.6]); cb.ax.tick_params(labelsize=6.8, length=2); cb.outline.set_visible(False)
cb.set_label(f'rg (bold: P < 0.05/{NT}; grey: P ≥ 0.05)', fontsize=7)
fig.text(0.02, 0.008, '*Campos 2023 includes FinnGen cases (consistency check, not independent). BMI-adj. GWAS are collider-prone for this contrast.',
         fontsize=6.5, color=fk.INK2)
out = fk.save_fig(fig, T + '/figures', 'Figure2v2_genetic_architecture')
fk.write_raw_xlsx(f'{T}/figures/raw_data/Fig2v2_raw.xlsx', {'b_c_rg': D, 'a_partition': A}, inputs, os.path.abspath(__file__),
                  {'rg_OSA': 'LDSC genetic correlation of OSA (FinnGen) with the trait', 'rg_F_BMI': 'Genomic SEM standardized covariance of F_BMI with the trait',
                   'rg_F_NonBMI': 'same for F_NonBMI', 'P_*': 'two-sided normal P from estimate / SE', 'P1': 'BMI-independent share of OSA SNP genetic variance'},
                  [(f'Bonferroni threshold 0.05/{NT}', BONF)] + [(f'{r.trait}: rg_OSA', round(r.rg_OSA, 4)) for r in D.head(3).itertuples()])
print(D[['trait', 'rg_OSA', 'rg_F_BMI', 'rg_F_NonBMI', 'P_F_NonBMI']].round(3).to_string(index=False)); print(out)
