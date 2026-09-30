#!/usr/bin/env python3
"""C9 Figure 3 (v2; 2026-09-27 main axis = upper airway; hypertension tests moved to Supplementary Figure 7; following the internal figure design specification): where traits sit on the two OSA components.
a  Axis map: rg with F_BMI (x) vs rg with F_NonBMI (y), cross 95% CI, y = x line, background tint by which component
   dominates; filled = rg with F_NonBMI Bonferroni-significant, hollow = not; square = FinnGen registry trait.
b  Prespecified two-axis equality tests (A1, A1b, A2) and the post hoc same-registry test: dumbbell from rg(F_NonBMI, T2D)
   to rg(F_NonBMI, second trait), text column with the difference and P.
c  Bidirectional MR between BMI and OSA (A3), both directions, shared method rows.
Raw data: figures/raw_data/Fig3v2_raw.xlsx"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
sys.path.insert(0, 'shared')
import figkit as fk

T = '.'
R = T + '/results'
C_BMI, C_NON = fk.ORANGE, fk.BLUE
MINUS = '−'
fmt = lambda v, d=2: f'{v:.{d}f}'.replace('-', MINUS)

# ---------- data a: reuse Figure 2 v2 raw (same estimates, one source of truth) ----------
raw2 = f'{T}/figures/raw_data/Fig2v2_raw.xlsx'
D = pd.read_excel(raw2, sheet_name='b_c_rg')
NT = len(D); BONF = 0.05 / NT
SHORT = {'HTN_PANUKB': 'Hypertension (UKB)', 'SBP_PANUKB': 'SBP (UKB)', 'DBP_PANUKB': 'DBP (UKB)', 'T2D_XUE': 'T2D', 'T2D_FG': 'T2D (FG)', 'HTN_FG': 'Hypertension (FG)', 'SBP_EVANGELOU': 'SBP*', 'DBP_EVANGELOU': 'DBP*',
         'CAD_ARAGAM': 'CAD', 'WHRADJBMI_PULIT': 'WHR*', 'SNORING_CAMPOS': 'Snoring', 'SNORING_CAMPOS_ADJBMI': 'Snoring*',
         'INSOMNIA_FG': 'Insomnia (FG)', 'SLEEPDUR_DASHTI': 'Sleep duration', 'CHRONOTYPE_JONES': 'Chronotype', 'FEV1FVC_SHRINE': 'FEV1/FVC',
         'ASTHMA_FG': 'Asthma (FG)', 'DEPRESSION_FG': 'Depression (FG)', 'INCOME_HILL': 'Income', 'COGNITION_LEE': 'Cognition'}
# domain hues avoid the component colours (orange F_BMI, blue F_NonBMI); airway domains share a purple family.
# dataviz validator 2026-09-27: CVD worst adjacent dE 9.9, normal-vision 17.1; grey is the intentional neutral; every point is direct-labelled
DOMCOL = {'Cardiometabolic': '#8c510a', 'Adiposity': fk.GREY, 'Sleep and upper airway': '#6a3d9a', 'Lower airway': '#a58bd1',
          'Neuropsychiatric': '#1f9e89', 'Socioeconomic and cognitive': fk.YELLOW}
# label offsets in points, set after visual check
OFF = {'DBP_PANUKB': (6, -16), 'T2D_XUE': (6, -9), 'T2D_FG': (6, 2), 'HTN_FG': (6, 4), 'CAD_ARAGAM': (6, -3), 'WHRADJBMI_PULIT': (-6, 5),
       'SNORING_CAMPOS': (-6, 4), 'SNORING_CAMPOS_ADJBMI': (6, 3)}
# labels placed in data coords (x, y, ha, leader); 2026-09-27 layout pass: labels of the crowded points near rg = 0 sit in the empty
# band below the x axis, joined by thin leaders that leave each point at an angle, so no label crosses a CI, the y = x line or a zero line
LEADER = {'HTN_PANUKB': (0.71, -0.21, 'right', True, (0.7, 1.0)), 'SBP_EVANGELOU': (-0.27, 0.03, 'left', True), 'DBP_EVANGELOU': (-0.27, 0.09, 'left', True),
          'SBP_PANUKB': (0.26, -0.12, 'left', True), 'FEV1FVC_SHRINE': (0.19, -0.165, 'left', True), 'CHRONOTYPE_JONES': (0.13, -0.21, 'left', True),
          'SLEEPDUR_DASHTI': (0.07, -0.255, 'left', True), 'COGNITION_LEE': (0.02, -0.30, 'left', True), 'INCOME_HILL': (0.012, -0.345, 'left', True),
          'DEPRESSION_FG': (0.02, 0.43, 'left', True), 'ASTHMA_FG': (0.24, 0.38, 'center', False), 'INSOMNIA_FG': (-0.06, 0.45, 'right', False)}
# 2026-09-28 (6.8 pt labels): band labels shifted left and the Hypertension (UKB) label moved to the right end of the -0.21 row so no two labels touch

# 2026-09-28 layout pass for journal artwork limits: every axes is placed in inches on a 6.8 in (172.7 mm) wide canvas and
# nothing is drawn outside it, so the final width stays <= 174 mm; all text >= 6.5 pt at final size. Panel c is stacked
# (BMI -> OSA above OSA -> BMI) so that the estimate and P columns sit in their own axes at a legible size.
FW, FH = 6.8, 8.15
fig = plt.figure(figsize=(FW, FH))
A_ = lambda x, y, w, h, **kw: fig.add_axes([x / FW, y / FH, w / FW, h / FH], **kw)   # inches from the bottom-left corner
FS_LAB, FS_TXT = 6.8, 7   # point labels in a; tick labels and numeric columns
lim, ylo = (-0.3, 0.72), -0.42    # y axis extended below to hold the label band
AW = 2.64; AH = AW * (lim[1] - ylo) / (lim[1] - lim[0]); A_TOP = FH - 0.22
ax = A_(0.6, A_TOP - AH, AW, AH)
ax.add_patch(Polygon([[lim[0], lim[0]], [lim[0], lim[1]], [lim[1], lim[1]]], color=C_NON, alpha=0.07, lw=0))
ax.add_patch(Polygon([[lim[0], ylo], [lim[1], ylo], [lim[1], lim[1]], [lim[0], lim[0]]], color=C_BMI, alpha=0.07, lw=0))
ax.plot(lim, lim, color=fk.INK2, lw=0.6, ls=(0, (3, 2)))
ax.axhline(0, color=fk.GREY, lw=0.5); ax.axvline(0, color=fk.GREY, lw=0.5)
ax.text(0.02, 0.69, 'closer to F$_{NonBMI}$', color=C_NON, fontsize=7, fontweight='bold', va='top')
ax.text(0.69, -0.405, 'closer to F$_{BMI}$', color=C_BMI, fontsize=7, fontweight='bold', ha='right')
for r in D.itertuples():
    c = DOMCOL[r.domain]; sig = r.P_F_NonBMI < BONF
    ax.errorbar(r.rg_F_BMI, r.rg_F_NonBMI, xerr=1.96 * r.se_F_BMI, yerr=1.96 * r.se_F_NonBMI, fmt='none', ecolor=c, elinewidth=0.6, alpha=0.55)
    ax.plot(r.rg_F_BMI, r.rg_F_NonBMI, 's' if r.trait_id.endswith('_FG') else 'o', ms=5.2, mfc=c if sig else 'white', mec=c, mew=1.1, zorder=3)
    if r.trait_id in LEADER:
        lx, ly, lha, lead = LEADER[r.trait_id][:4]
        rp = LEADER[r.trait_id][4] if len(LEADER[r.trait_id]) > 4 else ((0, 0.5) if lha == 'left' else (0.5, 0.5))
        ax.annotate(SHORT[r.trait_id], (r.rg_F_BMI, r.rg_F_NonBMI), xytext=(lx, ly), textcoords='data', fontsize=FS_LAB, ha=lha,
                    va='center', color=fk.INK, arrowprops=dict(arrowstyle='-', lw=0.4, color=fk.INK2, shrinkA=1, shrinkB=3, relpos=rp) if lead else None)
        continue
    dx, dy = OFF.get(r.trait_id, (5, 3))
    ax.annotate(SHORT[r.trait_id], (r.rg_F_BMI, r.rg_F_NonBMI), xytext=(dx, dy), textcoords='offset points', fontsize=FS_LAB,
                ha='left' if dx > 0 else ('right' if dx < 0 else 'center'), color=fk.INK)
ax.set_xlim(lim); ax.set_ylim(ylo, lim[1]); ax.set_aspect('equal')
ax.set_xticks([-0.2, 0, 0.2, 0.4, 0.6]); ax.set_yticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6])
ax.tick_params(labelsize=FS_TXT)
ax.set_xlabel('rg with F$_{BMI}$ (95% CI)'); ax.set_ylabel('rg with F$_{NonBMI}$ (95% CI)', labelpad=2)
for dom, c in DOMCOL.items():
    ax.plot([], [], 'o', color=c, ms=4.5, label=dom)
ax.plot([], [], 'o', mfc='white', mec=fk.INK2, ms=4.5, label=f'F$_{{NonBMI}}$ P ≥ 0.05/{NT}')
ax.plot([], [], 's', color=fk.INK2, ms=4.5, label='FinnGen (same registry as OSA)')
_h, _l = ax.get_legend_handles_labels()
LEG_TOP = A_TOP - AH - 0.42
fig.legend(_h, _l, loc='upper left', bbox_to_anchor=(0.06 / FW, LEG_TOP / FH), ncol=4, fontsize=6.8, handletextpad=0.2, columnspacing=0.9,
           labelspacing=0.35, borderaxespad=0, frameon=False)
fig.text(0.06 / FW, (LEG_TOP - 0.42) / FH, '*GWAS adjusted for BMI (collider-prone for this contrast); FG, FinnGen; UKB, Pan-UKB (not BMI-adjusted)',
         fontsize=6.5, color=fk.INK2, va='top')
fig.text(0.02 / FW, A_TOP / FH + 0.004, 'a', fontsize=9, fontweight='bold', va='bottom')

# ---------- b: two-axis tests ----------
TESTS = [('A1_HTN', 'Hypertension (FG)', 'A1'), ('A1b_CAD', 'CAD', 'A1b'), ('A2_SBP', 'SBP*', 'A2'), ('A2_DBP', 'DBP*', 'A2'),
         ('posthoc_HTN_vs_T2DFG', 'Hypertension vs\nT2D, both FG', 'post hoc'),
         ('A7_HTN_PANUKB', 'Hypertension (UKB)', 'A7'), ('A7_SBP_PANUKB', 'SBP (UKB)', 'A7'), ('A7_DBP_PANUKB', 'DBP (UKB)', 'A7')]
B, inputs = [], [raw2]
for lab, name, tag in TESTS:
    p = f'{R}/07_twoaxis_{lab}.tsv'; inputs.append(p)
    r = pd.read_csv(p, sep='\t').iloc[0]
    B.append(dict(test=lab, label=name, prereg=tag, metab=r.metab, cor_T2D=r.cor_NonBMI_metab, se_T2D=r.se_metab,
                  cor_axis2=r.cor_NonBMI_axis2, se_axis2=r.se_axis2, diff=r.cor_NonBMI_axis2 - r.cor_NonBMI_metab,
                  chisq=r.constrained_chisq, P=r.p_difference, grade=r.grade))
B = pd.DataFrame(B)
# ---------- b (Figure 3): upper-airway specificity of F_NonBMI ----------
SPEC = [('Upper airway', [('SNORING_CAMPOS', 'Snoring'), ('SNORING_CAMPOS_ADJBMI', 'Snoring, BMI-adj.\nGWAS (sensitivity)')]),
        ('Sleep', [('INSOMNIA_FG', 'Insomnia (FG)'), ('SLEEPDUR_DASHTI', 'Sleep duration'), ('CHRONOTYPE_JONES', 'Chronotype')]),
        ('Mood', [('DEPRESSION_FG', 'Depression (FG)')]),
        ('Lower airway', [('FEV1FVC_SHRINE', 'FEV1/FVC'), ('ASTHMA_FG', 'Asthma (FG)')])]
BX0 = 4.62
bx = A_(BX0, A_TOP - AH, FW - 0.22 - BX0, AH - 0.34)
BONF_A5 = 0.05 / 7   # panel b is the preregistered A5 mechanism panel: Bonferroni over its 7 traits (Stage 2 protocol), not over all 20
Dd = D.set_index('trait_id'); Sp, ys, y0, heads = [], [], 0.0, []
for grp_, items in SPEC:
    heads.append((grp_, y0)); y0 -= 0.7
    for tid, lab in items:
        r = Dd.loc[tid]
        for k, c, off, mk in (('F_BMI', C_BMI, 0.17, 'o'), ('F_NonBMI', C_NON, -0.17, 'D')):
            v, se_, pv = r[f'rg_{k}'], r[f'se_{k}'], r[f'P_{k}']
            sig = pv < BONF_A5
            bx.plot([v - 1.96 * se_, v + 1.96 * se_], [y0 + off] * 2, color=c, lw=1, alpha=1 if sig else 0.45)
            bx.plot(v, y0 + off, mk, ms=4.3, mfc=c if sig else 'white', mec=c, mew=0.9)
            Sp.append(dict(group=grp_, trait_id=tid, trait=lab.replace('\n', ' '), component=k, rg=v, se=se_, P=pv, bonferroni=sig))
        ys.append((y0, lab)); y0 -= 1.0
    y0 -= 0.3
for grp_, yy_ in heads:
    bx.text(-0.03, yy_, grp_, transform=bx.get_yaxis_transform(), ha='right', va='center', fontsize=7, fontweight='bold', color=fk.INK)
bx.set_yticks([y for y, _ in ys]); bx.set_yticklabels([l for _, l in ys], fontsize=FS_TXT); bx.tick_params(axis='y', length=0)
bx.tick_params(axis='x', labelsize=FS_TXT)
bx.axvline(0, color=fk.INK2, lw=0.5); bx.set_xlim(-0.3, 0.72); bx.set_ylim(y0 + 0.4, 0.5)
bx.spines['left'].set_visible(False); bx.set_xlabel('Genetic correlation (95% CI)')
bx.plot([], [], 'o', color=C_BMI, ms=4, label='F$_{BMI}$'); bx.plot([], [], 'D', color=C_NON, ms=4, label='F$_{NonBMI}$')
bx.plot([], [], 'o', mfc='white', mec=fk.INK2, ms=4, label='P ≥ 0.05/7')
bx.legend(loc='lower right', bbox_to_anchor=(1.0, 1.0), ncol=3, fontsize=6.8, handletextpad=0.2, columnspacing=0.8, borderaxespad=0.3)
fig.text(3.42 / FW, A_TOP / FH + 0.004, 'b', fontsize=9, fontweight='bold', va='bottom')
Sp = pd.DataFrame(Sp)

def pfig(p):
    """P for figure text: a x 10^-b below 0.001."""
    if p >= 0.001: return f'{p:.3g}' if p < 0.1 else f'{p:.2f}'   # 3 significant figures below 0.1 so 0.0497 is not shown as 0.05
    m, e = f'{p:.1e}'.split('e'); return f'{m}×10$^{{{int(e)}}}$'.replace('-', chr(0x2212))


# ---------- c: bidirectional MR ----------
mr = pd.read_csv(f'{R}/08_mr_A3.tsv', sep='\t'); inputs.append(f'{R}/08_mr_A3.tsv')
METH = [('Inverse variance weighted', 'IVW'), ('Weighted median', 'Weighted median'), ('MR Egger', 'MR-Egger'), ('Weighted mode', 'Weighted mode')]
C = []
for d, dlab, col in [('BMI_to_OSA', 'BMI → OSA', C_BMI), ('OSA_to_BMI', 'OSA → BMI', fk.INK2)]:
    x = mr[mr.direction == d].set_index('method')
    for m_, ml in METH:
        C.append(dict(direction=d, dlab=dlab, color=col, method=ml, nsnp=int(x.loc[m_, 'nsnp']), b=x.loc[m_, 'b'], lo=x.loc[m_, 'lo'],
                      hi=x.loc[m_, 'hi'], P=x.loc[m_, 'pval']))
    f0 = x.iloc[0]; bs, ss_ = f0.ivw_steiger_filtered_b, f0.ivw_steiger_filtered_se
    from scipy.stats import norm as _n
    C.append(dict(direction=d, dlab=dlab, color=col, method=f'IVW, Steiger-concordant ({int(f0.ivw_steiger_filtered_n)})', nsnp=int(f0.ivw_steiger_filtered_n),
                  b=bs, lo=bs - 1.96 * ss_, hi=bs + 1.96 * ss_, P=2 * _n.sf(abs(bs / ss_))))
C = pd.DataFrame(C)
NIV = {d: int(C[(C.direction == d) & (C.method == 'IVW')].nsnp.iloc[0]) for d in ('BMI_to_OSA', 'OSA_to_BMI')}
CX0, CXW, EX0, EXW, PX0, PXW = 1.62, 2.45, 4.2, 1.35, 5.62, 0.9   # forest, estimate column, P column (inches)
C_TOP = (LEG_TOP - 0.42) - 0.62; C_AH, C_GAP = 1.05, 0.62
for i, (d, lab, col, xl) in enumerate([('BMI_to_OSA', f'BMI → OSA ({NIV["BMI_to_OSA"]} SNPs)', C_BMI, 'log odds of OSA per SD BMI'),
                                       ('OSA_to_BMI', f'OSA → BMI ({NIV["OSA_to_BMI"]} SNPs)', fk.INK2, 'SD BMI per log odds of OSA')]):
    top_ = C_TOP - i * (C_AH + C_GAP)
    cx = A_(CX0, top_ - C_AH, CXW, C_AH)
    ex = A_(EX0, top_ - C_AH, EXW, C_AH, sharey=cx); px = A_(PX0, top_ - C_AH, PXW, C_AH, sharey=cx)
    cc = C[C.direction == d].reset_index(drop=True)
    yy = np.arange(len(cc))[::-1]
    for y, r in zip(yy, cc.itertuples()):
        sig = r.P < 0.05
        cx.plot([r.lo, r.hi], [y, y], color=col, lw=1, alpha=1 if sig else 0.45)
        cx.plot(r.b, y, 'D' if r.method.startswith('IVW') else 'o', ms=4.8 if r.method == 'IVW' else 3.9, mfc=col if sig else 'white', mec=col, mew=0.9)
        ex.text(0, y, f'{fmt(r.b)} ({fmt(r.lo)}{chr(0x2013)}{fmt(r.hi)})', va='center', fontsize=FS_TXT, color=fk.INK if sig else fk.GREY)
        px.text(0, y, pfig(r.P), va='center', fontsize=FS_TXT, color=fk.INK if sig else fk.GREY)
    for t_ in (ex, px):
        t_.set_xlim(0, 1); t_.axis('off')
    ex.text(0, len(cc) - 0.35, 'Estimate (95% CI)', fontsize=FS_TXT, fontweight='bold', color=fk.INK2, va='bottom')
    px.text(0, len(cc) - 0.35, 'P', fontsize=FS_TXT, fontweight='bold', color=fk.INK2, va='bottom')
    cx.set_yticks(yy); cx.set_yticklabels(cc.method, fontsize=FS_TXT); cx.tick_params(axis='y', length=0); cx.tick_params(axis='x', labelsize=FS_TXT)
    cx.axvline(0, color=fk.INK2, lw=0.5); cx.set_ylim(-0.6, len(cc) - 0.4)
    lo_, hi_ = min(0, cc.lo.min()), cc.hi.max()
    cx.set_xlim(lo_ - 0.05 * (hi_ - lo_), hi_ + 0.05 * (hi_ - lo_))
    cx.set_xlabel(xl, fontsize=7); cx.spines['left'].set_visible(False)
    cx.text(0.12 / FW, top_ / FH + 0.004, lab, transform=fig.transFigure, fontsize=7.5, fontweight='bold', color=col, va='bottom')
fig.text(0.02 / FW, (C_TOP + 0.22) / FH + 0.004, 'c', fontsize=9, fontweight='bold', va='bottom')
out = fk.save_fig(fig, T + '/figures', 'Figure3v2_upper_airway_axis')
fk.write_raw_xlsx(f'{T}/figures/raw_data/Fig3v2_raw.xlsx', {'a_axis_map': D, 'b_upper_airway': Sp, 'c_mr': C.drop(columns=['color'])}, inputs,
                  os.path.abspath(__file__),
                  {'rg_F_BMI': 'Genomic SEM correlation with F_BMI', 'rg_F_NonBMI': 'correlation with F_NonBMI',
                   'diff': 'rg(F_NonBMI, second trait) minus rg(F_NonBMI, T2D), both from the joint two-axis model',
                   'P (b)': '1-df chi-square P of the equality constraint', 'b': 'MR estimate'},
                  [(f'panel a Bonferroni threshold 0.05/{NT}', BONF), ('panel b Bonferroni threshold 0.05/7 (prespecified A5 panel)', BONF_A5)] + [(f'{r.test} diff', round(r.diff, 4)) for r in B.itertuples()])
print(out)
# ---------- Supplementary Figure 7: hypertension-related two-axis tests (A1, A1b, A2, A7; post hoc) ----------
fig, bx = plt.subplots(figsize=(5.8, 4.1))   # 2026-09-28: text raised to >= 6.5 pt
fig.subplots_adjust(left=0.34, right=0.97, top=0.86, bottom=0.2)
yb = np.arange(len(B))[::-1]
LIGHT = '#9cc3ee'
for y, r in zip(yb, B.itertuples()):
    sig = str(r.grade).startswith('PASS') and r.prereg != 'post hoc'   # preregistered grade (direction AND P); the post hoc row is never styled as a rule met
    bx.plot([r.cor_T2D, r.cor_axis2], [y, y], color=C_NON if sig else '#b9c9dc', lw=2.2 if sig else 1.4, solid_capstyle='round', zorder=1)
    bx.plot(r.cor_T2D, y, 'o', ms=5.5, mfc=LIGHT, mec=C_NON, mew=0.8, zorder=2)
    bx.plot(r.cor_axis2, y, 'D', ms=5, mfc=C_NON if sig else 'white', mec=C_NON, mew=1, zorder=2)
    bx.text(0.43, y, f'{"+" if r.diff >= 0 else MINUS}{abs(r.diff):.2f}', va='center', fontsize=7, color=fk.INK,
            fontweight='bold' if sig else 'normal')
    bx.text(0.56, y, pfig(r.P), va='center', fontsize=7,
            color=fk.INK if sig else fk.GREY, fontweight='bold' if sig else 'normal')
bx.text(0.43, len(B) - 0.35, 'Δ rg', fontsize=7, fontweight='bold', color=fk.INK2)
bx.text(0.56, len(B) - 0.35, 'P', fontsize=7, fontweight='bold', color=fk.INK2)
bx.set_yticks(yb); bx.set_yticklabels([r.label + ('  (post hoc)' if 'post' in str(r.prereg).lower() else '') for r in B.itertuples()], fontsize=7)
bx.set_xlim(-0.02, 0.66); bx.set_ylim(-0.6, len(B) - 0.1)
bx.set_xticks([0, 0.1, 0.2, 0.3]); bx.spines['left'].set_visible(False); bx.tick_params(axis='y', length=0)
bx.spines['bottom'].set_bounds(-0.02, 0.32)   # axis line only under the data, not under the text columns
bx.set_xlabel('rg with F$_{NonBMI}$'); bx.xaxis.set_label_coords((0.15 + 0.02) / 0.68, -0.09); bx.tick_params(axis='x', labelsize=7)
bx.plot([], [], 'o', mfc=LIGHT, mec=C_NON, label='with T2D'); bx.plot([], [], 'D', color=C_NON, label='with second trait, rule met')
bx.plot([], [], 'D', mfc='white', mec=C_NON, label='rule not met')
bx.legend(loc='lower left', bbox_to_anchor=(-0.02, 1.0), ncol=3, fontsize=6.8, handletextpad=0.2)
bx.set_title('Correlation of F$_{NonBMI}$ with T2D versus a second trait', loc='left', x=-0.62, fontsize=7.5, fontweight='bold', y=1.08)


fig.text(0.02, 0.015, '*GWAS adjusted for BMI (Evangelou 2018). FG, FinnGen (same registry as OSA); UKB, Pan-UKB (not BMI-adjusted).\nBold and filled = prespecified rule met; the post hoc row is descriptive only.', fontsize=6.5, color=fk.INK2)
fk.save_fig(fig, T + '/figures', 'SupplementaryFigure7_hypertension_tests')
fk.write_raw_xlsx(f'{T}/figures/raw_data/SFig7_raw.xlsx', B, [f'{R}/07_twoaxis_{t}.tsv' for t, _, _ in TESTS], os.path.abspath(__file__),
                  {'cor_T2D': 'rg(F_NonBMI, T2D) in the joint model', 'cor_axis2': 'rg(F_NonBMI, second trait)',
                   'diff': 'second trait minus T2D', 'P': '1-df chi-square P of the equality constraint', 'grade': 'prespecified rule (direction and P)'},
                  [(f'{r.test} diff', round(r.diff, 4)) for r in B.itertuples()],
                  notes=['A2 used BMI-adjusted GWAS (Evangelou 2018); A7 used Pan-UKB GWAS not adjusted for BMI and outside FinnGen. Post hoc row not prespecified.'])

