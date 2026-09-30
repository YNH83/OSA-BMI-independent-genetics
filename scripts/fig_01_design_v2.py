#!/usr/bin/env python3
"""C9 Figure 1 (v2, following the internal figure design specification).
a  Concept: OSA SNP genetic variance = part via BMI + BMI-independent part (Zhou 2026 Fig 1 'subtraction equation').
b  Genomic SEM Cholesky model in the path-diagram grammar of Woodward 2026 Fig 3/4a: SNP box pointing to both factors,
   factors as filled circles in their fixed colours, coefficients in factor-coloured labels, line width proportional to
   the standardized estimate, fixed-zero covariance as a grey dashed arc; T2D shown as the example external trait (P2).
Raw data: figures/raw_data/Fig1v2_raw.xlsx"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch
sys.path.insert(0, 'shared')
import figkit as fk

T = '.'
R = T + '/results'
C_OSA, C_BMI, C_NON = fk.INK2, fk.ORANGE, fk.BLUE
pn = f'{R}/02_numeric_T2D_XUE.tsv'
n = pd.read_csv(pn, sep='\t').iloc[0]
lb_osa, ln_osa, p1 = n.s_bmi, n.s_non, n.P1
c_bmi_t2d, c_non_t2d = n.rg_bmi_ext, n.rg_non_ext

fig = plt.figure(figsize=(7.2, 6.4))   # 2026-09-28: minimum text size raised to 6.8 pt (journal artwork check); final width kept <= 174 mm
# ---------- a: concept bar ----------
ax = fig.add_axes([0.04, 0.855, 0.92, 0.12]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
ax.text(0.0, 0.5, 'OSA SNP\ngenetic variance', ha='left', va='center', fontsize=7.5, fontweight='bold', color=C_OSA)
x0, w = 0.2, 0.62
ax.text(x0 - 0.03, 0.5, '=', ha='center', va='center', fontsize=12, color=C_OSA)
ax.add_patch(FancyBboxPatch((x0, 0.28), w * (1 - p1), 0.44, boxstyle='round,pad=0,rounding_size=0.02', color=C_BMI, lw=0))
ax.add_patch(FancyBboxPatch((x0 + w * (1 - p1) + 0.004, 0.28), w * p1 - 0.004, 0.44, boxstyle='round,pad=0,rounding_size=0.02', color=C_NON, lw=0))
ax.text(x0 + w * (1 - p1) / 2, 0.5, f'via BMI\n{1 - p1:.0%}', ha='center', va='center', color='white', fontsize=7, fontweight='bold')
ax.text(x0 + w * (1 - p1) + w * p1 / 2, 0.5, f'BMI-independent (F$_{{NonBMI}}$)  {p1:.0%}', ha='center', va='center', color='white',
        fontsize=7.5, fontweight='bold')
ax.text(x0 + w + 0.015, 0.5, 'Share of LDSC SNP\nheritability (Cholesky)', ha='left', va='center', fontsize=7, color=fk.INK2)
ax.text(-0.0, 1.05, 'a', fontsize=9, fontweight='bold', transform=ax.transAxes)

# ---------- b: path diagram ----------
bx = fig.add_axes([0.16, 0.03, 0.68, 0.80]); bx.set_xlim(0, 10); bx.set_ylim(-0.75, 9.65); bx.axis('off'); bx.set_aspect('equal')
fig.text(0.04 + 0.0 * 0.92, 0.83, 'b', fontsize=9, fontweight='bold')


def label(x, y, s, c, fs=7):
    bx.text(x, y, s, ha='center', va='center', fontsize=fs, color='white', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.25,rounding_size=0.3', fc=c, ec='none'), zorder=6)


def arrow(p, q, c, w, ls='-', style='-|>'):
    bx.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=7, lw=w, color=c, ls=ls, shrinkA=0, shrinkB=0, zorder=2))


# SNP box
bx.add_patch(FancyBboxPatch((3.9, 8.7), 2.2, 0.8, boxstyle='round,pad=0.05,rounding_size=0.2', fc='#e9e8e4', ec=fk.INK2, lw=0.8))
bx.text(5.0, 9.1, 'each SNP', ha='center', va='center', fontsize=7, color=fk.INK)
FB, FN = (2.6, 6.2), (7.4, 6.2)
arrow((4.4, 8.65), (FB[0] + 0.45, FB[1] + 0.95), fk.INK2, 0.8)
arrow((5.6, 8.65), (FN[0] - 0.45, FN[1] + 0.95), fk.INK2, 0.8)
bx.text(5.0, 7.45, 'SNP effect on\neach factor', ha='center', va='center', fontsize=6.8, color=fk.INK2)
for (x, y), c, lab in [(FB, C_BMI, 'F$_{BMI}$'), (FN, C_NON, 'F$_{NonBMI}$')]:
    bx.add_patch(Circle((x, y), 1.0, fc=c, ec='white', lw=1.2, zorder=4))
    bx.text(x, y, lab, ha='center', va='center', fontsize=8, color='white', fontweight='bold', zorder=5)
# fixed-zero covariance
bx.add_patch(FancyArrowPatch((FB[0] + 1.0, FB[1] + 0.1), (FN[0] - 1.0, FN[1] + 0.1), arrowstyle='<|-|>', mutation_scale=6,
                             connectionstyle='arc3,rad=-0.25', ls=(0, (3, 2)), color=fk.GREY, lw=0.8))
bx.text(5.0, 6.05, 'fixed 0', ha='center', fontsize=6.8, color=fk.GREY)
# observed indicators
BMI, OSA, T2D = (2.6, 2.9), (7.4, 2.9), (5.0, 0.45)
for (x, y), lab in [(BMI, 'BMI\n(GIANT)'), (OSA, 'OSA\n(FinnGen R12)')]:
    bx.add_patch(FancyBboxPatch((x - 1.05, y - 0.55), 2.1, 1.1, boxstyle='round,pad=0.02,rounding_size=0.25', fc='white', ec=fk.INK, lw=1.0, zorder=4))
    bx.text(x, y, lab, ha='center', va='center', fontsize=7, color=fk.INK, zorder=5)
W = lambda v: 0.6 + 3.2 * abs(v)     # line width proportional to the standardized estimate
arrow((FB[0], FB[1] - 1.0), (BMI[0], BMI[1] + 0.58), C_BMI, W(1.0))
arrow((FB[0] + 0.7, FB[1] - 0.72), (OSA[0] - 0.8, OSA[1] + 0.58), C_BMI, W(lb_osa))
arrow((FN[0], FN[1] - 1.0), (OSA[0], OSA[1] + 0.58), C_NON, W(ln_osa))
label(2.1, 4.55, '1', C_BMI); label(5.4, 4.9, f'{lb_osa:.2f}', C_BMI); label(7.95, 4.55, f'{ln_osa:.2f}', C_NON)
# external trait (P2 example): covariances with both factors
bx.add_patch(FancyBboxPatch((T2D[0] - 1.0, T2D[1] - 0.42), 2.0, 0.84, boxstyle='round,pad=0.02,rounding_size=0.25', fc='#f3f2ef', ec=fk.INK2, lw=0.8, zorder=4))
bx.text(T2D[0], T2D[1], 'T2D (external)', ha='center', va='center', fontsize=6.8, color=fk.INK, zorder=5)
bx.add_patch(FancyArrowPatch((FB[0] - 0.98, FB[1] - 0.15), (T2D[0] - 1.02, T2D[1]), arrowstyle='<|-|>', mutation_scale=6,
                             connectionstyle='arc3,rad=0.62', color=C_BMI, lw=W(c_bmi_t2d) * 0.7, alpha=0.85, zorder=1))
bx.add_patch(FancyArrowPatch((FN[0] + 0.98, FN[1] - 0.15), (T2D[0] + 1.02, T2D[1]), arrowstyle='<|-|>', mutation_scale=6,
                             connectionstyle='arc3,rad=-0.62', color=C_NON, lw=W(c_non_t2d) * 0.7, alpha=0.85, zorder=1))
label(0.55, 1.3, f'{c_bmi_t2d:.2f}', C_BMI); label(9.45, 1.3, f'{c_non_t2d:.2f}', C_NON)
bx.text(5.0, -0.45, f'P1 = {ln_osa:.2f}² / ({lb_osa:.2f}² + {ln_osa:.2f}²) = {p1:.2f}', ha='center', fontsize=7, color=fk.INK2)

out = fk.save_fig(fig, T + '/figures', 'Figure1v2_design')
raw = pd.DataFrame([dict(parameter='F_BMI =~ OSA (standardized)', value=lb_osa), dict(parameter='F_NonBMI =~ OSA (standardized)', value=ln_osa),
                    dict(parameter='F_BMI =~ BMI (standardized)', value=1.0), dict(parameter='P1', value=p1),
                    dict(parameter='F_BMI ~~ T2D (standardized)', value=c_bmi_t2d), dict(parameter='F_NonBMI ~~ T2D (standardized)', value=c_non_t2d)])
fk.write_raw_xlsx(f'{T}/figures/raw_data/Fig1v2_raw.xlsx', raw, [pn], os.path.abspath(__file__),
                  {'value': 'standardized Genomic SEM estimate from the primary model (FinnGen OSA, GIANT BMI, T2D Xue)'},
                  [('P1 = ln^2 / (lb^2 + ln^2)', ln_osa ** 2 / (lb_osa ** 2 + ln_osa ** 2))])
print(raw.round(4).to_string(index=False)); print(out)
