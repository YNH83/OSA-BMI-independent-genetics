#!/usr/bin/env python3
"""C9 Table 2: preregistered thresholds versus observed results (kill-switch, P1, P2, P3),
with raw-data Excel figures/raw_data/Table2_raw.xlsx. Missing results are shown as pending."""
import os, sys
import pandas as pd
sys.path.insert(0, 'shared')
import figkit

T = '.'
R = T + '/results'
PRE = 'analysis_plan/C9_預先登錄_Stage1.md'
inputs = [PRE]
rows = []

ks = pd.read_csv(f'{R}/01_killswitch_summary.tsv', sep='\t').set_index('metric')['value']; inputs.append(f'{R}/01_killswitch_summary.tsv')
rows.append(dict(item='Kill-switch 1', metric='OSA liability SNP heritability z', threshold='retire if z < 5',
                 observed=float(ks['OSA h2 z']), ci='', grade='PASS' if float(ks['OSA h2 z']) >= 5 else 'RETIRE'))
rg = float(ks['rg(OSA,BMI)']); rgse = float(ks['rg SE'])
rows.append(dict(item='Kill-switch 2', metric='rg(OSA, BMI)', threshold='retire if > 0.85', observed=rg,
                 ci=f'{rg - 1.96 * rgse:.3f} to {rg + 1.96 * rgse:.3f}', grade='PASS' if rg <= 0.85 else 'RETIRE'))

n = pd.read_csv(f'{R}/02_numeric_T2D_XUE.tsv', sep='\t').iloc[0]; inputs.append(f'{R}/02_numeric_T2D_XUE.tsv')
g = lambda x, a, b: 'PASS' if x >= a else ('PARTIAL' if x >= b else 'FAIL')
rows.append(dict(item='P1', metric='BMI-independent share of OSA SNP genetic variance', threshold='PASS >= 0.60; PARTIAL 0.40 to 0.60; FAIL < 0.40',
                 observed=n.P1, ci=f'{n.P1 - 1.96 * n.P1_se:.3f} to {n.P1 + 1.96 * n.P1_se:.3f}', grade=g(n.P1, 0.60, 0.40)))
rows.append(dict(item='P2', metric='Attenuation of rg with type 2 diabetes, 1 - rg(NonBMI,T2D)/rg(OSA,T2D)',
                 threshold='PASS >= 0.50; PARTIAL 0.25 to 0.50; FAIL < 0.25', observed=n.P2,
                 ci=f'{n.P2 - 1.96 * n.P2_se:.3f} to {n.P2 + 1.96 * n.P2_se:.3f} (delta method, independence approximation)', grade=g(n.P2, 0.50, 0.25)))
p3f = f'{R}/04_P3_summary_all_fixP.tsv'
if os.path.exists(p3f):
    p3 = pd.read_csv(p3f, sep='\t').set_index('item')['value']; inputs.append(p3f)
    v = float(p3['P3 share retained'])
    rows.append(dict(item='P3', metric='Share of FinnGen OSA loci retained at P < 5e-8 in NonBMI-OSA GWAS',
                     threshold='PASS >= 0.30; PARTIAL 0.15 to 0.30; FAIL < 0.15', observed=v,
                     ci=f"{p3['... retained in NonBMI-OSA (any GWS SNP within 500 kb)']} of {p3['FinnGen OSA loci (distance clumping, 500 kb)']} loci", grade=g(v, 0.30, 0.15)))
else:
    rows.append(dict(item='P3', metric='Share of FinnGen OSA loci retained at P < 5e-8 in NonBMI-OSA GWAS',
                     threshold='PASS >= 0.30; PARTIAL 0.15 to 0.30; FAIL < 0.15', observed=None, ci='', grade='pending (GWAS running)'))
df = pd.DataFrame(rows)
checks = [('P1 = (std loading NonBMI)^2', f'{n.s_non:.6f}^2 = {n.s_non ** 2:.6f}'),
          ('P2 = 1 - rg_NonBMI / rg_OSA', f'1 - {n.rg_non_ext:.6f} / {n.rg_osa_ext:.6f} = {1 - n.rg_non_ext / n.rg_osa_ext:.6f}'),
          ('model consistency: s_bmi*rg_bmi + s_non*rg_non = rg_OSA', f'{n.s_bmi:.4f}*{n.rg_bmi_ext:.4f} + {n.s_non:.4f}*{n.rg_non_ext:.4f} = {n.check_implied_rg:.4f} (LDSC {n.rg_osa_ext:.4f})')]
cols = {'threshold': 'copied verbatim from analysis_plan/C9_預先登錄_Stage1.md (SHA256 and time stamp in analysis_plan/SHA256_登錄時間戳.txt)',
        'observed': 'point estimate', 'ci': '95% confidence interval, or counts for P3', 'grade': 'decision under the preregistered rule'}
figkit.write_raw_xlsx(f'{T}/figures/raw_data/Table2_raw.xlsx', df, inputs, os.path.abspath(__file__), cols, checks,
                      notes=[f'Environment: {n.env}', 'P2 uses type 2 diabetes from Xue 2018 (GCST006867), as preregistered.'])
df.to_csv(f'{T}/figures/Table2_prereg.tsv', sep='\t', index=False)
print(df[['item', 'observed', 'ci', 'grade']].to_string(index=False))
