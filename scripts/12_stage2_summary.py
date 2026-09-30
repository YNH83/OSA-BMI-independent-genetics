#!/usr/bin/env python3
"""C9 Stage 2 (A1 to A5) results table generated directly from result files (no hand transcription).
Writes results/12_stage2_summary.tsv and figures/raw_data/STable_Stage2_raw.xlsx."""
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, 'shared')
import figkit as fk

T = '.'
R = T + '/results'
rows, inputs = [], []


def add(**k): rows.append(k)


for lab, f, status in [('A1 HTN vs T2D', '07_twoaxis_A1_HTN.tsv', 'preregistered primary'),
                       ('A1b CAD vs T2D', '07_twoaxis_A1b_CAD.tsv', 'preregistered secondary'),
                       ('A2 SBP vs T2D', '07_twoaxis_A2_SBP.tsv', 'preregistered (SBP GWAS adjusted for BMI; see note)'),
                       ('A2 DBP vs T2D', '07_twoaxis_A2_DBP.tsv', 'preregistered secondary (DBP adjusted for BMI)'),
                       ('Post hoc: FinnGen HTN vs FinnGen T2D', '07_twoaxis_posthoc_HTN_vs_T2DFG.tsv', 'post hoc sensitivity (same registry)'),
                       ('A7 HTN (Pan-UKB) vs T2D', '07_twoaxis_A7_HTN_PANUKB.tsv', 'addendum 03 primary: replication outside FinnGen, not BMI-adjusted'),
                       ('A7 SBP (Pan-UKB) vs T2D', '07_twoaxis_A7_SBP_PANUKB.tsv', 'addendum 03 secondary: continuous, not BMI-adjusted'),
                       ('A7 DBP (Pan-UKB) vs T2D', '07_twoaxis_A7_DBP_PANUKB.tsv', 'addendum 03 secondary: continuous, not BMI-adjusted')]:
    p = f'{R}/{f}'
    if not os.path.exists(p):
        add(item=lab, status=status, estimate='pending'); continue
    r = pd.read_csv(p, sep='\t').iloc[0]; inputs.append(p)
    add(item=lab, status=status, metab=r.metab, axis2=r.axis2,
        cor_NonBMI_metab=r.cor_NonBMI_metab, se_metab=r.se_metab, cor_NonBMI_axis2=r.cor_NonBMI_axis2, se_axis2=r.se_axis2,
        chisq=r.constrained_chisq, df=r.constrained_df, P=r.p_difference, grade=r.grade)
p = f'{R}/08_mr_A3.tsv'
if os.path.exists(p):
    m = pd.read_csv(p, sep='\t'); inputs.append(p)
    for r in m.itertuples():
        add(item=f'A3 MR {r.direction}', status='preregistered', method=r.method, nsnp=r.nsnp, b=r.b, lo=r.lo, hi=r.hi, P=r.pval,
            grade=f'order supported={r.A3_cholesky_order_supported}; reverse flag={r.A3_reverse_flag}')
    s = m.drop_duplicates('direction')
    for r in s.itertuples():
        add(item=f'A3 Steiger-filtered IVW {r.direction}', status='reported sensitivity', nsnp=r.ivw_steiger_filtered_n,
            b=r.ivw_steiger_filtered_b, lo=r.ivw_steiger_filtered_b - 1.96 * r.ivw_steiger_filtered_se,
            hi=r.ivw_steiger_filtered_b + 1.96 * r.ivw_steiger_filtered_se,
            P=float(__import__('math').erfc(abs(r.ivw_steiger_filtered_b / r.ivw_steiger_filtered_se) / 2 ** 0.5)))  # two-sided Wald P
panel = [('A4 household income', 'INCOME_HILL', 'preregistered (flag if |cor| >= 0.20 and CI excludes 0)'),
         ('A4 cognitive performance', 'COGNITION_LEE', 'preregistered secondary'),
         ('A5 snoring (GCST009760, main analysis)', 'SNORING_CAMPOS', 'preregistered panel (corrected source, addendum 02)'),
         ('A5 snoring adjusted for BMI (GCST009761)', 'SNORING_CAMPOS_ADJBMI', 'as registered; BMI-adjusted, collider-prone (addendum 02)'), ('A5 sleep duration', 'SLEEPDUR_DASHTI', 'preregistered panel'),
         ('A5 chronotype', 'CHRONOTYPE_JONES', 'preregistered panel'), ('A5 FEV1/FVC', 'FEV1FVC_SHRINE', 'preregistered panel'),
         ('A5 depression (FinnGen)', 'DEPRESSION_FG', 'preregistered panel (source changed, addendum 01)'),
         ('A5 insomnia (FinnGen)', 'INSOMNIA_FG', 'preregistered panel'), ('A5 asthma (FinnGen)', 'ASTHMA_FG', 'preregistered panel')]
for lab, tr, status in panel:
    p = f'{R}/02_numeric_{tr}.tsv'
    if not os.path.exists(p):
        add(item=lab, status=status, estimate='pending'); continue
    r = pd.read_csv(p, sep='\t').iloc[0]; inputs.append(p)
    lo, hi = r.rg_non_ext - 1.96 * r.rg_non_ext_se, r.rg_non_ext + 1.96 * r.rg_non_ext_se
    z = r.rg_non_ext / r.rg_non_ext_se
    from scipy.stats import norm
    pv = 2 * norm.sf(abs(z))
    g = ''
    if tr.startswith('INCOME') or tr.startswith('COGNITION'):
        g = 'FLAG' if abs(r.rg_non_ext) >= 0.20 and (lo > 0 or hi < 0) else 'no flag'
    else:
        g = ('Bonferroni-significant' if pv < 0.05 / 7 else 'not significant') + ('; |NonBMI| > |BMI|' if abs(r.rg_non_ext) > abs(r.rg_bmi_ext) else '; |NonBMI| <= |BMI|')
    add(item=lab, status=status, external=tr, rg_OSA=r.rg_osa_ext, cor_F_BMI=r.rg_bmi_ext, cor_F_NonBMI=r.rg_non_ext,
        lo=lo, hi=hi, P=pv, grade=g, same_registry_as_OSA=tr.endswith('_FG'))
p = f'{R}/11_lava_A6_verdicts.tsv'
if os.path.exists(p):
    inputs.append(p)
    for r in pd.read_csv(p, sep='\t').itertuples():
        add(item=f'A6 {r.rule}', status='addendum 03 (LAVA)', grade=r.result)
D = pd.DataFrame(rows)
D.to_csv(f'{R}/12_stage2_summary.tsv', sep='\t', index=False)
fk.write_raw_xlsx(f'{T}/figures/raw_data/STable_Stage2_raw.xlsx', D, inputs, os.path.abspath(__file__),
                  {'cor_NonBMI_metab': 'Genomic SEM correlation of F_NonBMI with the metabolic trait', 'P': 'chi-square difference P (A1/A2) or MR P or two-sided normal P from the correlation and its SE (A4/A5)',
                   'same_registry_as_OSA': 'True if the external trait is a FinnGen registry diagnosis like OSA (possible shared ascertainment)'},
                  [('A5 Bonferroni threshold', 0.05 / 7)],
                  notes=['A2 used Evangelou 2018 SBP and DBP, which were adjusted for BMI (verified in the Methods of PMC6284793).',
                         'Post hoc rows are not preregistered and are reported as sensitivity analyses.'])
print(D.drop(columns=[c for c in D.columns if c in ('status',)]).round(4).to_string(index=False))
