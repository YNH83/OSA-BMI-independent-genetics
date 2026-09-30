#!/usr/bin/env python3
"""Supplementary tables for the preregistered robustness analyses (C9 addendum 04), with graded read-outs.
STable14_split_half, STable15_multireg, STable16_coloc
Each table: figures/<name>.tsv + figures/raw_data/<name>_raw.xlsx. Usage (from the repository root): python3 shared/stab_robustness.py"""
import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figkit as fk

T = '.'
R, F = T + '/results', T + '/figures'
me = os.path.abspath(__file__)
N = ('STable14_split_half', 'STable15_multireg', 'STable16_coloc')

# ---- split half
s = pd.read_csv(f'{R}/16_split_half.tsv', sep='\t')
w = s.pivot_table(index='metric', columns='half', values=['estimate', 'lo', 'hi'], aggfunc='first')
rows = []
for m in dict.fromkeys(s.metric):
    rows.append(dict(metric=m, odd=w.loc[m, ('estimate', 'ODD')], odd_lo=w.loc[m, ('lo', 'ODD')], odd_hi=w.loc[m, ('hi', 'ODD')],
                     even=w.loc[m, ('estimate', 'EVEN')], even_lo=w.loc[m, ('lo', 'EVEN')], even_hi=w.loc[m, ('hi', 'EVEN')]))
S = pd.DataFrame(rows)
ov = lambda r: not (r.odd_hi < r.even_lo or r.even_hi < r.odd_lo)
S['CIs_overlap'] = S.apply(ov, axis=1)
sh = S[S.metric.str.startswith('share')].iloc[0]; p2 = S[S.metric.str.startswith('P2')].iloc[0]
S['readout'] = ''
S.loc[S.metric == sh.metric, 'readout'] = 'P1 stable in both halves' if (sh.odd >= 0.6 and sh.even >= 0.6 and ov(sh)) else 'P1 depends on chromosome half'
S.loc[S.metric == p2.metric, 'readout'] = 'P2 stable in both halves' if (p2.odd >= 0.5 and p2.even >= 0.5 and ov(p2)) else 'P2 depends on chromosome half'
S.to_csv(f'{F}/{N[0]}.tsv', sep='\t', index=False)
fk.write_raw_xlsx(f'{F}/raw_data/{N[0]}_raw.xlsx', S, [f'{R}/16_split_half.tsv'], me,
                  {'odd': 'estimate from odd chromosomes (GenomicSEM ldsc select = "ODD")', 'even': 'estimate from even chromosomes'},
                  [('rule', 'C9: P1 >= 0.60 and P2 >= 0.50 in both halves with overlapping CIs')])

# ---- multiple regression
if N[1]:
    M = pd.read_csv(f'{R}/17_multireg_A5.tsv', sep='\t')
    for sn, sub in M.groupby('snoring'):
        b = sub.set_index('quantity')
        bs, bi, bd = (b.loc[f'beta {x}'] for x in (sn, 'INSOMNIA_FG', 'DEPRESSION_FG'))
        verdict = ('upper-airway signal not explained by depression and insomnia' if (bs.estimate > 0 and bs.lo > 0 and bs.estimate > bi.estimate and bs.estimate > bd.estimate)
                   else 'shared; independent contributions not separable')
        M.loc[M.snoring == sn, 'readout'] = verdict
    M.to_csv(f'{F}/{N[1]}.tsv', sep='\t', index=False)
    fk.write_raw_xlsx(f'{F}/raw_data/{N[1]}_raw.xlsx', M, [f'{R}/17_multireg_A5.tsv'], me,
                      {'estimate': 'standardized regression coefficient of F_NonBMI on the trait, beta = Rxx^-1 rxy from one joint Genomic SEM fit',
                       'lo': '2.5th percentile of 1,000 parametric draws of S from N(S, V)', 'hi': '97.5th percentile'},
                      [('rule', 'snoring beta > 0, CI excludes 0 and larger than the insomnia and depression betas')])

# ---- coloc
C = pd.read_csv(f'{R}/18_coloc.tsv', sep='\t')
C.to_csv(f'{F}/{N[2]}.tsv', sep='\t', index=False)
fk.write_raw_xlsx(f'{F}/raw_data/{N[2]}_raw.xlsx', C, [f'{R}/18_coloc.tsv'], me,
                  {'PP.H3': 'posterior probability of distinct causal variants', 'PP.H4': 'posterior probability of a shared causal variant',
                   'call': 'PP.H4 >= 0.80 shared; PP.H3 >= 0.80 distinct; otherwise inconclusive'},
                  [('pairs', len(C)), ('shared', int((C.call == 'shared causal variant').sum())), ('distinct', int((C.call == 'distinct causal variants').sum()))])
print(S.round(3).to_string(index=False))
