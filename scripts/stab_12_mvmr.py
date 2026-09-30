#!/usr/bin/env python3
"""C9 Supplementary Table 12 (A8, addendum 03): multivariable MR (OSA and BMI jointly), univariable OSA robustness
(MR-PRESSO), CAUSE (causal versus sharing model) and MRlap (sample-overlap correction). Writes figures/STable12_mvmr.tsv
and figures/raw_data/STable12_raw.xlsx with one sheet per analysis."""
import glob, os, sys
import pandas as pd


def sci(p):
    """P in a x 10^b form for table notes."""
    if p >= 0.001: return f'{p:.3f}'
    m, e = f'{p:.1e}'.split('e'); return f'{m} × 10' + str(int(e)).translate(str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))
sys.path.insert(0, 'shared')
import figkit as fk

T = '.'
R = T + '/results'
mv = pd.read_csv(f'{R}/10_mvmr_A8.tsv', sep='\t')
un = pd.read_csv(f'{R}/10_mvmr_A8_univariable_OSA.tsv', sep='\t')
ml = pd.read_csv(f'{R}/10b_mrlap_A8.tsv', sep='\t')
cp, cq = [], []
for f in sorted(glob.glob(f'{R}/10b_cause_A8_p_*.tsv')):
    cp.append(pd.read_csv(f, sep='\t'))
for f in sorted(glob.glob(f'{R}/10b_cause_A8_posterior_*.tsv')):
    cq.append(pd.read_csv(f, sep='\t'))
cp = pd.concat(cp, ignore_index=True); cq = pd.concat(cq, ignore_index=True)
rows = []
for r in mv.itertuples():
    rows.append(dict(analysis='MVMR-IVW (OSA and BMI jointly)', exposure=r.exposure, outcome=r.outcome, role=r.role, n_instruments=r.nsnp,
                     estimate=r.b, lo=r.lo, hi=r.hi, P=r.p, note=f'conditional F: OSA {r.condF_OSA:.1f}, BMI {r.condF_BMI:.1f} (weak if < 10)'))
for o in cp.outcome:
    q = cq[cq.outcome == o]
    rows.append(dict(analysis='CAUSE (OSA total)', exposure='OSA', outcome=o, role='robustness', n_instruments=None,
                     estimate=float(q['gamma'].iloc[0]), lo=float(q['gamma'].iloc[1]), hi=float(q['gamma'].iloc[2]),   # rows: median, lower, upper
                     P=float(cp[cp.outcome == o].p_causal_vs_sharing.iloc[0]), note='P: causal model versus sharing-only model (one-sided); estimate = posterior median gamma'))
for r in ml.itertuples():
    rows.append(dict(analysis='MRlap (overlap-corrected)', exposure=r.exposure, outcome=r.outcome, role='robustness', n_instruments=None,
                     estimate=r.corrected_b, lo=r.corrected_b - 1.96 * r.corrected_se, hi=r.corrected_b + 1.96 * r.corrected_se, P=r.corrected_p,
                     note=f'observed {r.observed_b:.3f}; difference P {sci(r.diff_p)}; both traits standardized'))
for r in un.itertuples():
    rows.append(dict(analysis='Univariable OSA IVW (OSA-GWS instruments of the MVMR set)', exposure='OSA', outcome=r.outcome, role='robustness',
                     n_instruments=r.nsnp, estimate=r.ivw_b, lo=r.ivw_b - 1.96 * r.ivw_se, hi=r.ivw_b + 1.96 * r.ivw_se, P=r.ivw_p,
                     note=f'MR-PRESSO global P {r.presso_global_p:.2g}; outliers {int(r.presso_n_outliers)}'
                          + (f'; outlier-corrected {r.presso_corrected_b:.3f} (P {r.presso_corrected_p:.2g})' if pd.notna(r.presso_corrected_b) else '')))
D = pd.DataFrame(rows)
D.to_csv(f'{T}/figures/STable12_mvmr.tsv', sep='\t', index=False)
fk.write_raw_xlsx(f'{T}/figures/raw_data/STable12_raw.xlsx', {'summary': D, 'mvmr': mv, 'univariable_presso': un, 'cause_p': cp, 'cause_posterior': cq, 'mrlap': ml},
                  [f'{R}/10_mvmr_A8.tsv', f'{R}/10_mvmr_A8_univariable_OSA.tsv', f'{R}/10b_mrlap_A8.tsv'] + sorted(glob.glob(f'{R}/10b_cause_A8_*.tsv')),
                  os.path.abspath(__file__),
                  {'estimate': 'MVMR/IVW: log odds per log odds of OSA liability (binary outcomes) or SD (continuous); MRlap: standardized; CAUSE: gamma'},
                  [('MVMR OSA -> HTN (Pan-UKB) P', float(mv[(mv.outcome == 'HTN_PANUKB') & (mv.exposure == 'OSA')].p.iloc[0]))])
print(D[['analysis', 'exposure', 'outcome', 'estimate', 'P']].round(4).to_string(index=False))
