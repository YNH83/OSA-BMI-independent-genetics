#!/usr/bin/env python3
"""C9 Supplementary Table 7: null control. The identical subtraction pipeline with gout (FinnGen R12 GOUT_STRICT) in place of OSA,
side by side with OSA: BMI-independent share (P1 analogue), multivariable GWAS (N_hat, GWS SNPs), MAGMA tissues and genes.
Writes figures/STable7_null.tsv and figures/raw_data/STable7_raw.xlsx."""
import os, sys
import pandas as pd
sys.path.insert(0, 'shared')
import figkit as fk

T = '.'
R = T + '/results'
ins = [f'{R}/02_numeric_T2D_XUE.tsv', f'{R}/02_numeric_T2D_XUE__GOUT_FG__BMI_GIANT.tsv', f'{R}/09_postprocess_summary_all_fixP.tsv',
       f'{R}/09_postprocess_summary_null_gout_all.tsv', f'{R}/magma/F_NonBMI.gtex.gsa.out', f'{R}/magma/GOUT_NonBMI.gtex.gsa.out',
       f'{R}/magma/F_NonBMI.genes.out', f'{R}/magma/GOUT_NonBMI.genes.out',
       f'{T}/derived/gwas/F_NonBMI_QC.tsv.gz', f'{T}/derived/gwas/null_gout_F_NonBMI_QC.tsv.gz']
osa = pd.read_csv(ins[0], sep='\t').iloc[0]; gout = pd.read_csv(ins[1], sep='\t').iloc[0]
po = pd.read_csv(ins[2], sep='\t').query("factor == 'F_NonBMI'").set_index('item').value
pg = pd.read_csv(ins[3], sep='\t').query("factor == 'F_NonBMI'").set_index('item').value
def mag(t, g):
    d = pd.read_csv(t, sep=r'\s+', comment='#'); G = pd.read_csv(g, sep=r'\s+')
    return int((d.P < 0.05 / len(d)).sum()), d.sort_values('P').VARIABLE.iloc[0], float(d.P.min()), int((G.P < 0.05 / len(G)).sum())
def lam(f):  # lambda GC at full precision from the post-QC file (the 4-decimal value in the 09 summary rounds 1.23454 to 1.2345, then to 1.234)
    d = pd.read_csv(f, sep='\t', usecols=['BETA', 'SE'])
    return float(((d.BETA / d.SE) ** 2).median() / 0.4549364)
lo, lg = lam(ins[8]), lam(ins[9])
mo, mg = mag(ins[4], ins[6]), mag(ins[5], ins[7])
rows = [('Trait subtracted from', 'OSA (FinnGen R12)', 'Gout (FinnGen R12 GOUT_STRICT)'),
        ('rg with BMI', f'{osa.s_bmi:.3f}', f'{gout.s_bmi:.3f}'),
        ('BMI-independent share of SNP heritability', f'{osa.P1:.3f} ({osa.P1 - 1.96 * osa.P1_se:.3f} to {osa.P1 + 1.96 * osa.P1_se:.3f})',
         f'{gout.P1:.3f} ({gout.P1 - 1.96 * gout.P1_se:.3f} to {gout.P1 + 1.96 * gout.P1_se:.3f})'),
        ('Residual component: N_hat', f'{float(po["N_hat (MAF 0.1 to 0.4)"]):,.0f}', f'{float(pg["N_hat (MAF 0.1 to 0.4)"]):,.0f}'),
        ('Residual component: GWS SNPs', f'{int(po["GWS SNPs P < 5e-8"])}', f'{int(pg["GWS SNPs P < 5e-8"])}'),
        ('Residual component: lambda GC', f'{lo:.3f}', f'{lg:.3f}'),
        ('MAGMA GTEx tissues below Bonferroni', f'{mo[0]} (top: {mo[1]}, P {mo[2]:.2g})', f'{mg[0]} (top: {mg[1]}, P {mg[2]:.2g})'),
        ('MAGMA genes below Bonferroni', f'{mo[3]}', f'{mg[3]}')]
D = pd.DataFrame(rows, columns=['item', 'OSA minus BMI (F_NonBMI)', 'Gout minus BMI (null control)'])
D.to_csv(f'{T}/figures/STable7_null.tsv', sep='\t', index=False)
fk.write_raw_xlsx(f'{T}/figures/raw_data/STable7_raw.xlsx', D, ins, os.path.abspath(__file__),
                  {'Gout minus BMI (null control)': 'identical Cholesky subtraction, LDSC settings, userGWAS model and MAGMA pipeline, gout in place of OSA'},
                  [('OSA P1', float(osa.P1)), ('gout P1 analogue', float(gout.P1))])
print(D.to_string(index=False))
