#!/usr/bin/env python3
"""C9 P3 and lead loci (preregistered: Stage 1 section 5, P3; Stage 1 supplement 01, items 1 and 2).

Distance clumping: GWS SNPs (P < 5e-8) sorted by P; lead = smallest; SNPs within +/-500 kb (hg19, 1000G bim)
join its locus; repeat. P3 = share of FinnGen OSA loci whose lead +/-500 kb contains any NonBMI-OSA SNP
with P < 5e-8. Also flags, for each NonBMI-OSA locus, whether it is already GWS in the BMI or OSA input
(null 3: share of leads already significant in any single input GWAS).
Usage: 04_P3_loci.py [tag=all]
"""
import sys
import pandas as pd

TAG = sys.argv[1] if len(sys.argv) > 1 else 'all'
TOPIC = '.'
STD = 'data/hm3_std'
REF = 'tools/genetics_env/ref/reference.1000G.EUR.maf.txt'
GWS, WIN = 5e-8, 500_000

pos = pd.read_csv(REF, sep='\t', usecols=['SNP', 'CHR', 'BP'])
mhc = (pos.CHR == 6) & pos.BP.between(25e6, 35e6)
pos = pos[~mhc]


def load_std(t):
    d = pd.read_csv(f'{STD}/{t}.tsv.gz', sep='\t', usecols=['SNP', 'P', 'BETA', 'SE'])
    return d.merge(pos, on='SNP')


def clump(d, pcol='P'):
    g = d[d[pcol] < GWS].sort_values(pcol).copy()
    loci = []
    while len(g):
        lead = g.iloc[0]
        same = (g.CHR == lead.CHR) & ((g.BP - lead.BP).abs() <= WIN)
        loci.append(dict(lead_SNP=lead.SNP, CHR=int(lead.CHR), BP=int(lead.BP), P=lead[pcol], n_gws_snps=int(same.sum())))
        g = g[~same]
    return pd.DataFrame(loci)


def any_gws_near(loci, other, pcol='P'):
    o = other[other[pcol] < GWS]
    hit = []
    for _, l in loci.iterrows():
        near = o[(o.CHR == l.CHR) & ((o.BP - l.BP).abs() <= WIN)]
        hit.append(len(near) > 0)
    return hit


osa = load_std('OSA_FG'); bmi = load_std('BMI_GIANT')
nb = pd.read_csv(f'{TOPIC}/derived/gwas/userGWAS_F_NonBMI_{TAG}.tsv.gz', sep='\t')
nb = nb.rename(columns={'Pval_Estimate': 'P_nonbmi', 'est': 'BETA_nonbmi', 'SE': 'SE_nonbmi'})
nb = nb[['SNP', 'P_nonbmi', 'BETA_nonbmi', 'SE_nonbmi', 'Q_SNP', 'Q_SNP_pval'] if 'Q_SNP' in nb else ['SNP', 'P_nonbmi', 'BETA_nonbmi', 'SE_nonbmi']]
nb = nb.merge(pos, on='SNP')

L_osa = clump(osa)
L_osa['retained_in_NonBMI'] = any_gws_near(L_osa, nb, 'P_nonbmi')
P3 = L_osa['retained_in_NonBMI'].mean() if len(L_osa) else float('nan')
grade = 'PASS' if P3 >= 0.30 else ('PARTIAL' if P3 >= 0.15 else 'FAIL')

L_nb = clump(nb.rename(columns={'P_nonbmi': 'P'}))
if len(L_nb):
    L_nb['GWS_in_BMI_input'] = any_gws_near(L_nb, bmi)
    L_nb['GWS_in_OSA_input'] = any_gws_near(L_nb, osa)
    L_nb['GWS_in_any_input'] = L_nb.GWS_in_BMI_input | L_nb.GWS_in_OSA_input
    L_nb = L_nb.merge(nb[['SNP', 'BETA_nonbmi', 'SE_nonbmi']], left_on='lead_SNP', right_on='SNP', how='left').drop(columns='SNP')

L_osa.to_csv(f'{TOPIC}/results/04_OSA_FG_loci_{TAG}.tsv', sep='\t', index=False)
L_nb.to_csv(f'{TOPIC}/results/04_NonBMI_loci_{TAG}.tsv', sep='\t', index=False)
summ = pd.DataFrame([
    ('FinnGen OSA loci (distance clumping, 500 kb)', len(L_osa)),
    ('... retained in NonBMI-OSA (any GWS SNP within 500 kb)', int(L_osa['retained_in_NonBMI'].sum())),
    ('P3 share retained', round(P3, 4)),
    ('P3 grade (PASS >= 0.30; PARTIAL 0.15 to 0.30)', grade),
    ('NonBMI-OSA loci', len(L_nb)),
    ('NonBMI-OSA loci already GWS in BMI input', int(L_nb['GWS_in_BMI_input'].sum()) if len(L_nb) else 0),
    ('NonBMI-OSA loci already GWS in OSA input', int(L_nb['GWS_in_OSA_input'].sum()) if len(L_nb) else 0),
    ('NonBMI-OSA loci GWS in neither input (multivariate-only)', int((~L_nb['GWS_in_any_input']).sum()) if len(L_nb) else 0),
    ('GWAS tag', TAG)], columns=['item', 'value'])
summ.to_csv(f'{TOPIC}/results/04_P3_summary_{TAG}.tsv', sep='\t', index=False)
print(summ.to_string(index=False))
