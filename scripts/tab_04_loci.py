#!/usr/bin/env python3
"""C9 Table 4 (lead loci of F_NonBMI) and Supplementary Table 4 (all genome-wide significant F_NonBMI SNPs), from the
post-fix GWAS (tag all_fixP, deviation record 03). Lead loci: 500 kb distance clumping (Stage 1 supplement 01). For each lead: nearest gene
(MAGMA NCBI37.3, LOC/LINC/MIR/orf names skipped), effects on both factors, P in the OSA and BMI input GWAS, and whether the
locus is a retained OSA locus or new. Writes figures/Table4_loci.tsv, figures/STable4_gws.tsv and their raw-data Excel files."""
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, 'shared')
import figkit as fk

T = '.'
STD = 'data/hm3_std'
GL = 'tools/genetics_env/magma/NCBI37.3/NCBI37.3.gene.loc'
TAG = 'all_fixP'
nb = pd.read_csv(f'{T}/derived/gwas/userGWAS_F_NonBMI_{TAG}.tsv.gz', sep='\t', usecols=['SNP', 'CHR', 'BP', 'A1', 'A2', 'MAF', 'est', 'SE', 'Pval_Estimate'])
fb = pd.read_csv(f'{T}/derived/gwas/userGWAS_F_BMI_{TAG}.tsv.gz', sep='\t', usecols=['SNP', 'est', 'SE', 'Pval_Estimate'])
fb.columns = ['SNP', 'beta_F_BMI', 'se_F_BMI', 'P_F_BMI']
nb = nb.rename(columns={'est': 'beta_F_NonBMI', 'SE': 'se_F_NonBMI', 'Pval_Estimate': 'P_F_NonBMI'}).merge(fb, on='SNP')
osa = pd.read_csv(f'{STD}/OSA_FG.tsv.gz', sep='\t', usecols=['SNP', 'P']).rename(columns={'P': 'P_OSA_input'})
bmi = pd.read_csv(f'{STD}/BMI_GIANT.tsv.gz', sep='\t', usecols=['SNP', 'P']).rename(columns={'P': 'P_BMI_input'})
genes = pd.read_csv(GL, sep=r'\s+', header=None, names=['id', 'CHR', 'start', 'end', 'strand', 'gene'], dtype={'CHR': str})
genes = genes[~genes.gene.str.match(r'^(LOC|LINC|MIR|SNOR|FLJ|KIAA)|orf|-AS\d|^RNU')]
genes = genes[genes.CHR.str.isdigit()].assign(CHR=lambda g: g.CHR.astype(int))


def nearest(c, bp):
    g = genes[genes.CHR == c]
    dist = np.where((g.start <= bp) & (g.end >= bp), 0, np.minimum((g.start - bp).abs(), (g.end - bp).abs()))
    i = int(np.argmin(dist)); return g.gene.iloc[i], int(dist[i])


gws = nb[nb.P_F_NonBMI < 5e-8].merge(osa, on='SNP', how='left').merge(bmi, on='SNP', how='left').sort_values(['CHR', 'BP'])
L = pd.read_csv(f'{T}/results/04_NonBMI_loci_{TAG}.tsv', sep='\t')
Losa = pd.read_csv(f'{T}/results/04_OSA_FG_loci_{TAG}.tsv', sep='\t')
rows = []
for r in L.sort_values(['CHR', 'BP']).itertuples():
    s = gws[gws.SNP == r.lead_SNP].iloc[0]
    g, dist = nearest(r.CHR, r.BP)
    near_osa = Losa[(Losa.CHR == r.CHR) & ((Losa.BP - r.BP).abs() <= 500_000)]
    rows.append(dict(lead_SNP=r.lead_SNP, CHR=r.CHR, BP_hg19=r.BP, nearest_gene=g, distance_bp=dist, A1=s.A1, A2=s.A2, MAF=s.MAF,
                     beta_F_NonBMI=s.beta_F_NonBMI, se_F_NonBMI=s.se_F_NonBMI, P_F_NonBMI=s.P_F_NonBMI,
                     beta_F_BMI=s.beta_F_BMI, se_F_BMI=s.se_F_BMI, P_F_BMI=s.P_F_BMI, P_OSA_input=s.P_OSA_input, P_BMI_input=s.P_BMI_input,
                     n_gws_snps=r.n_gws_snps, status='retained OSA locus' if len(near_osa) else 'new',
                     GWS_in_BMI_input=bool(r.GWS_in_BMI_input)))
D = pd.DataFrame(rows)
D.to_csv(f'{T}/figures/Table4_loci.tsv', sep='\t', index=False)
fk.write_raw_xlsx(f'{T}/figures/raw_data/Table4_raw.xlsx', D, [f'{T}/results/04_NonBMI_loci_{TAG}.tsv', f'{T}/derived/gwas/userGWAS_F_NonBMI_{TAG}.tsv.gz', GL],
                  os.path.abspath(__file__),
                  {'beta_F_NonBMI': 'effect of the A1 allele on the standardized BMI-independent factor', 'beta_F_BMI': 'effect on the BMI factor, same model',
                   'status': 'retained = within 500 kb of a FinnGen OSA lead (P3 rule); new = no OSA GWS lead within 500 kb',
                   'GWS_in_BMI_input': 'any BMI (GIANT) SNP with P < 5e-8 within 500 kb'},
                  [('F_NonBMI lead loci', len(D)), ('retained OSA loci', int((D.status == 'retained OSA locus').sum()))])
gws.to_csv(f'{T}/figures/STable4_gws.tsv', sep='\t', index=False)
fk.write_raw_xlsx(f'{T}/figures/raw_data/STable4_raw.xlsx', gws, [f'{T}/derived/gwas/userGWAS_F_NonBMI_{TAG}.tsv.gz'], os.path.abspath(__file__),
                  {'P_F_NonBMI': 'Genomic SEM F_NonBMI ~ SNP P'}, [('GWS SNPs', len(gws))])
print(D[['lead_SNP', 'CHR', 'BP_hg19', 'nearest_gene', 'beta_F_NonBMI', 'P_F_NonBMI', 'beta_F_BMI', 'P_F_BMI', 'status', 'GWS_in_BMI_input']].to_string(index=False))
print('GWS SNPs', len(gws))
