#!/usr/bin/env python3
"""MAGMA gene-based and GTEx v8 tissue gene-property analysis (FUMA-equivalent settings).

Steps
 1. SNP -> gene annotation: NCBI37.3 gene locations, window 0 kb (FUMA default), hg19 positions from the 1000G bim.
 2. Gene analysis: SNP-wise mean model with the 1000G EUR reference (g1000_eur).
    N: per-SNP N for real GWAS; for Genomic SEM latent factors the implied effective N is
    mean over SNPs with 0.1 <= MAF <= 0.4 of 1 / (2 p q SE^2)  (Mallard et al. 2022 convention).
 3. Tissue enrichment: gene-property analysis of GTEx v8 median expression log2(TPM + 1) per tissue,
    conditioning on the average across tissues, one-sided (direction=greater), as in FUMA.
    GTEx Ensembl IDs are mapped to NCBI Entrez via gene symbol; the mapping rate is reported.

Usage: magma_pipeline.py OUTDIR LABEL GWAS_FILE SNPCOL PCOL [NCOL|--latent BETA SE MAF]
"""
import os, subprocess, sys
import numpy as np
import pandas as pd

M = 'tools/genetics_env/magma'
BIN = os.path.join(M, 'magma')
GENELOC = os.path.join(M, 'NCBI37.3', 'NCBI37.3.gene.loc')
BFILE = os.path.join(M, 'g1000_eur', 'g1000_eur')
REF = 'tools/genetics_env/ref/reference.1000G.EUR.maf.txt'
GTEX = os.path.join(M, 'GTEx_v8_gene_median_tpm.gct.gz')
COVAR = os.path.join(M, 'gtex_v8_log2tpm_entrez.txt')


def run(cmd, log):
    with open(log, 'a') as f:
        f.write('$ ' + ' '.join(cmd) + '\n')
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if p.returncode:
        raise SystemExit(f'MAGMA failed ({p.returncode}); see {log}')


def build_covar():
    if os.path.exists(COVAR):
        return
    g = pd.read_csv(GTEX, sep='\t', skiprows=2)
    loc = pd.read_csv(GENELOC, sep='\t', header=None, names=['entrez', 'chr', 'start', 'end', 'strand', 'symbol'])
    expr = g.drop(columns=['Name']).groupby('Description').mean()
    expr = np.log2(expr + 1)
    m = loc[['entrez', 'symbol']].merge(expr, left_on='symbol', right_index=True, how='inner')
    m = m.drop(columns='symbol').drop_duplicates('entrez')
    m.columns = ['GENE'] + [c.replace(' ', '_').replace('(', '').replace(')', '').replace('-', '_') for c in m.columns[1:]]
    m['Average'] = m.iloc[:, 1:].mean(axis=1)
    m.to_csv(COVAR, sep='\t', index=False)
    with open(COVAR + '.mapping.txt', 'w') as f:
        f.write(f'GTEx genes with symbol: {expr.shape[0]}\nNCBI37.3 genes: {len(loc)}\nmapped by symbol: {len(m)}\n'
                f'mapping rate among NCBI genes: {len(m) / len(loc):.3f}\n')


def main(outdir, label, gwas, snpcol, pcol, rest):
    os.makedirs(outdir, exist_ok=True)
    log = os.path.join(outdir, f'{label}.magma.log')
    d = pd.read_csv(gwas, sep='\t', low_memory=False)
    d = d.rename(columns={snpcol: 'SNP', pcol: 'P'})
    if rest[0] == '--latent':
        b, se, maf = rest[1:4]
        sub = d[(d[maf] >= 0.1) & (d[maf] <= 0.4)]
        neff = float(np.mean(1 / (2 * sub[maf] * (1 - sub[maf]) * sub[se] ** 2)))
        d['N'] = round(neff)
        nnote = f'latent-factor implied Neff = {neff:.0f} from {len(sub)} SNPs with MAF 0.1 to 0.4'
    else:
        d = d.rename(columns={rest[0]: 'N'}); nnote = f'per-SNP N from column {rest[0]}'
    pos = pd.read_csv(REF, sep='\t', usecols=['SNP', 'CHR', 'BP'])
    d = d[['SNP', 'P', 'N']].dropna().merge(pos, on='SNP')
    d = d[(d.P > 0) & (d.P <= 1)]
    nfloor = int((d.P < 1e-300).sum())
    d['P'] = d.P.clip(lower=1e-300)   # MAGMA cannot parse subnormal doubles (e.g. 6.7e-315); floor at 1e-300
    loc = os.path.join(outdir, f'{label}.snploc'); pv = os.path.join(outdir, f'{label}.pval')
    d[['SNP', 'CHR', 'BP']].to_csv(loc, sep='\t', index=False, header=False)
    d[['SNP', 'P', 'N']].to_csv(pv, sep='\t', index=False)
    pre = os.path.join(outdir, label)
    run([BIN, '--annotate', '--snp-loc', loc, '--gene-loc', GENELOC, '--out', pre], log)
    run([BIN, '--bfile', BFILE, '--pval', pv, 'ncol=N', '--gene-annot', pre + '.genes.annot', '--out', pre], log)
    build_covar()
    run([BIN, '--gene-results', pre + '.genes.raw', '--gene-covar', COVAR, '--model', 'condition-hide=Average',
         'direction=greater', '--out', pre + '.gtex'], log)
    with open(os.path.join(outdir, f'{label}.N_note.txt'), 'w') as f:
        f.write(nnote + f'\nSNPs passed to MAGMA: {len(d)}\nP values floored at 1e-300: {nfloor}\n')
    print(label, nnote, 'SNPs', len(d))


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], a[2], a[3], a[4], a[5:])
