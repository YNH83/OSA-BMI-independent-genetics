#!/usr/bin/env python3
"""Standardize raw GWAS summary statistics to one format restricted to HapMap3 SNPs.

Output: data/hm3_std/<TRAIT>.tsv.gz with columns SNP A1 A2 BETA SE P N MAF (A1 = effect allele).
Also writes data/hm3_std/<TRAIT>.meta.json (source, counts, Neff, prevalences) and appends to data/standardize_log.tsv.
Binary traits: N = effective sample size 4 / (1/Ncase + 1/Ncontrol), constant across SNPs (per-SNP N not provided).
Raw files are read only, never modified.

Usage (from the repository root): python3 data_prep/standardize.py TRAIT [TRAIT ...]   (or ALL)
Raw files are expected under data/raw/ in the folder layout of data_prep/download_manifest.tsv.
"""
import gzip, json, os, sys, time
import numpy as np
import pandas as pd

RAW = 'data/raw'
OUT = os.environ.get('STD_OUT', 'data/hm3_std')
HM3 = 'ref/ldsc_ref/w_hm3.snplist'
REF = 'tools/genetics_env/ref/reference.1000G.EUR.maf.txt'   # 1000G Phase 3 EUR SNP CHR BP MAF A1 A2
FG = RAW + '/FinnGen R12/finngen_R12_{}.gz'
GC = RAW + '/NHGRI-EBI GWAS Catalog/{}/{}'


def neff(ca, co):
    return 4.0 / (1.0 / ca + 1.0 / co)


FINNGEN_COLS = dict(SNP='rsids', A1='alt', A2='ref', BETA='beta', SE='sebeta', P='pval', FRQ='af_alt')
GCH_COLS = dict(SNP='rsid', A1='effect_allele', A2='other_allele', BETA='beta', SE='standard_error', P='p_value', FRQ='effect_allele_frequency')
GCH_OLD = dict(SNP='hm_rsid', A1='hm_effect_allele', A2='hm_other_allele', BETA='hm_beta', SE='standard_error', P='p_value', FRQ='hm_effect_allele_frequency')

TRAITS = {
 # ---- FinnGen R12 (hg38; beta for alt allele)
 'OSA_FG':          dict(path=FG.format('G6_SLEEPAPNO'), cols=FINNGEN_COLS, case=56885, ctrl=441137),
 'OSA_FG_INCLAVO':  dict(path=FG.format('G6_SLEEPAPNO_INCLAVO'), cols=FINNGEN_COLS, case=62722, ctrl=437626),
 'HTN_FG':          dict(path=FG.format('I9_HYPTENSESS'), cols=FINNGEN_COLS, case=132515, ctrl=345634),
 'GOUT_FG':         dict(path=FG.format('GOUT_STRICT'), cols=FINNGEN_COLS, case=6715, ctrl=490619),
 'T2D_FG':          dict(path=FG.format('T2D'), cols=FINNGEN_COLS, case=82878, ctrl=403489),
 'INSOMNIA_FG':    dict(path=FG.format('F5_INSOMNIA'), cols=FINNGEN_COLS, case=6776, ctrl=490763),
 'ASTHMA_FG':      dict(path=FG.format('J10_ASTHMA_EXMORE'), cols=FINNGEN_COLS, case=59100, ctrl=259839),
 'DEPRESSION_FG':  dict(path=FG.format('F5_DEPRESSIO'), cols=FINNGEN_COLS, case=59333, ctrl=434831),
 # ---- GWAS Catalog
 'BMI_GIANT':       dict(path=GC.format('GCST006900', 'Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz'),
                         cols=dict(SNP='SNP', A1='Tested_Allele', A2='Other_Allele', BETA='BETA', SE='SE', P='P', FRQ='Freq_Tested_Allele_in_HRS', N='N')),
 'OSA_CAMPOS':      dict(path=GC.format('GCST90267361', 'GCST90267361.h.tsv.gz'), cols=GCH_COLS, case=25008, ctrl=337630),
 'OSA_CAMPOS_ADJBMI': dict(path=GC.format('GCST90267363', 'GCST90267363.h.tsv.gz'), cols=GCH_COLS, case=25008, ctrl=337630),
 'OSA_STRAUSZ_BMILT30': dict(path=GC.format('GCST90624797', 'GCST90624797.h.tsv.gz'), cols=GCH_COLS, case=20413, ctrl=443463),
 'T2D_XUE':         dict(path=GC.format('GCST006867', '30054458-GCST006867-EFO_0001360.h.tsv.gz'), cols=GCH_OLD, case=61714, ctrl=593952),
 'CAD_ARAGAM':      dict(path=GC.format('GCST90132314', 'GCST90132314.h.tsv.gz'), cols=GCH_COLS, case=181522, ctrl=984168),
 # Evangelou 2018 blood pressure (METAL output, whitespace-separated): MarkerName 'chr:pos:SNP' (GRCh37) mapped to
 # rsID via 1000G bim positions; Allele1 = effect allele (METAL convention); per-SNP N_effective used as N.
 'SBP_EVANGELOU': dict(path=GC.format('GCST006624', 'Evangelou_30224653_SBP.txt.gz'), sep=r'\s+', marker_chrpos=True,
                       cols=dict(SNP='MarkerName', A1='Allele1', A2='Allele2', BETA='Effect', SE='StdErr', P='P', FRQ='Freq1', N='N_effective')),
 'DBP_EVANGELOU': dict(path=GC.format('GCST006630', 'Evangelou_30224653_DBP.txt.gz'), sep=r'\s+', marker_chrpos=True,
                       cols=dict(SNP='MarkerName', A1='Allele1', A2='Allele2', BETA='Effect', SE='StdErr', P='P', FRQ='Freq1', N='N_effective')),
 # Hill 2019 household income (non-MTAG; readme: Effect_Allele = effect allele, GRCh37, N = 286,301, no allele frequency)
 'INCOME_HILL': dict(path=GC.format('GCST009523', 'HillWD_31844048_household_Income.txt.gz'), sep=r'\s+', n=286301,
                     cols=dict(SNP='SNP', A1='Effect_Allele', A2='Non_effect_Allele', BETA='Beta', SE='Standard_Error_of_Beta', P='P')),
 # Campos 2020 snoring (UKB, BOLT-LMM on a binary trait; old-style harmonised files); LDSC uses Z, Neff from case/control counts.
 # GCST009760 = main analysis (NOT BMI-adjusted); GCST009761 = BMI-adjusted sensitivity GWAS (raw file name ...Snoring-adjustedbyBMI).
 'SNORING_CAMPOS': dict(path=GC.format('GCST009760', '32060260-GCST009760-EFO_0008341.h.tsv.gz'), cols=GCH_OLD, case=152302, ctrl=256015),
 'SNORING_CAMPOS_ADJBMI': dict(path=GC.format('GCST009761', '32060260-GCST009761-EFO_0008341.h.tsv.gz'), cols=GCH_OLD, case=151836, ctrl=255230),
 # Shrine 2019 FEV1/FVC (old-style harmonised; per-SNP n; adjusted for age, age2, sex, height, smoking; not BMI)
 'FEV1FVC_SHRINE': dict(path=GC.format('GCST007431', '30804560-GCST007431-EFO_0004713.h.tsv.gz'), cols=dict(GCH_OLD, N='n'),
                        hm_fallback=dict(SNP='variant_id', A1='effect_allele', A2='other_allele', BETA='beta', FRQ='effect_allele_frequency')),   # 13,434 hm_code 15 rows (harmonisation failed) kept from the original columns; audit A W2, 2026-09-27
 # Pulit 2019 WHRadjBMI (raw GIANT/UKB format, md5 verified against the publisher's md5sum.txt): SNP column is
 # 'rsID:ref:alt', Allele1 is the effect allele (lower case), no per-SNP N -> constant N = 694,649 (GWAS Catalog).
 'WHRADJBMI_PULIT': dict(path=GC.format('GCST008999', 'fat-distn.giant.ukbb.meta-analysis.whradjbmi.combined.tbl.gz'),
                         cols=dict(SNP='SNP', A1='Allele1', A2='Allele2', BETA='Effect', SE='StdErr', P='P-value', FRQ='Freq1'),
                         n=694649, snp_sep=':'),
 # Dashti 2019 sleep duration (UKB, BOLT-LMM, raw zip with one txt): ALLELE1 = effect allele, A1FREQ = ALLELE1 frequency;
 # no per-SNP N -> constant N = 446,118 (PMID 30846698 abstract, checked 2026-09-27); INFO kept out (LDSC applies its own filter).
 'SLEEPDUR_DASHTI': dict(path=GC.format('GCST007561', 'Dashti_30846698_sleepdurationsumstats.txt.zip'), compression='zip', n=446118,
                         cols=dict(SNP='SNP', A1='ALLELE1', A2='ALLELE0', BETA='BETA_SLEEPDURATION', SE='SE_SLEEPDURATION',
                                   P='P_SLEEPDURATION', FRQ='A1FREQ')),
 # Jones 2019 chronotype (UKB-only BOLT-LMM raw file; GWAS Catalog GCST007576 initialSampleSize 449,734 European, checked
 # via the Catalog REST API 2026-09-27): ALLELE1 = effect allele, A1FREQ = ALLELE1 frequency; treated as a quantitative
 # trait (only rg is used downstream, so the observed-scale h2 is not interpreted).
 'CHRONOTYPE_JONES': dict(path=GC.format('GCST007576', 'chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.gz'),
                          n=449734, cols=dict(SNP='SNP', A1='ALLELE1', A2='ALLELE0', BETA='BETA', SE='SE', P='P_BOLT_LMM', FRQ='A1FREQ')),
 # Lee 2018 cognitive performance (COGENT + UKB, old-style harmonised; hm_beta is flipped consistently with hm_effect_allele,
 # checked on rs3740422); no per-SNP N -> constant N = 257,841 (GWAS Catalog REST API initialSampleSize, checked 2026-09-27).
 'COGNITION_LEE': dict(path=GC.format('GCST006572', '30038396-GCST006572-EFO_0008354.h.tsv.gz'), cols=GCH_OLD, n=257841),
 # Pan-UKB (EUR columns, GRCh37; SAIGE; covariates age, sex, age*sex, age2, age2*sex, 10 PCs, no BMI, per pan.ukbb.broadinstitute.org/docs/qc).
 # alt = effect allele; P from neglog10_pval_EUR; variants flagged low_confidence_EUR dropped; rsID by hg19 chr:pos (1000G bim);
 # EUR N from the Pan-UKB phenotype manifest (2026-09-27): SBP 396,663; DBP 396,667; hypertension phecode 401 94,311 / 325,488.
 'SBP_PANUKB': dict(path=f'{RAW}/Pan-UK Biobank/continuous-4080-both_sexes-irnt.tsv.bgz', n=396663, p_neglog10=True,
                    cols=dict(CHR='chr', BP='pos', A2='ref', A1='alt', FRQ='af_EUR', BETA='beta_EUR', SE='se_EUR', P='neglog10_pval_EUR', LOWCONF='low_confidence_EUR')),
 'DBP_PANUKB': dict(path=f'{RAW}/Pan-UK Biobank/continuous-4079-both_sexes-irnt.tsv.bgz', n=396667, p_neglog10=True,
                    cols=dict(CHR='chr', BP='pos', A2='ref', A1='alt', FRQ='af_EUR', BETA='beta_EUR', SE='se_EUR', P='neglog10_pval_EUR', LOWCONF='low_confidence_EUR')),
 'HTN_PANUKB': dict(path=f'{RAW}/Pan-UK Biobank/phecode-401-both_sexes.tsv.bgz', case=94311, ctrl=325488, p_neglog10=True,
                    # binary phenotypes have no af_EUR column; allele frequency in EUR controls is used (MAF filter only)
                    cols=dict(CHR='chr', BP='pos', A2='ref', A1='alt', FRQ='af_controls_EUR', BETA='beta_EUR', SE='se_EUR', P='neglog10_pval_EUR', LOWCONF='low_confidence_EUR')),
}


_posmap = None
def POSMAP():
    global _posmap
    if _posmap is None:
        r = pd.read_csv(REF, sep='\t', usecols=['SNP', 'CHR', 'BP'])
        r = r[r.SNP.isin(set(pd.read_csv(HM3, sep='\t').SNP))]
        _posmap = dict(zip('chr' + r.CHR.astype(str) + ':' + r.BP.astype(str), r.SNP))
    return _posmap


def load_hm3():
    s = pd.read_csv(HM3, sep='\t')
    return s.set_index('SNP')[['A1', 'A2']]


def run(trait, hm3):
    cfg = TRAITS[trait]
    c = cfg['cols']
    t0 = time.time()
    use = [v for v in c.values()]
    fb = cfg.get('hm_fallback', {})
    use += [v for v in fb.values() if v not in use]
    if 'CHR' in c and 'SNP' not in c:   # position-keyed files (Pan-UKB): SNP is built from chr:pos below
        pass
    SEP = cfg.get('sep', '\t')
    head = pd.read_csv(cfg['path'], sep=SEP, nrows=5, compression=cfg.get('compression', 'gzip'))
    missing = [v for v in use if v not in head.columns]
    if missing:
        raise SystemExit(f'{trait}: columns not found {missing}; available {list(head.columns)}')
    parts, n_raw = [], 0
    for ch in pd.read_csv(cfg['path'], sep=SEP, usecols=use, compression=cfg.get('compression', 'gzip'), chunksize=2_000_000,
                          dtype={c['SNP']: str} if 'SNP' in c else None, low_memory=False):
        n_raw += len(ch)
        ch = ch.rename(columns={v: k for k, v in c.items()})
        for k, v in fb.items():   # harmonised field missing (e.g. hm_code 15): use the author-reported field of the same row
            miss = ch[k].isna() | (ch[k].astype(str).isin(['nan', 'NA', '']))
            ch.loc[miss, k] = ch.loc[miss, v]
        if 'LOWCONF' in ch:
            ch = ch[ch['LOWCONF'].astype(str).str.lower() != 'true']
        if 'CHR' in ch and 'SNP' not in ch:
            ch['SNP'] = ('chr' + ch['CHR'].astype(str) + ':' + ch['BP'].astype(str)).map(POSMAP())
            ch = ch.dropna(subset=['SNP'])
        if cfg.get('p_neglog10'):
            ch['P'] = 10 ** (-pd.to_numeric(ch['P'], errors='coerce'))
        ch['SNP'] = ch['SNP'].astype(str).str.split(cfg.get('snp_sep', ',')).str[0]   # FinnGen: 'rs1,rs2'; Pulit: 'rs:ref:alt'
        if cfg.get('marker_chrpos'):
            mk = ch['SNP'].str.split(':')   # do not reuse the name 'parts' (the chunk list)
            ch['SNP'] = ('chr' + mk.str[0] + ':' + mk.str[1]).map(POSMAP())
            ch = ch.dropna(subset=['SNP'])
        ch = ch[ch['SNP'].isin(hm3.index)]
        parts.append(ch)
    d = pd.concat(parts, ignore_index=True)
    d['A1'] = d['A1'].astype(str).str.upper(); d['A2'] = d['A2'].astype(str).str.upper()
    ref = hm3.loc[d['SNP']]
    same = (d['A1'].values == ref['A1'].values) & (d['A2'].values == ref['A2'].values)
    flip = (d['A1'].values == ref['A2'].values) & (d['A2'].values == ref['A1'].values)
    d = d[same | flip].copy()
    if 'FRQ' not in d:
        d['FRQ'] = float('nan')   # e.g. Hill 2019 income has no allele frequency column
    for k in ('BETA', 'SE', 'P', 'FRQ'):
        d[k] = pd.to_numeric(d[k], errors='coerce')
    d = d.dropna(subset=['BETA', 'SE', 'P'])
    # P reported as 0 is numerical underflow at the strongest loci (e.g. FTO in Yengo 2018 BMI: 0.0e+00). Until 2026-09-27 these
    # rows were dropped by 'P > 0', which silently removed the top signals (deviation record 03). Recover P from BETA/SE, floored at 1e-300.
    from scipy.stats import norm as _norm
    under = (d['P'] <= 0) & (d['SE'] > 0)
    n_under = int(under.sum())
    if n_under:
        z = (d.loc[under, 'BETA'] / d.loc[under, 'SE']).abs()
        d.loc[under, 'P'] = np.maximum(2 * _norm.sf(z), 1e-300)
    d = d[(d['SE'] > 0) & (d['P'] > 0) & (d['P'] <= 1)]
    d = d.drop_duplicates('SNP', keep=False)
    d['MAF'] = d['FRQ'].where(d['FRQ'] <= 0.5, 1 - d['FRQ'])
    if d['MAF'].isna().all():   # no frequency in the source file: use 1000G Phase3 EUR MAF (only used for the MAF >= 0.01 filter)
        ref = pd.read_csv(REF, sep='\t', usecols=['SNP', 'MAF'])
        d = d.drop(columns='MAF').merge(ref, on='SNP', how='left')
    if 'N' in c:
        d['N'] = pd.to_numeric(d['N'], errors='coerce')
        nmode = 'per-SNP N from file'
    elif 'n' in cfg:
        d['N'] = cfg['n']; nmode = f"constant N {cfg['n']}"
    else:
        d['N'] = round(neff(cfg['case'], cfg['ctrl']), 1)
        nmode = f"Neff {d['N'].iloc[0]} from {cfg['case']} cases / {cfg['ctrl']} controls"
    d = d[['SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'N', 'MAF']]
    os.makedirs(OUT, exist_ok=True)
    d.to_csv(f'{OUT}/{trait}.tsv.gz', sep='\t', index=False, compression='gzip', float_format='%.6g')
    meta = dict(trait=trait, source=cfg['path'], n_raw_rows=n_raw, n_hm3_kept=len(d),
                n_mode=nmode, binary='case' in cfg,
                sample_prev=(cfg['case'] / (cfg['case'] + cfg['ctrl'])) if 'case' in cfg else None,
                min_p=float(d['P'].min()), n_gws=int((d['P'] < 5e-8).sum()), n_p_underflow_recovered=n_under,
                seconds=round(time.time() - t0), created=time.strftime('%Y-%m-%d %H:%M:%S'))
    json.dump(meta, open(f'{OUT}/{trait}.meta.json', 'w'), indent=1, ensure_ascii=False)
    with open(OUT + '/../standardize_log.tsv', 'a') as f:
        f.write('\t'.join(str(meta[k]) for k in ('created', 'trait', 'n_raw_rows', 'n_hm3_kept', 'n_mode', 'n_gws')) + '\n')
    print(trait, 'raw', n_raw, 'hm3 kept', len(d), nmode, 'GWS', meta['n_gws'], f"{meta['seconds']}s", flush=True)


if __name__ == '__main__':
    want = list(TRAITS) if sys.argv[1:] == ['ALL'] else sys.argv[1:]
    hm3 = load_hm3()
    for t in want:
        if not os.path.exists(TRAITS[t]['path']):
            print(t, 'RAW FILE NOT YET AVAILABLE', flush=True); continue
        run(t, hm3)
