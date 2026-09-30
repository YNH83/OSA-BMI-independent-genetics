#!/usr/bin/env python3
"""C9 Table 1: input GWAS summary statistics, with raw-data Excel (figures/raw_data/Table1_raw.xlsx)."""
import json, os, sys
import pandas as pd
sys.path.insert(0, 'shared')
import figkit

TOPIC = '.'
STD = 'data/hm3_std'
import importlib.util
_sp = importlib.util.spec_from_file_location('standardize', 'data_prep/standardize.py')
STDCFG = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(STDCFG)
# case/control counts and constant N are read from the standardization config (the numbers actually used), and must agree
# with the values typed below; a mismatch stops the script.
ROLE = [
 ('OSA_FG', 'Obstructive sleep apnoea (primary)', 'FinnGen R12 G6_SLEEPAPNO', 'Finnish', 'hg38', 56885, 441137, None),
 ('BMI_GIANT', 'Body mass index (primary)', 'GIANT, Yengo 2018; GCST006900', 'European', 'hg19', None, None, 'per-SNP'),
 ('T2D_XUE', 'Type 2 diabetes (external correlate)', 'Xue 2018; GCST006867', 'European', 'hg38 harmonised', 61714, 593952, None),
 ('HTN_FG', 'Essential hypertension (external correlate)', 'FinnGen R12 I9_HYPTENSESS', 'Finnish', 'hg38', 132515, 345634, None),
 ('CAD_ARAGAM', 'Coronary artery disease (external correlate)', 'Aragam 2022; GCST90132314', 'European', 'hg38 harmonised', 181522, 984168, None),
 ('GOUT_FG', 'Gout (null control)', 'FinnGen R12 GOUT_STRICT', 'Finnish', 'hg38', 6715, 490619, None),
 ('OSA_CAMPOS_ADJBMI', 'Sleep apnoea, adjusted for BMI (triangulation; mtCOJO per the publication, Catalog label "adjusted for BMI")', 'Campos 2023; GCST90267363', 'European', 'hg38 harmonised', 25008, 337630, None),
 ('OSA_STRAUSZ_BMILT30', 'Sleep apnoea in BMI < 30 (triangulation)', 'Strausz 2025; GCST90624797', 'Finnish', 'hg38 harmonised', 20413, 443463, None),
 ('OSA_CAMPOS', 'Sleep apnoea, consistency check (overlaps FinnGen: 9,096 of 25,008 cases)', 'Campos 2023; GCST90267361', 'European', 'hg38 harmonised', 25008, 337630, None),
 ('OSA_FG_INCLAVO', 'Sleep apnoea incl. primary care (sensitivity)', 'FinnGen R12 G6_SLEEPAPNO_INCLAVO', 'Finnish', 'hg38', 62722, 437626, None),
 ('T2D_FG', 'Type 2 diabetes (sensitivity)', 'FinnGen R12 T2D', 'Finnish', 'hg38', 82878, 403489, None),
 ('WHRADJBMI_PULIT', 'Waist-to-hip ratio adjusted for BMI (H3; BMI-adjusted)', 'Pulit 2019; GCST008999', 'European', 'hg19', None, None, 'const'),
 # Stage 2 (preregistered 2026-09-26; A1 to A5)
 ('SBP_EVANGELOU', 'Systolic blood pressure (A2; BMI-adjusted)', 'Evangelou 2018; GCST006624', 'European', 'hg19', None, None, 'per-SNP'),
 ('DBP_EVANGELOU', 'Diastolic blood pressure (A2; BMI-adjusted)', 'Evangelou 2018; GCST006630', 'European', 'hg19', None, None, 'per-SNP'),
 ('INCOME_HILL', 'Household income (A4)', 'Hill 2019; GCST009523', 'European', 'hg19', None, None, 'const'),
 ('COGNITION_LEE', 'Cognitive performance (A4, secondary)', 'Lee 2018; GCST006572', 'European', 'hg38 harmonised', None, None, 'const'),
 ('SNORING_CAMPOS', 'Snoring (A5)', 'Campos 2020; GCST009760', 'European', 'hg38 harmonised', 152302, 256015, None),
 ('SNORING_CAMPOS_ADJBMI', 'Snoring adjusted for BMI (A5; BMI-adjusted)', 'Campos 2020; GCST009761', 'European', 'hg38 harmonised', 151836, 255230, None),
 ('INSOMNIA_FG', 'Insomnia (A5)', 'FinnGen R12 F5_INSOMNIA', 'Finnish', 'hg38', 6776, 490763, None),
 ('SLEEPDUR_DASHTI', 'Sleep duration (A5)', 'Dashti 2019; GCST007561', 'European', 'hg19', None, None, 'const'),
 ('CHRONOTYPE_JONES', 'Chronotype, positive = morningness (A5)', 'Jones 2019; GCST007576', 'European', 'hg19', None, None, 'const'),
 ('FEV1FVC_SHRINE', 'FEV1/FVC ratio (A5)', 'Shrine 2019; GCST007431', 'European', 'hg38 harmonised', None, None, 'per-SNP'),
 ('DEPRESSION_FG', 'Depression (A5)', 'FinnGen R12 F5_DEPRESSIO', 'Finnish', 'hg38', 59333, 434831, None),
 ('ASTHMA_FG', 'Asthma (A5, lower-airway contrast)', 'FinnGen R12 J10_ASTHMA_EXMORE', 'Finnish', 'hg38', 59100, 259839, None),
 # Stage 2 addendum 03 (A7, A8): Pan-UKB, not BMI-adjusted, outside FinnGen
 ('HTN_PANUKB', 'Hypertension (A7 primary; A8 outcome)', 'Pan-UKB phecode 401', 'European', 'hg19', 94311, 325488, None),
 ('SBP_PANUKB', 'Systolic blood pressure (A7)', 'Pan-UKB 4080 (irnt)', 'European', 'hg19', None, None, 'const'),
 ('DBP_PANUKB', 'Diastolic blood pressure (A7)', 'Pan-UKB 4079 (irnt)', 'European', 'hg19', None, None, 'const'),
]
for tr, *_r, ca, co, nnote in ROLE:
    c = STDCFG.TRAITS[tr]
    assert (c.get('case'), c.get('ctrl')) == (ca, co), f'{tr}: case/control differs from standardize.py {c.get("case"), c.get("ctrl")}'
rows, inputs = [], []
for tr, label, src, anc, build, ca, co, nnote in ROLE:
    mp = f'{STD}/{tr}.meta.json'
    if not os.path.exists(mp):
        rows.append(dict(trait_id=tr, trait=label, source=src, ancestry=anc, build=build, cases=ca, controls=co,
                         N_effective_or_total=None, raw_rows=None, hapmap3_snps_used=None, raw_file=None, raw_sha256=None, status='not yet standardized'))
        continue
    m = json.load(open(mp))
    raw = m['source']   # path relative to the repository root (data/raw/...)
    inputs.append(raw)
    neff = None if nnote else round(4 / (1 / ca + 1 / co))
    nval = neff if neff else (STDCFG.TRAITS[tr]['n'] if nnote == 'const' else 'per-SNP N in file')
    rows.append(dict(trait_id=tr, trait=label, source=src, ancestry=anc, build=build, cases=ca, controls=co,
                     N_effective_or_total=nval, raw_rows=m['n_raw_rows'],
                     hapmap3_snps_used=m['n_hm3_kept'], raw_file=m['source'], raw_sha256=figkit.sha256(raw), status='standardized'))
df = pd.DataFrame(rows)
checks = [(f"{r.trait_id}: Neff = 4 / (1/{r.cases} + 1/{r.controls})", r.N_effective_or_total)
          for r in df.itertuples() if pd.notna(r.cases) and r.status == 'standardized']
cols = {'trait_id': 'internal identifier used in scripts', 'cases': 'number of cases (binary traits)',
        'controls': 'number of controls', 'N_effective_or_total': 'binary traits: effective sample size 4/(1/Ncase+1/Ncontrol) used in LDSC and Genomic SEM; continuous traits: per-SNP N from the file, or the constant total N from the publication when the file has none',
        'raw_rows': 'variant rows in the raw file', 'hapmap3_snps_used': 'HapMap3 SNPs retained after allele matching, deduplication and QC',
        'raw_sha256': 'SHA256 of the raw downloaded file (see data/raw/SHA256SUMS.txt)', 'build': 'genome build of the raw file; merging is by rsID'}
out = figkit.write_raw_xlsx(f'{TOPIC}/figures/raw_data/Table1_raw.xlsx', df, inputs, os.path.abspath(__file__), cols, checks)
df.drop(columns=['raw_sha256', 'raw_file']).to_csv(f'{TOPIC}/figures/Table1_inputs.tsv', sep='\t', index=False)
print(df[['trait_id', 'cases', 'controls', 'N_effective_or_total', 'hapmap3_snps_used', 'status']].to_string(index=False))
print('wrote', out)
