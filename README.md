# BMI-independent genetic liability to obstructive sleep apnoea

This repository contains the analysis code for partitioning the genetic liability to obstructive sleep apnoea (OSA) into a BMI-related component (F_BMI) and a BMI-independent component (F_NonBMI) by GWAS-by-subtraction in Genomic SEM. OSA from FinnGen R12 and BMI from GIANT (Yengo 2018) enter a Cholesky model estimated on the LDSC genetic covariance matrix; the BMI-independent share of OSA SNP heritability, the attenuation of genetic correlations with cardiometabolic, sleep, respiratory and psychiatric traits, a multivariable GWAS of F_NonBMI, locus-level triangulation, MAGMA/GTEx tissue enrichment, LAVA local genetic correlation, bidirectional and multivariable Mendelian randomization, colocalisation and a gout minus BMI null control follow from it. All inputs are public GWAS summary statistics; no individual-level data are used or distributed.

## Repository layout

```
analysis_plan/   prespecified analysis plan (Stage 1, Stage 2 and addenda) and deviation records, in Chinese;
                 SHA256_登錄時間戳.txt lists the SHA256 and registration time of each file. Byte-for-byte copies; do not edit.
data_prep/       download, verification and HapMap3 standardization of the raw summary statistics
                 (download_manifest.tsv lists every file with its URL)
scripts/         topic analysis scripts (R and Python), numbered in run order, plus figure and table scripts
shared/          helpers and pipeline steps shared by several scripts (figure kit, chunk merge, MAGMA, coloc,
                 split-half, robustness tables, QQ, enrichment and sensitivity figures, orchestration shells)
env/             R environment: setup_env.R (sourced by every R script) and package installers
```

All scripts are run from the repository root. They create and use these folders, which are not tracked:

```
data/raw/        raw downloads (layout as in data_prep/download_manifest.tsv, e.g. data/raw/FinnGen R12/...)
data/hm3_std/    standardized HapMap3 summary statistics, <TRAIT>.tsv.gz and <TRAIT>.meta.json
ref/             LDSC HapMap3 SNP list
tools/genetics_env/   R library, LDSC LD scores, 1000G reference, MAGMA, PLINK and LAVA support files
derived/         munged sumstats, multivariable GWAS output, LAVA and MR intermediates
results/         model output tables and logs
figures/, tables/   display items (TIFF, PDF, TSV) and raw-data Excel files (figures/raw_data/)
```

## Data sources

Trait codes are the identifiers used in the scripts. FinnGen R12 files are identified by endpoint; GWAS Catalog files by study accession; Pan-UK Biobank files by phenotype code (EUR columns). URLs are in `data_prep/download_manifest.tsv`; per-trait parsing rules and sample sizes are in `data_prep/standardize.py`.

| Trait code | Trait | Source | Accession or endpoint | Role |
|---|---|---|---|---|
| OSA_FG | Sleep apnoea | FinnGen R12 | G6_SLEEPAPNO | primary |
| BMI_GIANT | Body mass index | GIANT and UK Biobank, Yengo 2018 | GCST006900 | primary |
| T2D_XUE | Type 2 diabetes | Xue 2018 | GCST006867 | external correlate (P2) |
| HTN_FG | Essential hypertension | FinnGen R12 | I9_HYPTENSESS | external correlate |
| CAD_ARAGAM | Coronary artery disease | Aragam 2022 | GCST90132314 | external correlate |
| GOUT_FG | Gout, strict definition | FinnGen R12 | GOUT_STRICT | null control |
| OSA_FG_INCLAVO | Sleep apnoea including primary care | FinnGen R12 | G6_SLEEPAPNO_INCLAVO | sensitivity |
| T2D_FG | Type 2 diabetes | FinnGen R12 | T2D | sensitivity |
| OSA_CAMPOS | Sleep apnoea | Campos 2023 | GCST90267361 | consistency check (overlaps FinnGen) |
| OSA_CAMPOS_ADJBMI | Sleep apnoea adjusted for BMI | Campos 2023 | GCST90267363 | triangulation |
| OSA_STRAUSZ_BMILT30 | Sleep apnoea in BMI < 30 | Strausz 2025 | GCST90624797 | triangulation |
| WHRADJBMI_PULIT | Waist-to-hip ratio adjusted for BMI | Pulit 2019 | GCST008999 | H3 (fat distribution) |
| SBP_EVANGELOU | Systolic blood pressure (BMI-adjusted GWAS) | Evangelou 2018 | GCST006624 | Stage 2 A2 |
| DBP_EVANGELOU | Diastolic blood pressure (BMI-adjusted GWAS) | Evangelou 2018 | GCST006630 | Stage 2 A2 |
| INCOME_HILL | Household income | Hill 2019 | GCST009523 | Stage 2 A4 |
| COGNITION_LEE | Cognitive performance | Lee 2018 | GCST006572 | Stage 2 A4 |
| SNORING_CAMPOS | Snoring (main analysis) | Campos 2020 | GCST009760 | Stage 2 A5 |
| SNORING_CAMPOS_ADJBMI | Snoring adjusted for BMI | Campos 2020 | GCST009761 | Stage 2 A5 |
| SLEEPDUR_DASHTI | Sleep duration | Dashti 2019 | GCST007561 | Stage 2 A5 |
| CHRONOTYPE_JONES | Chronotype | Jones 2019 | GCST007576 | Stage 2 A5 |
| FEV1FVC_SHRINE | FEV1/FVC ratio | Shrine 2019 | GCST007431 | Stage 2 A5 |
| INSOMNIA_FG | Insomnia | FinnGen R12 | F5_INSOMNIA | Stage 2 A5 |
| DEPRESSION_FG | Depression | FinnGen R12 | F5_DEPRESSIO | Stage 2 A5 |
| ASTHMA_FG | Asthma | FinnGen R12 | J10_ASTHMA_EXMORE | Stage 2 A5 |
| HTN_PANUKB | Hypertension | Pan-UK Biobank | phecode-401 | Stage 2 A7, A8 |
| SBP_PANUKB | Systolic blood pressure (irnt) | Pan-UK Biobank | continuous-4080 | Stage 2 A7, A8 |
| DBP_PANUKB | Diastolic blood pressure (irnt) | Pan-UK Biobank | continuous-4079 | Stage 2 A7 |

The Pan-UK Biobank phenotype manifest (`phenotype_manifest.tsv.bgz`) is downloaded for the EUR case and control counts.

## Software

- R 4.5.3 with GenomicSEM 0.0.5 and lavaan 0.6-19. lavaan is pinned: lavaan 0.7.x rejects DWLS on continuous data as GenomicSEM 0.0.5 calls it. `env/install_R_packages.R` installs lavaan from CRAN, so install 0.6-19 afterwards, for example `remotes::install_version("lavaan", "0.6-19", lib = "tools/genetics_env/R_library")`.
- R packages: TwoSampleMR 0.7.11, ieugwasr 1.2.0, coloc 5.2.3, MRPRESSO, LAVA, MRlap, cause, data.table, jsonlite (installers in `env/`; LAVA, MRlap and cause are GitHub packages, see `env/install_addendum03.R`).
- MAGMA v1.10 with the 1000 Genomes EUR reference (g1000_eur), NCBI37.3 gene locations and GTEx v8 gene median TPM.
- PLINK 1.9 (LD clumping for MR and CAUSE; LD for regional plots).
- Python 3 with pandas, numpy, scipy, matplotlib, openpyxl and Pillow.

`env/setup_env.R` puts `tools/genetics_env/R_library` first on the library path and fixes `parallel::detectCores()` at 8 (some sandboxed systems return NA, which breaks lavaan); this does not change any estimate.

## External reference files

| Path | Content |
|---|---|
| `ref/ldsc_ref/w_hm3.snplist` | LDSC HapMap3 SNP list (SNP, A1, A2) |
| `tools/genetics_env/ldsc/eur_ld/` | LDSC European LD scores (per chromosome) |
| `tools/genetics_env/ldsc/eur_wld_noMHC25_35/` | LDSC European regression weights with the MHC (chr6 25 to 35 Mb) removed |
| `tools/genetics_env/ref/reference.1000G.EUR.maf.txt` | 1000 Genomes Phase 3 EUR reference with columns SNP CHR BP MAF A1 A2 (hg19; from the .frq and .bim files) for GenomicSEM `sumstats()` and position mapping |
| `tools/genetics_env/magma/` | `magma` binary, `NCBI37.3/NCBI37.3.gene.loc`, `g1000_eur/g1000_eur.{bed,bim,fam}`, `GTEx_v8_gene_median_tpm.gct.gz` |
| `tools/genetics_env/bin/plink19/plink` | PLINK 1.9 binary |
| `tools/genetics_env/src_pkgs/LAVA-main/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile` | LAVA genomic blocks (from the LAVA source) |
| `tools/genetics_env/R_library/` | project R library (created by the installers) |

## Run order

All commands are run from the repository root. Numbered prefixes follow the order of the analysis; `shared/rerun_all_v2.sh`, `shared/gwas_chain_v2.sh` and `shared/gwas_chain_v3_fixP.sh` are the orchestration shells used for the final run.

1. **Data preparation**
   `python3 data_prep/download.py` (or `data_prep/seg_download.py --manifest-ids ...` for throttled servers), `python3 data_prep/verify.py`, then `python3 data_prep/standardize.py ALL`. `data_prep/compare_fixP.py` compares standardized files before and after the P-value underflow fix (deviation record 03).
2. **Kill-switch, LDSC and Cholesky models** (`shared/rerun_all_v2.sh` runs the core set)
   `Rscript scripts/01_killswitch_ldsc.R`;
   `Rscript scripts/02_cholesky_P1_P2.R <EXT> [OSA] [BMI]` for EXT = T2D_XUE, T2D_FG, HTN_FG, CAD_ARAGAM, WHRADJBMI_PULIT, OSA_CAMPOS_ADJBMI, OSA_STRAUSZ_BMILT30 and every Stage 2 trait (HTN_PANUKB, SBP_EVANGELOU, DBP_EVANGELOU, SBP_PANUKB, DBP_PANUKB, INCOME_HILL, COGNITION_LEE, SNORING_CAMPOS, SNORING_CAMPOS_ADJBMI, SLEEPDUR_DASHTI, CHRONOTYPE_JONES, FEV1FVC_SHRINE, DEPRESSION_FG, INSOMNIA_FG, ASTHMA_FG); sensitivity and null runs `T2D_XUE OSA_FG_INCLAVO BMI_GIANT`, `T2D_XUE OSA_CAMPOS BMI_GIANT`, `T2D_XUE GOUT_FG BMI_GIANT`;
   `Rscript scripts/05_ldsc_matrix.R`.
3. **Multivariable GWAS (userGWAS, in chunks)** (`shared/gwas_chain_v2.sh`, `shared/gwas_chain_v3_fixP.sh`)
   `Rscript scripts/03_userGWAS_subtraction.R all 1 <CHUNK> 8` for CHUNK 1 to 8; `Rscript scripts/03b_userGWAS_null_gout.R all 1 <CHUNK> 8` for the gout null control. `scripts/03c_userGWAS_add_recovered.R` is the record of the one-off repair that added the 14 SNPs recovered by the P-underflow fix; it is not meant to be rerun.
4. **Merge and QC**
   `python3 shared/merge_gwas.py derived/gwas userGWAS_ 8 F_BMI F_NonBMI` and `python3 shared/merge_gwas.py derived/gwas null_gout_userGWAS_ 8 F_BMI F_NonBMI`;
   `python3 scripts/09_postprocess_factor_gwas.py all_fixP` and `python3 scripts/09_postprocess_factor_gwas.py all null_gout_userGWAS_`.
   The downstream scripts read the main GWAS under the tag `all_fixP` (the merged run plus the recovered SNPs). On a fresh run from corrected inputs, the merged `userGWAS_F_*_all.tsv.gz` files already contain those SNPs and are used under that tag.
5. **Loci**
   `python3 scripts/04_P3_loci.py all_fixP`; `python3 scripts/tab_04_loci.py`; `Rscript shared/ldsc_factor_h2.R`.
6. **Stage 2: mechanistic panel, two-axis tests, LAVA, MR and MVMR**
   `Rscript scripts/07_twoaxis_test.R <AXIS2> <LABEL> [METAB]` for `HTN_FG A1_HTN`, `CAD_ARAGAM A1b_CAD`, `SBP_EVANGELOU A2_SBP`, `DBP_EVANGELOU A2_DBP`, `HTN_FG posthoc_HTN_vs_T2DFG T2D_FG`, `HTN_PANUKB A7_HTN_PANUKB`, `SBP_PANUKB A7_SBP_PANUKB`, `DBP_PANUKB A7_DBP_PANUKB`;
   `Rscript scripts/08_mr_bidirectional.R`;
   `Rscript scripts/10_mvmr_A8.R`; `Rscript scripts/10b_mr_robust_A8.R cause HTN_PANUKB`, `... cause T2D_XUE`, `... mrlap`;
   `Rscript scripts/11_lava_A6.R prep`, then `Rscript scripts/11_lava_A6.R run <CHUNK> 4` for CHUNK 1 to 4; `python3 scripts/12b_lava_summary_A6.py`;
   `Rscript scripts/17_multireg_A5.R`.
7. **MAGMA and GTEx tissue enrichment**
   `python3 shared/magma_pipeline.py results/magma <LABEL> <GWAS_FILE> SNP P N` for OSA_FG (`data/hm3_std/OSA_FG.tsv.gz`), BMI_GIANT (`data/hm3_std/BMI_GIANT.tsv.gz`), F_NonBMI (`derived/gwas/F_NonBMI_QC.tsv.gz`) and GOUT_NonBMI (`derived/gwas/null_gout_F_NonBMI_QC.tsv.gz`).
8. **Colocalisation**
   `Rscript shared/coloc_loci.R` (needs `figures/Table4_loci.tsv` from step 5).
9. **Robustness**
   `Rscript shared/split_half.R`; `python3 shared/stab_robustness.py`; `python3 shared/build_sensitivity.py`.
10. **Figures and tables**
    `python3 scripts/10_table1_inputs.py`, `scripts/11_table2_prereg.py`, `scripts/12_stage2_summary.py`, `shared/build_stables_ldsc_sem.py`, `scripts/stab_07_null.py`, `scripts/stab_12_mvmr.py`;
    `scripts/fig_01_design_v2.py`, `scripts/fig_02_architecture_v2.py` (before `fig_03_axes_v2.py`, which reads its raw-data file), `scripts/fig_03_axes_v2.py`, `scripts/fig_04_miami_v2.py`, `scripts/fig_05_triangulation.py`, `shared/fig_enrichment_v2.py`, `shared/fig_qq.py`, `scripts/sfig_02_regional.py`, `scripts/fig_08_twoaxis_mr.py`.

## Analysis plan

`analysis_plan/` holds the prespecified analysis plan (Stage 1, its supplement, Stage 2 and addenda 01 to 04) and the deviation records 01 to 06, written in Chinese. `SHA256_登錄時間戳.txt` lists the SHA256 and the time stamp at which each file was registered, and whether any result existed at that time. The files are kept byte-for-byte so the hashes can be checked with `shasum -a 256 analysis_plan/*.md`. Scripts refer to thresholds and deviation records by these file names.

## Citation

Manuscript under review.
