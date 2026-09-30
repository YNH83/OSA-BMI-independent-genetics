#!/bin/bash
# v2 GWAS chain: NonBMI-OSA multivariable GWAS, then the gout null-control GWAS. Run from the repository root after
# shared/rerun_all_v2.sh has finished (the main GWAS needs results/02_ldsc_T2D_FG.rds and the gout GWAS needs
# results/02_ldsc_T2D_XUE__GOUT_FG__BMI_GIANT.rds). Each step runs 8 single-thread chunks; the main GWAS is then merged.
# Old chunk outputs are archived first, never overwritten.
C9="."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
chunks() { local dir="$1" script="$2" prefix="$3"
  mkdir -p "$dir/derived/gwas/_v1_sepweightsFALSE"
  mv "$dir"/derived/gwas/${prefix}*_all_chunk* "$dir/derived/gwas/_v1_sepweightsFALSE/" 2>/dev/null
  echo "$(date '+%H:%M') start $script"
  for c in 1 2 3 4 5 6 7 8; do (cd "$dir" && Rscript "scripts/$script" all 1 $c 8 > "results/$(basename $script .R)_chunk${c}of8.log" 2>&1) & done; wait
  n=$(ls "$dir"/derived/gwas/${prefix}*_all_chunk*of08.tsv.gz 2>/dev/null | wc -l | tr -d ' ')
  echo "$(date '+%H:%M') done $script: $n chunk files (expect 16)"; }

chunks "$C9" 03_userGWAS_subtraction.R userGWAS_F_
(cd "$C9" && python3 shared/merge_gwas.py derived/gwas userGWAS_ 8 F_BMI F_NonBMI) && echo "C9 merged"
chunks "$C9" 03b_userGWAS_null_gout.R null_gout_userGWAS_F_
echo "GWAS_CHAIN_V2_DONE $(date '+%H:%M')"
