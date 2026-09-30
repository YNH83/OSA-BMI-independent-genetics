#!/bin/bash
# v3 chain after deviation record 03 (P-underflow fix): the gout null-control GWAS. Run from the repository root.
# The F_NonBMI GWAS is not rerun (LDSC unchanged; the 14 recovered SNPs were added by scripts/03c_userGWAS_add_recovered.R).
C9="."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
chunks() { local dir="$1" script="$2"
  echo "$(date '+%H:%M') start $script"
  for c in 1 2 3 4 5 6 7 8; do (cd "$dir" && Rscript "scripts/$script" all 1 $c 8 > "results/$(basename $script .R)_chunk${c}of8.log" 2>&1) & done; wait
  echo "$(date '+%H:%M') done $script"; }
# first build the sumstats cache once (single process) so the 8 chunk processes do not race to create it
(cd "$C9" && Rscript scripts/03b_userGWAS_null_gout.R 1 1 > results/03b_prebuild_sumstats.log 2>&1) && echo "gout sumstats built"
chunks "$C9" 03b_userGWAS_null_gout.R
ls "$C9"/derived/gwas/null_gout_userGWAS_F_*_all_chunk*of08.tsv.gz | wc -l | xargs echo "gout chunk files (expect 16):"
echo "GWAS_CHAIN_V3_DONE $(date '+%H:%M')"
