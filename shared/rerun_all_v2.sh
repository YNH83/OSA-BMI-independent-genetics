#!/bin/bash
# v2 rerun after the audit (sep_weights = TRUE, MHC chr6 25 to 35 Mb excluded from LDSC weights). Run from the repository root.
# Models first; stops on the first failure so no display item is built from incomplete results.
set -u
C9="."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
run() { local d="$1"; shift; local log="$1"; shift
  echo "$(date '+%H:%M:%S') START $log"; (cd "$d" && "$@" > "results/$log.log" 2>&1) || { echo "FAILED $log"; exit 1; }
  echo "$(date '+%H:%M:%S') OK    $log"; }

# ---------------- models
run "$C9" 01_run Rscript scripts/01_killswitch_ldsc.R
for e in T2D_XUE T2D_FG HTN_FG CAD_ARAGAM WHRADJBMI_PULIT OSA_CAMPOS_ADJBMI OSA_STRAUSZ_BMILT30; do
  run "$C9" "02_run_$e" Rscript scripts/02_cholesky_P1_P2.R "$e"; done
run "$C9" 02_run_T2D_XUE__OSA_FG_INCLAVO__BMI_GIANT Rscript scripts/02_cholesky_P1_P2.R T2D_XUE OSA_FG_INCLAVO BMI_GIANT
run "$C9" 02_run_T2D_XUE__OSA_CAMPOS__BMI_GIANT Rscript scripts/02_cholesky_P1_P2.R T2D_XUE OSA_CAMPOS BMI_GIANT
run "$C9" 02_run_T2D_XUE__GOUT_FG__BMI_GIANT Rscript scripts/02_cholesky_P1_P2.R T2D_XUE GOUT_FG BMI_GIANT
run "$C9" 05_run Rscript scripts/05_ldsc_matrix.R
# ---------------- display items
for s in scripts/10_table1_inputs.py scripts/11_table2_prereg.py; do
  python3 "$s" > /dev/null 2>&1 && echo "display OK $(basename "$s")" || echo "display FAILED $s"; done
for s in build_stables_ldsc_sem.py build_sensitivity.py fig_qq.py; do
  python3 "shared/$s" > /dev/null 2>&1 && echo "display OK $s" || echo "display FAILED $s"; done
echo "RERUN_V2_DONE $(date '+%H:%M:%S')"
