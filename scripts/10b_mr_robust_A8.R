# C9 Stage 2 addendum 03, A8 robustness (preregistered 2026-09-27 02:32):
#  CAUSE for OSA -> HTN (Pan-UKB, primary) and OSA -> T2D (Xue): does a causal model fit better than shared (correlated
#  pleiotropy) only? LD pruning of CAUSE variants with local plink (r2 < 0.01, 500 kb, 1000G EUR) via ieugwasr::ld_clump.
#  MRlap for BMI -> OSA and BMI -> HTN (Pan-UKB) / T2D (Xue): corrects for sample overlap (GIANT includes UKB), weak
#  instruments and winner's curse using cross-trait LDSC.
# Usage: Rscript 10b_mr_robust_A8.R [cause|mrlap|all]
source("env/setup_env.R")
A_ <- commandArgs(trailingOnly = TRUE); MODE <- A_[1]; if (is.na(MODE)) MODE <- "all"
CAUSE_OUT <- if (length(A_) >= 2) A_[-1] else c("HTN_PANUKB", "T2D_XUE")   # optional: outcomes for CAUSE
TOPIC <- normalizePath(".")
STD   <- normalizePath("data/hm3_std", mustWork = FALSE)
PLINK <- normalizePath("tools/genetics_env/bin/plink19/plink", mustWork = FALSE); BFILE <- file.path(normalizePath("tools/genetics_env/magma/g1000_eur", mustWork = FALSE), "g1000_eur")
LDDIR <- paste0(normalizePath("tools/genetics_env/ldsc/eur_ld", mustWork = FALSE), "/")   # MRlap reads l2 and M_5_50 from one folder (it has no separate weights or MHC option)
HM3   <- normalizePath("ref/ldsc_ref/w_hm3.snplist", mustWork = FALSE)
rd <- function(t) data.table::fread(file.path(STD, paste0(t, ".tsv.gz")), data.table = FALSE)
dir.create(file.path(TOPIC, "derived/mr"), showWarnings = FALSE)
set.seed(20260927)
if (MODE %in% c("cause", "all")) {
  suppressMessages(library(cause))
  osa <- rd("OSA_FG")
  res <- list()
  for (o in CAUSE_OUT) {
    y <- rd(o)
    X <- gwas_merge(osa, y, snp_name_cols = c("SNP", "SNP"), beta_hat_cols = c("BETA", "BETA"), se_cols = c("SE", "SE"),
                    A1_cols = c("A1", "A1"), A2_cols = c("A2", "A2"), pval_cols = c("P", "P"))
    vs <- sample(X$snp, min(1e6, nrow(X)))
    params <- est_cause_params(X, vs)
    top <- data.frame(rsid = X$snp, pval = X$p1)
    top <- top[top$pval < 1e-3, ]
    cl <- ieugwasr::ld_clump(top, clump_kb = 500, clump_r2 = 0.01, clump_p = 1e-3, plink_bin = PLINK, bfile = BFILE)
    fit <- cause(X = X, variants = cl$rsid, param_ests = params)
    s <- summary(fit)
    tab <- as.data.frame(s$tab); tab$outcome <- o; tab$n_variants <- nrow(cl)
    q <- as.data.frame(s$quants[[2]]); q$outcome <- o          # posterior medians and CIs of the causal model
    res[[o]] <- list(tab = tab, quants = q, p = s$p)
    cat(o, ": CAUSE causal vs sharing one-sided p =", s$p, "\n"); print(s$tab)
  }
  data.table::fwrite(do.call(rbind, lapply(res, `[[`, "tab")), file.path(TOPIC, paste0("results/10b_cause_A8_elpd_", paste(CAUSE_OUT, collapse = "_"), ".tsv")), sep = "\t")
  data.table::fwrite(do.call(rbind, lapply(res, `[[`, "quants")), file.path(TOPIC, paste0("results/10b_cause_A8_posterior_", paste(CAUSE_OUT, collapse = "_"), ".tsv")), sep = "\t")
  data.table::fwrite(data.frame(outcome = names(res), p_causal_vs_sharing = sapply(res, `[[`, "p")), file.path(TOPIC, paste0("results/10b_cause_A8_p_", paste(CAUSE_OUT, collapse = "_"), ".tsv")), sep = "\t")
}
if (MODE %in% c("mrlap", "all")) {
  suppressMessages(library(MRlap))
  pos <- data.table::fread(normalizePath("tools/genetics_env/ref/reference.1000G.EUR.maf.txt", mustWork = FALSE), select = c("SNP", "CHR", "BP"), data.table = FALSE)
  fmt <- function(t) { x <- merge(rd(t), pos, by = "SNP")
    data.frame(snp = x$SNP, chr = x$CHR, pos = x$BP, a1 = x$A1, a2 = x$A2, beta = x$BETA, se = x$SE, N = x$N) }
  bmi <- fmt("BMI_GIANT"); out <- list()
  owd <- setwd(file.path(TOPIC, "derived/mr"))
  for (o in c("OSA_FG", "HTN_PANUKB", "T2D_XUE")) {
    if (!file.exists(file.path(STD, paste0(o, ".tsv.gz")))) { cat("skip", o, "\n"); next }
    r <- tryCatch(MRlap(exposure = bmi, exposure_name = "BMI", outcome = fmt(o), outcome_name = o, ld = LDDIR, hm3 = HM3,
                        MR_plink = PLINK, MR_bfile = BFILE, save_logfiles = FALSE, verbose = FALSE),
                  error = function(e) { cat(o, "MRlap error:", conditionMessage(e), "\n"); NULL })
    if (is.null(r)) next
    saveRDS(r, file.path(TOPIC, "derived/mr", paste0("mrlap_BMI_", o, ".rds")))
    m <- r$MRcorrection
    g <- function(L, k) if (!is.null(L[[k]])) L[[k]] else NA
    out[[o]] <- data.frame(exposure = "BMI", outcome = o, observed_b = g(m, "observed_effect"), observed_se = g(m, "observed_effect_se"),
                           corrected_b = g(m, "corrected_effect"), corrected_se = g(m, "corrected_effect_se"), corrected_p = g(m, "corrected_effect_p"),
                           diff_stat = g(m, "test_difference"), diff_p = g(m, "p_difference"), n_IVs = length(g(m, "IVs")),
                           rg = g(r$LDSC, "rg"), int_crosstrait = g(r$LDSC, "int_crosstrait"), int_crosstrait_se = g(r$LDSC, "int_crosstrait_se"),
                           scale = "both traits standardized (MRlap)")
    print(out[[o]])
  }
  setwd(owd)
  if (length(out)) data.table::fwrite(do.call(rbind, out), file.path(TOPIC, "results/10b_mrlap_A8.tsv"), sep = "\t")
}
cat("done", ENV_VERSIONS, "\n")
