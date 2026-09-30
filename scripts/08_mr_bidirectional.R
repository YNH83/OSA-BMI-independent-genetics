# C9 Stage 2 A3 (preregistered): bidirectional MR between BMI (GIANT, Yengo 2018) and OSA (FinnGen R12).
# Instruments: P < 5e-8, plink LD clumping r2 < 0.001, 10,000 kb, 1000G Phase3 EUR (MAGMA g1000_eur bfile).
# Palindromic SNPs are removed (standardized files carry MAF, not effect-allele frequency, so strand cannot be inferred).
# Methods: IVW (primary), MR-Egger, weighted median, weighted mode; Cochran Q; Egger intercept; Steiger direction.
source("env/setup_env.R")
suppressMessages(library(TwoSampleMR))
TOPIC <- normalizePath(".")
STD   <- normalizePath("data/hm3_std", mustWork = FALSE)
PLINK <- normalizePath("tools/genetics_env/bin/plink19/plink", mustWork = FALSE)
BFILE <- file.path(normalizePath("tools/genetics_env/magma/g1000_eur", mustWork = FALSE), "g1000_eur")
dir.create(file.path(TOPIC, "derived/mr"), showWarnings = FALSE)
rd <- function(t) data.table::fread(file.path(STD, paste0(t, ".tsv.gz")), data.table = FALSE)
pal <- function(a1, a2) paste0(a1, a2) %in% c("AT", "TA", "CG", "GC")

clump <- function(d, tag) {
  g <- d[d$P < 5e-8 & !pal(d$A1, d$A2), c("SNP", "P")]
  f <- file.path(TOPIC, "derived/mr", paste0(tag, "_gws.txt"))
  write.table(g, f, sep = "\t", quote = FALSE, row.names = FALSE)
  out <- file.path(TOPIC, "derived/mr", paste0(tag, "_clump"))
  # every flag and value is its own element; paths quoted (the project path contains spaces and CJK characters)
  system2(PLINK, c("--bfile", shQuote(BFILE), "--clump", shQuote(f), "--clump-snp-field", "SNP", "--clump-field", "P",
                   "--clump-p1", "5e-8", "--clump-r2", "0.001", "--clump-kb", "10000", "--out", shQuote(out)),
          stdout = FALSE, stderr = FALSE)
  cl <- read.table(paste0(out, ".clumped"), header = TRUE)
  cat(tag, ": GWS non-palindromic", nrow(g), "-> clumped", nrow(cl), "\n")
  cl$SNP
}
fmt <- function(d, snps, type, name, ncase = NULL, nctrl = NULL) {
  x <- d[d$SNP %in% snps, ]
  y <- data.frame(SNP = x$SNP, beta = x$BETA, se = x$SE, effect_allele = x$A1, other_allele = x$A2, pval = x$P,
                  samplesize = x$N, eaf = x$MAF, units = if (is.null(ncase)) "SD" else "log odds")
  # eaf = MAF: used only for variance explained, which depends on p(1-p) and is symmetric in allele orientation;
  # HapMap3 excludes palindromic SNPs, so eaf is never used for strand resolution
  if (!is.null(ncase)) { y$ncase <- ncase; y$ncontrol <- nctrl }
  f <- format_data(y, type = type, snp_col = "SNP", beta_col = "beta", se_col = "se", effect_allele_col = "effect_allele",
                   other_allele_col = "other_allele", pval_col = "pval", samplesize_col = "samplesize", eaf_col = "eaf",
                   units_col = "units", ncase_col = "ncase", ncontrol_col = "ncontrol")
  # prevalence for Steiger on the liability scale (binary traits only); this TwoSampleMR version has no prevalence_col
  f[[paste0("prevalence.", type)]] <- if (is.null(ncase)) NA else ncase / (ncase + nctrl)
  f
}
bmi <- rd("BMI_GIANT"); osa <- rd("OSA_FG")
osa_case <- 56885; osa_ctrl <- 441137
res_all <- list(); inst_all <- list()
for (dir in c("BMI_to_OSA", "OSA_to_BMI")) {
  if (dir == "BMI_to_OSA") {
    s <- clump(bmi, "BMI"); ex <- fmt(bmi, s, "exposure", "BMI"); ou <- fmt(osa, s, "outcome", "OSA", osa_case, osa_ctrl)
  } else {
    s <- clump(osa, "OSA"); ex <- fmt(osa, s, "exposure", "OSA", osa_case, osa_ctrl); ou <- fmt(bmi, s, "outcome", "BMI")
  }
  ex$exposure <- ex$id.exposure <- sub("_to_.*", "", dir); ou$outcome <- ou$id.outcome <- sub(".*_to_", "", dir)
  h <- harmonise_data(ex, ou, action = 3)
  h <- h[h$mr_keep, ]
  st <- tryCatch(directionality_test(h), error = function(e) NULL)
  h <- tryCatch(steiger_filtering(h), error = function(e) { h$steiger_dir <- NA; h })
  set.seed(20260926)   # weighted median and weighted mode SEs are bootstrapped; fixed seed for reproducibility (added 2026-09-27)
  m <- mr(h, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
  m$direction <- dir
  m$lo <- m$b - 1.96 * m$se; m$hi <- m$b + 1.96 * m$se
  het <- mr_heterogeneity(h); pl <- mr_pleiotropy_test(h)
  m$Q_ivw <- het$Q[het$method == "Inverse variance weighted"]; m$Q_p_ivw <- het$Q_pval[het$method == "Inverse variance weighted"]
  m$egger_intercept <- pl$egger_intercept; m$egger_intercept_p <- pl$pval
  m$n_steiger_correct <- sum(h$steiger_dir, na.rm = TRUE); m$n_instruments <- nrow(h)
  if (!is.null(st)) { m$steiger_overall_correct <- st$correct_causal_direction; m$steiger_p <- st$steiger_pval }
  res_all[[dir]] <- m; inst_all[[dir]] <- cbind(direction = dir, h[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
    "beta.exposure", "se.exposure", "pval.exposure", "beta.outcome", "se.outcome", "pval.outcome", "steiger_dir")])
}
R <- do.call(rbind, lapply(res_all, function(x) x[, intersect(c("direction", "method", "nsnp", "b", "se", "pval", "lo", "hi", "Q_ivw", "Q_p_ivw",
  "egger_intercept", "egger_intercept_p", "n_instruments", "n_steiger_correct", "steiger_overall_correct", "steiger_p"), names(x))]))
R$env <- ENV_VERSIONS
b2o <- R[R$direction == "BMI_to_OSA", ]; o2b <- R[R$direction == "OSA_to_BMI", ]
ok_order <- all(b2o$b[b2o$method %in% c("Inverse variance weighted", "Weighted median")] > 0) &&
            all(b2o$pval[b2o$method %in% c("Inverse variance weighted", "Weighted median")] < 0.05)
rev_flag <- all(o2b$pval[o2b$method %in% c("Inverse variance weighted", "Weighted median")] < 0.05) &&
            length(unique(sign(o2b$b[o2b$method %in% c("Inverse variance weighted", "Weighted median")]))) == 1
R$A3_cholesky_order_supported <- ok_order; R$A3_reverse_flag <- rev_flag
# sensitivity (reported, not a preregistered decision rule): IVW restricted to Steiger-concordant instruments
for (dir in names(inst_all)) {
  d <- inst_all[[dir]]; d <- d[isTRUE(TRUE) & d$steiger_dir %in% TRUE, ]
  if (nrow(d) >= 3) {   # same variance model as the primary IVW (TwoSampleMR mr_ivw: multiplicative random effects, SE_lm / min(1, sigma)); corrected 2026-09-27
    f <- summary(lm(d$beta.outcome ~ 0 + d$beta.exposure, weights = 1 / d$se.outcome^2))
    R[R$direction == dir, "ivw_steiger_filtered_b"] <- coef(f)[1, 1]; R[R$direction == dir, "ivw_steiger_filtered_se"] <- coef(f)[1, 2] / min(1, f$sigma)
    R[R$direction == dir, "ivw_steiger_filtered_n"] <- nrow(d) }
}
write.table(R, file.path(TOPIC, "results/08_mr_A3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, inst_all), file.path(TOPIC, "results/08_mr_A3_instruments.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(R[, c("direction", "method", "nsnp", "b", "lo", "hi", "pval")], digits = 3)
cat("Cholesky order supported:", ok_order, "; reverse-effect flag:", rev_flag, "\n")
