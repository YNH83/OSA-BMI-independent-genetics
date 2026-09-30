# C9 Stage 2 addendum 03, A8 (preregistered 2026-09-27 02:32): multivariable MR with OSA (FinnGen R12) and BMI (GIANT) as
# joint exposures; outcomes HTN (Pan-UKB, primary), SBP (Pan-UKB), T2D (Xue), CAD (Aragam), HTN (FinnGen, sensitivity only).
# Instruments: union of P < 5e-8 SNPs of either exposure, palindromic removed, clumped on the smaller P (r2 < 0.001, 10,000 kb,
# 1000G EUR). MVMR-IVW = weighted regression of outcome betas on both exposure betas without intercept, residual SE scaled
# (as TwoSampleMR mv_multiple). Conditional F statistics per Sanderson 2021 with zero exposure covariance (two samples).
# Univariable OSA -> outcome robustness: MR-PRESSO (outlier-corrected), Steiger. CAUSE and MRlap run in 10b when installed.
source("env/setup_env.R")
suppressMessages({ library(TwoSampleMR); library(MRPRESSO) })
TOPIC <- normalizePath(".")
STD   <- normalizePath("data/hm3_std", mustWork = FALSE)
PLINK <- normalizePath("tools/genetics_env/bin/plink19/plink", mustWork = FALSE); BFILE <- file.path(normalizePath("tools/genetics_env/magma/g1000_eur", mustWork = FALSE), "g1000_eur")
OUTS  <- c(HTN_PANUKB = "primary", SBP_PANUKB = "secondary", T2D_XUE = "secondary", CAD_ARAGAM = "secondary", HTN_FG = "sensitivity (same registry as OSA)")
rd <- function(t) data.table::fread(file.path(STD, paste0(t, ".tsv.gz")), data.table = FALSE)
pal <- function(a1, a2) paste0(a1, a2) %in% c("AT", "TA", "CG", "GC")
osa <- rd("OSA_FG"); bmi <- rd("BMI_GIANT")
# ---- instruments ----
m <- merge(osa[, c("SNP", "A1", "A2", "BETA", "SE", "P", "MAF")], bmi[, c("SNP", "A1", "A2", "BETA", "SE", "P")], by = "SNP", suffixes = c(".osa", ".bmi"))
flip <- m$A1.bmi == m$A2.osa & m$A2.bmi == m$A1.osa
m$BETA.bmi[flip] <- -m$BETA.bmi[flip]
m <- m[(m$A1.bmi == m$A1.osa & m$A2.bmi == m$A2.osa) | flip, ]
m$Pmin <- pmin(m$P.osa, m$P.bmi)
g <- m[m$Pmin < 5e-8 & !pal(m$A1.osa, m$A2.osa), ]
dir.create(file.path(TOPIC, "derived/mr"), showWarnings = FALSE)
f <- file.path(TOPIC, "derived/mr/mvmr_union_gws.txt"); write.table(data.frame(SNP = g$SNP, P = g$Pmin), f, sep = "\t", quote = FALSE, row.names = FALSE)
out <- file.path(TOPIC, "derived/mr/mvmr_union_clump")
system2(PLINK, c("--bfile", shQuote(BFILE), "--clump", shQuote(f), "--clump-snp-field", "SNP", "--clump-field", "P", "--clump-p1", "5e-8",
                 "--clump-r2", "0.001", "--clump-kb", "10000", "--out", shQuote(out)), stdout = FALSE, stderr = FALSE)
iv <- g[g$SNP %in% read.table(paste0(out, ".clumped"), header = TRUE)$SNP, ]
cat("instruments: union GWS", nrow(g), "-> clumped", nrow(iv), "(OSA-GWS", sum(iv$P.osa < 5e-8), "; BMI-GWS", sum(iv$P.bmi < 5e-8), ")\n")
# conditional F (Sanderson 2021), exposure covariance 0
condF <- function(bx, sx, bo, so) { d <- coef(lm(bx ~ 0 + bo)); Q <- sum((bx - d * bo)^2 / (sx^2 + d^2 * so^2)); Q / (length(bx) - 1) }
F_osa <- condF(iv$BETA.osa, iv$SE.osa, iv$BETA.bmi, iv$SE.bmi); F_bmi <- condF(iv$BETA.bmi, iv$SE.bmi, iv$BETA.osa, iv$SE.osa)
res <- list(); uni <- list()
for (o in names(OUTS)) {
  fo <- file.path(STD, paste0(o, ".tsv.gz"))
  if (!file.exists(fo)) { cat("skip", o, "(not standardized yet)\n"); next }
  y <- rd(o)[, c("SNP", "A1", "A2", "BETA", "SE", "P")]
  d <- merge(iv, y, by = "SNP")
  fl <- d$A1 == d$A2.osa & d$A2 == d$A1.osa; d$BETA[fl] <- -d$BETA[fl]
  d <- d[(d$A1 == d$A1.osa & d$A2 == d$A2.osa) | fl, ]
  w <- 1 / d$SE^2
  fit <- lm(d$BETA ~ 0 + d$BETA.osa + d$BETA.bmi, weights = w)
  s <- summary(fit); sig <- s$sigma
  cf <- coef(s); se <- cf[, 2] / min(1, sig)                       # multiplicative random effects: lm SE (= sigma x fixed SE) kept when sigma > 1, fixed-effect SE when sigma < 1 (TwoSampleMR mv_multiple). Corrected 2026-09-27, deviation record 04
  for (k in 1:2) {
    b <- cf[k, 1]; z <- b / se[k]
    res[[paste(o, k)]] <- data.frame(outcome = o, role = OUTS[[o]], exposure = c("OSA", "BMI")[k], method = "MVMR-IVW", nsnp = nrow(d),
                                     b = b, se = se[k], lo = b - 1.96 * se[k], hi = b + 1.96 * se[k], p = 2 * pnorm(-abs(z)),
                                     condF_OSA = F_osa, condF_BMI = F_bmi, residual_sigma = sig)
  }
  # univariable OSA -> outcome on OSA-GWS instruments of the same set, with MR-PRESSO
  u <- d[d$P.osa < 5e-8, ]
  if (nrow(u) >= 4) {
    ivw <- summary(lm(u$BETA ~ 0 + u$BETA.osa, weights = 1 / u$SE^2)); b <- coef(ivw)[1, 1]; se1 <- coef(ivw)[1, 2] / min(1, ivw$sigma)   # as above (TwoSampleMR mr_ivw)
    pr <- tryCatch(mr_presso(BetaOutcome = "BETA", BetaExposure = "BETA.osa", SdOutcome = "SE", SdExposure = "SE.osa", OUTLIERtest = TRUE,
                             DISTORTIONtest = TRUE, data = u, NbDistribution = 3000, SignifThreshold = 0.05), error = function(e) NULL)
    mp <- if (!is.null(pr)) pr$`Main MR results` else NULL
    uni[[o]] <- data.frame(outcome = o, exposure = "OSA", nsnp = nrow(u), ivw_b = b, ivw_se = se1, ivw_p = 2 * pnorm(-abs(b / se1)),
                           presso_global_p = if (!is.null(pr)) pr$`MR-PRESSO results`$`Global Test`$Pvalue else NA,
                           presso_corrected_b = if (!is.null(mp)) mp$`Causal Estimate`[2] else NA,
                           presso_corrected_p = if (!is.null(mp)) mp$`P-value`[2] else NA,
                           presso_n_outliers = if (!is.null(pr) && !is.null(pr$`MR-PRESSO results`$`Distortion Test`)) length(pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`) else 0)
  }
}
R <- do.call(rbind, res); U <- do.call(rbind, uni)
R$env <- ENV_VERSIONS; R$run_time <- format(Sys.time())
write.table(R, file.path(TOPIC, "results/10_mvmr_A8.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(U, file.path(TOPIC, "results/10_mvmr_A8_univariable_OSA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(iv, file.path(TOPIC, "results/10_mvmr_A8_instruments.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# preregistered reading
h <- R[R$outcome == "HTN_PANUKB" & R$exposure == "OSA", ]; t2 <- R[R$outcome == "T2D_XUE" & R$exposure == "OSA", ]
if (nrow(h)) {
  direct <- h$b > 0 && h$p < 0.05
  contrast <- direct && nrow(t2) && t2$lo <= 0 && t2$hi >= 0
  cat("A8 OSA direct effect on HTN (Pan-UKB):", ifelse(direct, "SUPPORTED", "NOT supported"), "; axis contrast:", ifelse(contrast, "SUPPORTED", "not supported"),
      "; conditional F OSA", round(F_osa, 1), ifelse(F_osa < 10, "(WEAK)", ""), "\n")
}
print(R[, c("outcome", "exposure", "nsnp", "b", "lo", "hi", "p")], digits = 3); print(U, digits = 3)
