# C9 Stage 2 addendum 04 S2 (preregistered 2026-09-27 17:11): standardized multiple regression of F_NonBMI on snoring, insomnia and depression.
# Implementation: the regression form (F_NonBMI ~ X1 + X2 + X3 with the residual variance fixed to 1) returned standardized coefficients
# inconsistent with the implied correlations (R2 > 0.9 against a zero-order maximum of 0.27), so the same estimand is computed from the
# joint correlation model: beta = Rxx^-1 rxy, where Rxx are the genetic correlations among the three traits and rxy their correlations with
# F_NonBMI, all from one Genomic SEM fit of the primary Cholesky with the three traits added. SE and 95% CI: 1,000 parametric draws of the
# LDSC S matrix from N(S, V), refit, recompute (percentile CI). Sensitivity: BMI-adjusted snoring. Output results/17_multireg_A5.tsv.
source("env/setup_env.R")
T   <- normalizePath(".")
STD <- normalizePath("data/hm3_std", mustWork = FALSE)
LD  <- paste0(normalizePath("tools/genetics_env/ldsc/eur_ld", mustWork = FALSE), "/"); WLD <- paste0(normalizePath("tools/genetics_env/ldsc/eur_wld_noMHC25_35", mustWork = FALSE), "/")
setwd(file.path(T, "derived/munged"))
meta <- function(t) jsonlite::fromJSON(file.path(STD, paste0(t, ".meta.json")))
NDRAW <- 1000; set.seed(20260927)
res <- list()
for (SN in c("SNORING_CAMPOS", "SNORING_CAMPOS_ADJBMI")) {
  traits <- c("BMI_GIANT", "OSA_FG", SN, "INSOMNIA_FG", "DEPRESSION_FG"); X <- c(SN, "INSOMNIA_FG", "DEPRESSION_FG")
  prev <- sapply(traits, function(t) { m <- meta(t); if (isTRUE(m$binary)) m$sample_prev else NA }); sp <- ifelse(is.na(prev), NA, 0.5)
  cov <- ldsc(sep_weights = TRUE, traits = paste0(traits, ".sumstats.gz"), sample.prev = sp, population.prev = prev, ld = LD, wld = WLD,
              trait.names = traits, ldsc.log = file.path(T, paste0("results/17_multireg_", SN)))
  m <- sprintf("
F_BMI =~ NA*BMI_GIANT + OSA_FG
F_NonBMI =~ NA*OSA_FG
F_BMI ~~ 1*F_BMI
F_NonBMI ~~ 1*F_NonBMI
F_BMI ~~ 0*F_NonBMI
BMI_GIANT ~~ 0*BMI_GIANT
OSA_FG ~~ 0*OSA_FG
F_NonBMI ~~ %s + INSOMNIA_FG + DEPRESSION_FG
F_BMI ~~ %s + INSOMNIA_FG + DEPRESSION_FG
%s ~~ INSOMNIA_FG + DEPRESSION_FG
INSOMNIA_FG ~~ DEPRESSION_FG
", SN, SN, SN)
  beta_of <- function(cv) {
    r <- quiet(usermodel(cv, estimation = "DWLS", model = m, std.lv = FALSE, imp_cov = FALSE))$results
    g <- function(a, b) as.numeric(r$STD_Genotype[r$op == "~~" & ((r$lhs == a & r$rhs == b) | (r$lhs == b & r$rhs == a))])
    Rxx <- outer(X, X, Vectorize(function(a, b) if (a == b) 1 else g(a, b))); rxy <- sapply(X, function(a) g("F_NonBMI", a))
    b <- solve(Rxx, rxy); c(b, R2 = sum(b * rxy), rxy)
  }
  quiet <- function(expr) { out <- NULL; invisible(capture.output(out <- suppressWarnings(suppressMessages(expr)))); out }
  est <- beta_of(cov)
  # parametric draws of the vectorised S (lower triangle) from N(s, V)
  k <- length(traits); lt <- lower.tri(cov$S, diag = TRUE); s <- cov$S[lt]; Lc <- t(chol(cov$V + diag(1e-12, nrow(cov$V))))
  draws <- t(sapply(seq_len(NDRAW), function(i) {
    cv <- cov; S <- matrix(0, k, k); S[lt] <- s + as.vector(Lc %*% rnorm(length(s))); S <- S + t(S) - diag(diag(S))
    dimnames(S) <- dimnames(cov$S); cv$S <- S
    tryCatch(beta_of(cv), error = function(e) rep(NA, length(est)))
  }))
  ok <- stats::complete.cases(draws)
  lab <- c(paste0("beta ", X), "R2", paste0("zero-order rg with F_NonBMI: ", X))
  res[[SN]] <- data.frame(snoring = SN, quantity = lab, estimate = est, se = apply(draws[ok, ], 2, sd),
                          lo = apply(draws[ok, ], 2, quantile, 0.025), hi = apply(draws[ok, ], 2, quantile, 0.975), n_draws_ok = sum(ok))
  print(res[[SN]][, -1], digits = 3)
}
out <- do.call(rbind, res); out$env <- ENV_VERSIONS
write.table(out, file.path(T, "results/17_multireg_A5.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
