# C9 Stage 2 A1 / A1b / A2 (preregistered: analysis_plan/C9_預先登錄_Stage2_代謝軸對高血壓軸.md)
# Test cor(F_NonBMI, F_METAB) = cor(F_NonBMI, F_AXIS2), where METAB = T2D (Xue 2018) and AXIS2 = HTN / CAD / SBP / DBP.
# T2D and AXIS2 enter as single-indicator factors with unit variance and zero residual, so all factor covariances are
# genetic correlations. The unconstrained model is saturated (df 0); the chi-square of the constrained model (df 1)
# is therefore the chi-square difference test.
# Usage: Rscript 07_twoaxis_test.R AXIS2 [LABEL] [METAB=T2D_XUE]
#   METAB = T2D_FG gives a same-registry comparison with FinnGen HTN (post-hoc sensitivity, not preregistered).
source("env/setup_env.R")
a <- commandArgs(trailingOnly = TRUE)
AX <- a[1]; LABEL <- if (length(a) >= 2) a[2] else AX
METAB <- if (length(a) >= 3) a[3] else "T2D_XUE"
TOPIC <- normalizePath(".")
STD   <- normalizePath("data/hm3_std", mustWork = FALSE)
HM3   <- normalizePath("ref/ldsc_ref/w_hm3.snplist", mustWork = FALSE)
LD    <- paste0(normalizePath("tools/genetics_env/ldsc/eur_ld", mustWork = FALSE), "/")
WLD   <- paste0(normalizePath("tools/genetics_env/ldsc/eur_wld_noMHC25_35", mustWork = FALSE), "/")
setwd(file.path(TOPIC, "derived/munged"))
meta <- function(t) jsonlite::fromJSON(file.path(STD, paste0(t, ".meta.json")))
traits <- c("BMI_GIANT", "OSA_FG", METAB, AX)
for (t in traits) if (!file.exists(paste0(t, ".sumstats.gz")))
  munge(files = file.path(STD, paste0(t, ".tsv.gz")), hm3 = HM3, trait.names = t, info.filter = 0.9, maf.filter = 0.01)
prev <- sapply(traits, function(t) { m <- meta(t); if (isTRUE(m$binary)) m$sample_prev else NA })
cov <- ldsc(sep_weights = TRUE, traits = paste0(traits, ".sumstats.gz"), sample.prev = ifelse(is.na(prev), NA, 0.5),
            population.prev = prev, ld = LD, wld = WLD, trait.names = traits,
            ldsc.log = file.path(TOPIC, paste0("results/07_ldsc_twoaxis_", LABEL)))
# Marker loadings are labelled and constrained > 0: single-indicator factors are otherwise sign-indeterminate, and a
# negative marker loading flips the sign of the unstandardized covariance that the equality constraint acts on
# (observed 2026-09-26 in the first A1b run: T2D loading -0.31 produced a spurious chi-square of 28.6).
base <- sprintf("
F_BMI =~ NA*BMI_GIANT + lb*BMI_GIANT + OSA_FG
F_NonBMI =~ NA*OSA_FG + ln*OSA_FG
F_METAB =~ NA*%s + lm*%s
F_AXIS2 =~ NA*%s + la*%s
lb > 0.001
ln > 0.001
lm > 0.001
la > 0.001
F_BMI ~~ 1*F_BMI
F_NonBMI ~~ 1*F_NonBMI
F_METAB ~~ 1*F_METAB
F_AXIS2 ~~ 1*F_AXIS2
F_BMI ~~ 0*F_NonBMI
BMI_GIANT ~~ 0*BMI_GIANT
OSA_FG ~~ 0*OSA_FG
%s ~~ 0*%s
%s ~~ 0*%s
F_METAB ~~ F_BMI + F_AXIS2
F_AXIS2 ~~ F_BMI
", METAB, METAB, AX, AX, METAB, METAB, AX, AX)
free <- paste0(base, "F_METAB ~~ F_NonBMI\nF_AXIS2 ~~ F_NonBMI\n")
cons <- paste0(base, "F_METAB ~~ eq*F_NonBMI\nF_AXIS2 ~~ eq*F_NonBMI\n")
f1 <- usermodel(cov, estimation = "DWLS", model = free, std.lv = FALSE, imp_cov = FALSE)
f0 <- usermodel(cov, estimation = "DWLS", model = cons, std.lv = FALSE, imp_cov = FALSE)
e <- f1$results
g <- function(x, y, col) as.numeric(e[[col]][e$op == "~~" & ((e$lhs == x & e$rhs == y) | (e$lhs == y & e$rhs == x))])
rM <- g("F_METAB", "F_NonBMI", "STD_Genotype"); sM <- g("F_METAB", "F_NonBMI", "STD_Genotype_SE")
rA <- g("F_AXIS2", "F_NonBMI", "STD_Genotype"); sA <- g("F_AXIS2", "F_NonBMI", "STD_Genotype_SE")
chisq <- as.numeric(f0$modelfit$chisq); df <- as.numeric(f0$modelfit$df); p <- as.numeric(f0$modelfit$p_chisq)
grade <- if (rA > rM && p < 0.05) "PASS (two axes)" else "FAIL"
out <- data.frame(label = LABEL, metab = METAB, axis2 = AX,
                  cor_NonBMI_metab = rM, se_metab = sM, cor_NonBMI_axis2 = rA, se_axis2 = sA,
                  cor_BMI_metab = g("F_METAB", "F_BMI", "STD_Genotype"), cor_BMI_axis2 = g("F_AXIS2", "F_BMI", "STD_Genotype"),
                  cor_metab_axis2 = g("F_METAB", "F_AXIS2", "STD_Genotype"),
                  free_model_df = as.numeric(f1$modelfit$df), constrained_chisq = chisq, constrained_df = df, p_difference = p,
                  grade = grade, env = ENV_VERSIONS, run_time = format(Sys.time()))
write.table(out, file.path(TOPIC, paste0("results/07_twoaxis_", LABEL, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rbind(cbind(model = "free", e), cbind(model = "equal", f0$results)),
            file.path(TOPIC, paste0("results/07_twoaxis_estimates_", LABEL, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
print(t(out))
