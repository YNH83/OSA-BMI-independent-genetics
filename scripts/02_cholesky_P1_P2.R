# C9 primary thresholds P1 and P2 (preregistered: analysis_plan/C9_預先登錄_Stage1.md, section 5)
# Cholesky: F_BMI loads on BMI and OSA; F_NonBMI loads on OSA only; factors orthogonal; DWLS.
# P1 = share of OSA genetic variance from F_NonBMI.
# P2 = attenuation of rg with T2D: 1 - rg(F_NonBMI, T2D) / rg(OSA, T2D). T2D = Xue 2018 (preregistered).
# Usage: Rscript 02_cholesky_P1_P2.R [EXT=T2D_XUE] [OSA=OSA_FG] [BMI=BMI_GIANT]
#   OSA may be replaced for sensitivity (OSA_FG_INCLAVO) or for null control 1 (GOUT_FG: gout minus BMI).
source("env/setup_env.R")
args <- commandArgs(trailingOnly = TRUE)
T2D <- if (length(args) >= 1) args[1] else "T2D_XUE"
OSA <- if (length(args) >= 2) args[2] else "OSA_FG"
BMI <- if (length(args) >= 3) args[3] else "BMI_GIANT"
RUN <- if (OSA == "OSA_FG" && BMI == "BMI_GIANT") T2D else paste(T2D, OSA, BMI, sep = "__")

TOPIC <- normalizePath(".")
STD   <- normalizePath("data/hm3_std", mustWork = FALSE)
HM3   <- normalizePath("ref/ldsc_ref/w_hm3.snplist", mustWork = FALSE)
LD    <- paste0(normalizePath("tools/genetics_env/ldsc/eur_ld", mustWork = FALSE), "/")
WLD   <- paste0(normalizePath("tools/genetics_env/ldsc/eur_wld_noMHC25_35", mustWork = FALSE), "/")
setwd(file.path(TOPIC, "derived/munged"))
meta <- function(t) jsonlite::fromJSON(file.path(STD, paste0(t, ".meta.json")))

traits <- c(BMI, OSA, T2D)
for (t in traits) if (!file.exists(paste0(t, ".sumstats.gz")))
  munge(files = file.path(STD, paste0(t, ".tsv.gz")), hm3 = HM3, trait.names = t, info.filter = 0.9, maf.filter = 0.01)
prev <- sapply(traits, function(t) { m <- meta(t); if (isTRUE(m$binary)) m$sample_prev else NA })
sp <- ifelse(is.na(prev), NA, 0.5)

cov <- ldsc(sep_weights = TRUE, traits = paste0(traits, ".sumstats.gz"), sample.prev = sp, population.prev = prev,
            ld = LD, wld = WLD, trait.names = traits, ldsc.log = file.path(TOPIC, paste0("results/02_ldsc_", RUN)))
saveRDS(cov, file.path(TOPIC, paste0("results/02_ldsc_", RUN, ".rds")))

# ---- P1: two-indicator Cholesky (identified with the fixed-variance parameterisation)
m1 <- sprintf("
F_BMI =~ NA*%s + %s
F_NonBMI =~ NA*%s
F_BMI ~~ 1*F_BMI
F_NonBMI ~~ 1*F_NonBMI
F_BMI ~~ 0*F_NonBMI
%s ~~ 0*%s
%s ~~ 0*%s
", BMI, OSA, OSA, BMI, BMI, OSA, OSA)
pick <- function(est, a, b, op, col) {
  v <- est[[col]][est$op == op & ((est$lhs == a & est$rhs == b) | (est$lhs == b & est$rhs == a))]
  if (length(v) != 1) stop(sprintf("parameter %s %s %s not found uniquely", a, op, b))
  as.numeric(v)
}
fit1 <- usermodel(cov, estimation = "DWLS", model = m1, std.lv = FALSE, imp_cov = FALSE)
est1 <- fit1$results
s_bmi <- pick(est1, "F_BMI", OSA, "=~", "STD_Genotype"); s_bmi_se <- pick(est1, "F_BMI", OSA, "=~", "STD_Genotype_SE")
s_non <- pick(est1, "F_NonBMI", OSA, "=~", "STD_Genotype"); s_non_se <- pick(est1, "F_NonBMI", OSA, "=~", "STD_Genotype_SE")
# standardized loadings: s_bmi^2 + s_non^2 = 1, so P1 = s_non^2; delta-method SE = 2 * s_non * SE(s_non)
P1 <- s_non^2; P1_se <- 2 * s_non * s_non_se

m2 <- paste0(m1, "\n", T2D, " ~~ F_BMI + F_NonBMI\n")
fit2 <- usermodel(cov, estimation = "DWLS", model = m2, std.lv = FALSE, imp_cov = FALSE)
est2 <- fit2$results
rg_non <- pick(est2, "F_NonBMI", T2D, "~~", "STD_Genotype"); rg_non_se <- pick(est2, "F_NonBMI", T2D, "~~", "STD_Genotype_SE")
rg_bmi <- pick(est2, "F_BMI", T2D, "~~", "STD_Genotype")
std_cov <- ldsc(sep_weights = TRUE, traits = paste0(traits, ".sumstats.gz"), sample.prev = sp, population.prev = prev, ld = LD, wld = WLD,
                trait.names = traits, stand = TRUE, ldsc.log = file.path(TOPIC, paste0("results/02_ldsc_std_", RUN)))
iO <- which(traits == OSA); iT <- which(traits == T2D)
rg_osa <- std_cov$S_Stand[iO, iT]
k <- length(traits); idx <- matrix(0, k, k); idx[lower.tri(idx, diag = TRUE)] <- seq_len(k * (k + 1) / 2)
rg_osa_se <- sqrt(diag(std_cov$V_Stand))[idx[max(iO, iT), min(iO, iT)]]
P2 <- 1 - rg_non / rg_osa
# delta method treating the two rg estimates as independent (approximation, stated in output)
P2_se <- sqrt((rg_non_se / rg_osa)^2 + (rg_non * rg_osa_se / rg_osa^2)^2)
check <- s_bmi * rg_bmi + s_non * rg_non   # should reproduce the model-implied rg(OSA, T2D)

ci <- function(x, se) sprintf("%.4f (%.4f to %.4f)", x, x - 1.96 * se, x + 1.96 * se)
grade <- function(x, pass, partial) if (x >= pass) "PASS" else if (x >= partial) "PARTIAL" else "FAIL"
out <- data.frame(
  item = c("environment", "T2D source", "rg(OSA,BMI) [std loading F_BMI -> OSA]", "std loading F_NonBMI -> OSA",
           "P1 NonBMI share of OSA genetic variance (95% CI)", "P1 grade (PASS >= 0.60; PARTIAL 0.40 to 0.60)",
           "rg(OSA,T2D) LDSC (95% CI)", "rg(F_BMI,T2D)", "rg(F_NonBMI,T2D) (95% CI)",
           "check: s_bmi*rg_bmi + s_non*rg_non", "P2 attenuation 1 - rg_non/rg_osa (95% CI, delta method, independence approx.)",
           "P2 grade (PASS >= 0.50; PARTIAL 0.25 to 0.50)"),
  value = c(ENV_VERSIONS, T2D, ci(s_bmi, s_bmi_se), ci(s_non, s_non_se), ci(P1, P1_se), grade(P1, 0.60, 0.40),
            ci(rg_osa, rg_osa_se), sprintf("%.4f", rg_bmi), ci(rg_non, rg_non_se), sprintf("%.4f", check),
            ci(P2, P2_se), grade(P2, 0.50, 0.25)))
write.table(out, file.path(TOPIC, paste0("results/02_P1_P2_", RUN, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
num <- data.frame(external = T2D, osa_trait = OSA, bmi_trait = BMI, s_bmi = s_bmi, s_bmi_se = s_bmi_se, s_non = s_non, s_non_se = s_non_se,
                  P1 = P1, P1_se = P1_se, rg_osa_ext = rg_osa, rg_osa_ext_se = rg_osa_se, rg_bmi_ext = rg_bmi,
                  rg_non_ext = rg_non, rg_non_ext_se = rg_non_se, P2 = P2, P2_se = P2_se, check_implied_rg = check,
                  env = ENV_VERSIONS, run_time = format(Sys.time()))
write.table(num, file.path(TOPIC, paste0("results/02_numeric_", RUN, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rbind(cbind(model = "P1_cholesky", est1), cbind(model = "P2_with_T2D", est2)), file.path(TOPIC, paste0("results/02_model_estimates_", RUN, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(fit1 = fit1, fit2 = fit2, std_cov = std_cov), file.path(TOPIC, paste0("results/02_fits_", RUN, ".rds")))
print(out, right = FALSE)
