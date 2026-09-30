# Odd/even chromosome split-half of the Cholesky subtraction (C9 addendum 04 S1).
# Same LDSC settings as the primary analyses (sep_weights = TRUE, MHC excluded via the weights folder), select = "ODD" / "EVEN".
# C9 : BMI_GIANT, OSA_FG, T2D_XUE   -> P1 = std loading(F_NonBMI -> OSA)^2 ; P2 = 1 - rg(F_NonBMI, T2D) / rg(OSA, T2D)
# Usage (from the repository root): Rscript shared/split_half.R    Output: results/16_split_half.tsv
source("env/setup_env.R")
topic <- "C9"
STD <- normalizePath("data/hm3_std", mustWork = FALSE)
LD  <- paste0(normalizePath("tools/genetics_env/ldsc/eur_ld", mustWork = FALSE), "/"); WLD <- paste0(normalizePath("tools/genetics_env/ldsc/eur_wld_noMHC25_35", mustWork = FALSE), "/")
T <- normalizePath("."); X <- "BMI_GIANT"; Y <- "OSA_FG"; F1 <- "F_BMI"; F2 <- "F_NonBMI"
E <- "T2D_XUE"; traits <- c(X, Y, E)
setwd(file.path(T, "derived/munged"))
meta <- function(t) jsonlite::fromJSON(file.path(STD, paste0(t, ".meta.json")))
prev <- sapply(traits, function(t) { m <- meta(t); if (isTRUE(m$binary)) m$sample_prev else NA }); sp <- ifelse(is.na(prev), NA, 0.5)
m1 <- sprintf("%s =~ NA*%s + %s\n%s =~ NA*%s\n%s ~~ 1*%s\n%s ~~ 1*%s\n%s ~~ 0*%s\n%s ~~ 0*%s\n%s ~~ 0*%s\n",
              F1, X, Y, F2, Y, F1, F1, F2, F2, F1, F2, X, X, Y, Y)
m2 <- paste0(m1, E, " ~~ ", F1, " + ", F2, "\n")
pick <- function(est, a, b, op, col) as.numeric(est[[col]][est$op == op & ((est$lhs == a & est$rhs == b) | (est$lhs == b & est$rhs == a))])
rows <- list()
for (half in c("ODD", "EVEN")) {
  cov <- ldsc(sep_weights = TRUE, traits = paste0(traits, ".sumstats.gz"), sample.prev = sp, population.prev = prev, ld = LD, wld = WLD,
              trait.names = traits, select = half, ldsc.log = file.path(T, paste0("results/16_split_", half)))
  std <- ldsc(sep_weights = TRUE, traits = paste0(traits, ".sumstats.gz"), sample.prev = sp, population.prev = prev, ld = LD, wld = WLD,
              trait.names = traits, select = half, stand = TRUE, ldsc.log = file.path(T, paste0("results/16_split_std_", half)))
  e1 <- usermodel(cov, estimation = "DWLS", model = m1, std.lv = FALSE, imp_cov = FALSE)$results
  e2 <- usermodel(cov, estimation = "DWLS", model = m2, std.lv = FALSE, imp_cov = FALSE)$results
  s_non <- pick(e1, F2, Y, "=~", "STD_Genotype"); s_non_se <- pick(e1, F2, Y, "=~", "STD_Genotype_SE")
  r_non <- pick(e2, F2, E, "~~", "STD_Genotype"); r_non_se <- pick(e2, F2, E, "~~", "STD_Genotype_SE")
  k <- 3; idx <- matrix(0, k, k); idx[lower.tri(idx, diag = TRUE)] <- 1:6; se_of <- function(i, j) sqrt(diag(std$V_Stand))[idx[max(i, j), min(i, j)]]
  rg_YE <- std$S_Stand[2, 3]; rg_YE_se <- se_of(2, 3); rg_XE <- std$S_Stand[1, 3]; rg_XE_se <- se_of(1, 3); rg_XY <- std$S_Stand[1, 2]
  share <- s_non^2; share_se <- 2 * abs(s_non) * s_non_se
  add <- function(metric, est, se) rows[[length(rows) + 1]] <<- data.frame(topic = topic, half = half, metric = metric, estimate = est, se = se,
                                                                            lo = est - 1.96 * se, hi = est + 1.96 * se)
  add(sprintf("rg(%s, %s)", Y, X), rg_XY, se_of(1, 2))
  add(sprintf("share of %s genetic variance not via %s (std loading^2)", Y, X), share, share_se)
  add(sprintf("rg(%s, %s)", F2, E), r_non, r_non_se)
  add(sprintf("rg(%s, %s)", Y, E), rg_YE, rg_YE_se)
  add(sprintf("rg(%s, %s)", X, E), rg_XE, rg_XE_se)
  add("P2 = 1 - rg(F_NonBMI,T2D)/rg(OSA,T2D) (delta method, independence approx.)", 1 - r_non / rg_YE,
                         sqrt((r_non_se / rg_YE)^2 + (r_non * rg_YE_se / rg_YE^2)^2))
}
out <- do.call(rbind, rows); out$env <- ENV_VERSIONS
write.table(out, file.path(T, "results/16_split_half.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(out[, 2:7], digits = 3)
