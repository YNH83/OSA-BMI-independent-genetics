# C9 kill-switch (preregistered: analysis_plan/C9_預先登錄_Stage1.md, section 4)
# Retire if OSA liability SNP-h2 z < 5, or rg(OSA, BMI) > 0.85.
# Binary trait: munged with Neff, so sample.prev = 0.5; population prevalence = FinnGen sample prevalence (preregistered).
source("env/setup_env.R")

TOPIC <- normalizePath(".")
STD   <- normalizePath("data/hm3_std", mustWork = FALSE)
HM3   <- normalizePath("ref/ldsc_ref/w_hm3.snplist", mustWork = FALSE)
LD    <- paste0(normalizePath("tools/genetics_env/ldsc/eur_ld", mustWork = FALSE), "/")
WLD   <- paste0(normalizePath("tools/genetics_env/ldsc/eur_wld_noMHC25_35", mustWork = FALSE), "/")
dir.create(file.path(TOPIC, "derived/munged"), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(TOPIC, "derived/munged"))

traits <- c("OSA_FG", "BMI_GIANT")
for (t in traits) {
  if (!file.exists(paste0(t, ".sumstats.gz")))
    munge(files = file.path(STD, paste0(t, ".tsv.gz")), hm3 = HM3, trait.names = t, info.filter = 0.9, maf.filter = 0.01)
}

osa_prev <- 56885 / (56885 + 441137)
res <- ldsc(sep_weights = TRUE, traits = paste0(traits, ".sumstats.gz"),
            sample.prev = c(0.5, NA), population.prev = c(osa_prev, NA),
            ld = LD, wld = WLD, trait.names = traits,
            ldsc.log = file.path(TOPIC, "results/01_killswitch"))
saveRDS(res, file.path(TOPIC, "results/01_killswitch_ldsc.rds"))

k <- nrow(res$S)
SE <- matrix(0, k, k); SE[lower.tri(SE, diag = TRUE)] <- sqrt(diag(res$V)); SE[upper.tri(SE)] <- t(SE)[upper.tri(SE)]
h2_osa <- res$S[1, 1]; h2_osa_se <- SE[1, 1]
rg <- cov2cor(res$S)[1, 2]
# delta-method-free SE for rg via standardized ldsc run
res_std <- ldsc(sep_weights = TRUE, traits = paste0(traits, ".sumstats.gz"), sample.prev = c(0.5, NA), population.prev = c(osa_prev, NA),
                ld = LD, wld = WLD, trait.names = traits, stand = TRUE, ldsc.log = file.path(TOPIC, "results/01_killswitch_std"))
rg_se <- sqrt(diag(res_std$V_Stand))[2]

z <- h2_osa / h2_osa_se
verdict <- if (z < 5) "RETIRE: OSA h2 z < 5" else if (rg > 0.85) "RETIRE: rg(OSA,BMI) > 0.85" else "PASS: proceed to Cholesky"
out <- data.frame(metric = c("OSA liability h2", "OSA h2 SE", "OSA h2 z", "rg(OSA,BMI)", "rg SE", "OSA LDSC intercept", "BMI h2 (observed)", "verdict"),
                  value = c(round(h2_osa, 4), round(h2_osa_se, 4), round(z, 2), round(rg, 4), round(rg_se, 4),
                            round(res$I[1, 1], 4), round(res$S[2, 2], 4), verdict))
write.table(out, file.path(TOPIC, "results/01_killswitch_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(out)
