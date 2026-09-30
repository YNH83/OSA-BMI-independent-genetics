# Single-trait LDSC (h2 and intercept) for the post-QC multivariable GWAS of the residual component F_NonBMI, for the QQ figure.
# N = N_hat (Mallard); observed scale. Writes results/15_ldsc_F_NonBMI_ldsc.log and results/15_ldsc_factor_h2.tsv.
# Usage (from the repository root): Rscript shared/ldsc_factor_h2.R
source("env/setup_env.R")
HM3 <- normalizePath("ref/ldsc_ref/w_hm3.snplist", mustWork = FALSE)
LD  <- paste0(normalizePath("tools/genetics_env/ldsc/eur_ld", mustWork = FALSE), "/")
WLD <- paste0(normalizePath("tools/genetics_env/ldsc/eur_wld_noMHC25_35", mustWork = FALSE), "/")
for (x in list(c(".", "F_NonBMI"))) {
  T <- normalizePath(x[1]); f <- x[2]
  dir.create(file.path(T, "derived/munged"), showWarnings = FALSE); setwd(file.path(T, "derived/munged"))
  munge(files = file.path(T, "derived/gwas", paste0(f, "_QC.tsv.gz")), hm3 = HM3, trait.names = f, info.filter = 0.9, maf.filter = 0.01)
  r <- ldsc(traits = paste0(f, ".sumstats.gz"), sample.prev = NA, population.prev = NA, ld = LD, wld = WLD, trait.names = f,
            sep_weights = TRUE, ldsc.log = file.path(T, paste0("results/15_ldsc_", f)))
  out <- data.frame(trait = f, h2_observed = r$S[1, 1], h2_SE = sqrt(r$V[1, 1]), intercept = r$I[1, 1])
  write.table(out, file.path(T, "results/15_ldsc_factor_h2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); print(out)
}
