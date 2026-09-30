# C9 deviation record 03 (P-value underflow fix): run userGWAS only for SNPs recovered by the P-underflow fix, then merge into the main GWAS output.
# STATUS 2026-09-27: already applied once. The pre-fix sumstats_BMI_OSA.rds it compares against was deleted afterwards
# (only the latest outputs were kept), so this script is kept as the record of the repair and is not meant to be rerun.
# Valid only if the LDSC covariance structure is unchanged by the fix (checked by 01_killswitch before running this):
# sumstats() standardizes each SNP from its own BETA/SE/MAF, so existing SNP results are unaffected; the model, covstruc,
# estimator, and options are identical to 03_userGWAS_subtraction.R.
source("env/setup_env.R")
TOPIC <- normalizePath(".")
STD   <- normalizePath("data/hm3_std", mustWork = FALSE)
REF   <- normalizePath("tools/genetics_env/ref/reference.1000G.EUR.maf.txt", mustWork = FALSE)
setwd(file.path(TOPIC, "derived/gwas"))
src <- readLines(file.path(TOPIC, "scripts/03_userGWAS_subtraction.R"))
# reuse subset_cov, keep, model and guard exactly as in the main script
eval(parse(text = src[grep("^cov <- readRDS", src):grep("^stopifnot\\(abs", src)]))
eval(parse(text = src[grep("^model <- \"", src):(grep("^model <- \"", src) + 11)]))
old <- readRDS("sumstats_BMI_OSA.rds")
new <- sumstats(files = file.path(STD, paste0(keep, ".tsv.gz")), ref = REF, trait.names = keep,
                se.logit = c(FALSE, TRUE), OLS = c(TRUE, FALSE), linprob = c(FALSE, FALSE),
                N = c(NA, NA), maf.filter = 0.01, keep.indel = FALSE, parallel = FALSE)
new <- new[!(new$CHR == 6 & new$BP >= 25e6 & new$BP <= 35e6), ]
saveRDS(new, "sumstats_BMI_OSA_fixP.rds")
# existing SNPs must be identical between the two sumstats builds
common <- intersect(old$SNP, new$SNP)
o <- old[match(common, old$SNP), ]; n <- new[match(common, new$SNP), ]
num <- sapply(o, is.numeric); dif <- max(abs(as.matrix(o[, num]) - as.matrix(n[, num])), na.rm = TRUE)
cat("common SNPs", length(common), "; max abs difference in sumstats columns", dif, "\n"); stopifnot(dif < 1e-8)
add <- new[!new$SNP %in% old$SNP, ]
cat("recovered SNPs to run:", nrow(add), "\n")
res <- userGWAS(covstruc = cov, SNPs = add, estimation = "DWLS", model = model, printwarn = TRUE,
                sub = c("F_BMI~SNP", "F_NonBMI~SNP"), cores = 1, toler = 1e-50, SNPSE = FALSE,
                parallel = FALSE, GC = "standard", MPI = FALSE, smooth_check = TRUE)
for (i in seq_along(res)) {
  nm <- c("F_BMI", "F_NonBMI")[i]
  f <- file.path(TOPIC, paste0("derived/gwas/userGWAS_", nm, "_recovered_fixP.tsv"))
  data.table::fwrite(res[[i]], f, sep = "\t"); system2("gzip", c("-f", shQuote(f)))
}
cat("done", ENV_VERSIONS, "\n")
