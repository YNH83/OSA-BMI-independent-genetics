# C9 null control 1 (preregistered, Stage 1 section 7): gout minus BMI, same GWAS-by-subtraction pipeline as OSA
# Model after Demange et al. 2021: SNP predicts both F_BMI and F_NonBMI.
# Usage: Rscript 03b_userGWAS_null_gout.R [NSNP=all] [CORES=1] [CHUNK=1] [NCHUNK=1]
#   NSNP > 0 runs the first NSNP SNPs only (runtime stress test before the main run).
#   GenomicSEM's own parallel mode needs local sockets (blocked in the original sandboxed environment); instead the SNP set is
#   split into NCHUNK contiguous chunks run as independent single-core processes, then merged (shared/merge_gwas.py).
source("env/setup_env.R")
a <- commandArgs(trailingOnly = TRUE)
NSNP  <- if (length(a) >= 1 && a[1] != "all") as.integer(a[1]) else NA
CORES <- if (length(a) >= 2) as.integer(a[2]) else 1
CHUNK <- if (length(a) >= 3) as.integer(a[3]) else 1
NCHUNK <- if (length(a) >= 4) as.integer(a[4]) else 1

TOPIC <- normalizePath(".")
STD   <- normalizePath("data/hm3_std", mustWork = FALSE)
REF   <- normalizePath("tools/genetics_env/ref/reference.1000G.EUR.maf.txt", mustWork = FALSE)
dir.create(file.path(TOPIC, "derived/gwas"), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(TOPIC, "derived/gwas"))

cov <- readRDS(file.path(TOPIC, "results/02_ldsc_T2D_XUE__GOUT_FG__BMI_GIANT.rds"))
subset_cov <- function(cv, keep) {
  # GenomicSEM stores S with colnames only; V follows vech order (lower triangle, column-major)
  nm <- colnames(cv$S); i <- match(keep, nm); stopifnot(!anyNA(i))
  k <- length(nm); pos <- matrix(0, k, k); pos[lower.tri(pos, diag = TRUE)] <- seq_len(k * (k + 1) / 2)
  sel <- c(); for (b in seq_along(i)) for (a in b:length(i)) sel <- c(sel, pos[max(i[a], i[b]), min(i[a], i[b])])
  S <- cv$S[i, i, drop = FALSE]; colnames(S) <- keep
  I <- cv$I[i, i, drop = FALSE]
  # element ORDER matters: userGWAS reads covstruc[[1]] as V and covstruc[[2]] as S (positional)
  list(V = cv$V[sel, sel, drop = FALSE], S = S, I = I, N = cv$N[, sel, drop = FALSE], m = cv$m)
}
keep <- c("BMI_GIANT", "GOUT_FG")
cov <- subset_cov(cov, keep)
g_rg <- read.delim(file.path(TOPIC, "results/02_numeric_T2D_XUE__GOUT_FG__BMI_GIANT.tsv"))$s_bmi
stopifnot(abs(cov2cor(cov$S)[1, 2] - g_rg) < 2e-3)   # guard: must reproduce rg(gout, BMI) from the current 02 null run

ss_file <- "sumstats_BMI_GOUT.rds"
if (!file.exists(ss_file)) {
  ss <- sumstats(files = file.path(STD, paste0(keep, ".tsv.gz")), ref = REF, trait.names = keep,
                 se.logit = c(FALSE, TRUE), OLS = c(TRUE, FALSE), linprob = c(FALSE, FALSE),
                 N = c(NA, NA), maf.filter = 0.01, keep.indel = FALSE, parallel = FALSE)
  # exclude MHC (hg19 chr6 25 to 35 Mb)
  ss <- ss[!(ss$CHR == 6 & ss$BP >= 25e6 & ss$BP <= 35e6), ]
  saveRDS(ss, ss_file)
} else ss <- readRDS(ss_file)
if (!is.na(NSNP)) ss <- ss[seq_len(min(NSNP, nrow(ss))), ]
if (NCHUNK > 1) { grp <- cut(seq_len(nrow(ss)), NCHUNK, labels = FALSE); ss <- ss[grp == CHUNK, ] }

model <- "
F_BMI =~ NA*BMI_GIANT + start(0.4)*GOUT_FG
F_NonBMI =~ NA*GOUT_FG
F_NonBMI ~~ 1*F_NonBMI
F_BMI ~~ 1*F_BMI
F_BMI ~~ 0*F_NonBMI
GOUT_FG ~~ 0*BMI_GIANT
GOUT_FG ~~ 0*GOUT_FG
BMI_GIANT ~~ 0*BMI_GIANT
F_BMI ~ SNP
F_NonBMI ~ SNP
"
t0 <- Sys.time()
res <- userGWAS(covstruc = cov, SNPs = ss, estimation = "DWLS", model = model, printwarn = TRUE,
                sub = c("F_BMI~SNP", "F_NonBMI~SNP"), cores = CORES, toler = 1e-50, SNPSE = FALSE,
                parallel = CORES > 1, GC = "standard", MPI = FALSE, smooth_check = TRUE)
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
tag <- paste0(if (is.na(NSNP)) "all" else paste0("n", NSNP), if (NCHUNK > 1) sprintf("_chunk%02dof%02d", CHUNK, NCHUNK) else "")
for (i in seq_along(res)) {
  nm <- c("F_BMI", "F_NonBMI")[i]
  f <- file.path(TOPIC, paste0("derived/gwas/null_gout_userGWAS_", nm, "_", tag, ".tsv"))
  data.table::fwrite(res[[i]], f, sep = "\t")   # this data.table build lacks zlib; gzip afterwards
  system2("gzip", c("-f", shQuote(f)))
}
msg <- sprintf("%s  SNPs=%d cores=%d elapsed=%.1f s (%.4f s/SNP)  %s", format(Sys.time()), nrow(ss), CORES, el, el / nrow(ss), ENV_VERSIONS)
cat(msg, "\n", file = file.path(TOPIC, "results/03b_null_gout_runtime.log"), append = TRUE); cat(msg, "\n")
