# C9 Stage 2 addendum 03, A6 (preregistered 2026-09-27 02:32): LAVA local genetic correlation.
# Loci: LAVA blocks_s2500_m25_f1_w200 (GRCh37), MHC-overlapping blocks removed. LD reference: 1000G EUR (g1000_eur).
# Sample overlap: LDSC cross-trait intercept matrix (GenomicSEM ldsc, same settings as the main analysis), converted to a
# correlation matrix as in the LAVA sample-overlap vignette. Univariate filter: local h2 P < 1e-4 for every trait in a test.
# A6a/A6c: bivariate local rg of F_NonBMI (and OSA for reference) with each trait. A6b: partial local rg(OSA, X | BMI).
# Usage: Rscript 11_lava_A6.R prep            -> writes LAVA inputs and the overlap matrix
#        Rscript 11_lava_A6.R run CHUNK NCHUNK -> analyses a contiguous chunk of loci (single-thread processes)
source("env/setup_env.R")
suppressMessages(library(LAVA))
a <- commandArgs(trailingOnly = TRUE); MODE <- a[1]
TOPIC <- normalizePath(".")
STD   <- normalizePath("data/hm3_std", mustWork = FALSE)
HM3   <- normalizePath("ref/ldsc_ref/w_hm3.snplist", mustWork = FALSE)
BFILE <- file.path(normalizePath("tools/genetics_env/magma/g1000_eur", mustWork = FALSE), "g1000_eur")
LOC   <- normalizePath("tools/genetics_env/src_pkgs/LAVA-main/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile", mustWork = FALSE)
D     <- file.path(TOPIC, "derived/lava"); dir.create(D, showWarnings = FALSE)
meta  <- function(t) jsonlite::fromJSON(file.path(STD, paste0(t, ".meta.json")))
EXT   <- c("T2D_XUE", "HTN_FG", "HTN_PANUKB", "SNORING_CAMPOS", "INSOMNIA_FG", "DEPRESSION_FG", "FEV1FVC_SHRINE", "SLEEPDUR_DASHTI", "CHRONOTYPE_JONES")
BASE  <- c("OSA_FG", "BMI_GIANT")
if (MODE == "prep") {
  traits <- c(BASE, EXT[file.exists(file.path(STD, paste0(EXT, ".tsv.gz")))])
  info <- data.frame(phenotype = character(), cases = numeric(), controls = numeric(), filename = character())
  for (t in traits) {
    x <- data.table::fread(file.path(STD, paste0(t, ".tsv.gz")))
    y <- data.frame(SNP = x$SNP, A1 = x$A1, A2 = x$A2, N = x$N, Z = x$BETA / x$SE)
    f <- file.path(D, paste0(t, ".lava.tsv")); data.table::fwrite(y, f, sep = "\t")
    m <- meta(t); ca <- if (isTRUE(m$binary)) round(m$sample_prev * 1e6) else NA
    info <- rbind(info, data.frame(phenotype = t, cases = if (isTRUE(m$binary)) m$sample_prev else NA,
                                   controls = if (isTRUE(m$binary)) 1 - m$sample_prev else NA, filename = f))
  }
  # F_NonBMI (A9 QC output, N = N_hat)
  x <- data.table::fread(file.path(TOPIC, "derived/gwas/F_NonBMI_QC.tsv.gz"))
  f <- file.path(D, "NONBMI.lava.tsv"); data.table::fwrite(data.frame(SNP = x$SNP, A1 = x$A1, A2 = x$A2, N = x$N, Z = x$BETA / x$SE), f, sep = "\t")
  info <- rbind(info, data.frame(phenotype = "NONBMI", cases = NA, controls = NA, filename = f))
  info$filename <- basename(info$filename)   # LAVA reads the info file with read.table: the project path contains spaces
  write.table(info, file.path(D, "input.info.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  # sample overlap from LDSC intercepts (munge LAVA-format files with GenomicSEM)
  setwd(D)
  # munge needs P: use the standardized files (and the A9 QC file for F_NonBMI) for the overlap LDSC, not the Z-only LAVA files
  srcs <- ifelse(info$phenotype == "NONBMI", file.path(TOPIC, "derived/gwas/F_NonBMI_QC.tsv.gz"), file.path(STD, paste0(info$phenotype, ".tsv.gz")))
  munge(files = srcs, hm3 = HM3, trait.names = info$phenotype, info.filter = 0.9, maf.filter = 0.01)
  cov <- ldsc(traits = paste0(info$phenotype, ".sumstats.gz"), sample.prev = rep(NA, nrow(info)), population.prev = rep(NA, nrow(info)),
              ld = paste0(normalizePath("tools/genetics_env/ldsc/eur_ld", mustWork = FALSE), "/"), wld = paste0(normalizePath("tools/genetics_env/ldsc/eur_wld_noMHC25_35", mustWork = FALSE), "/"),
              sep_weights = TRUE, trait.names = info$phenotype, ldsc.log = file.path(D, "lava_overlap"))
  I <- cov$I; dimnames(I) <- list(info$phenotype, info$phenotype)
  O <- round(cov2cor(I), 5)                                # LAVA vignette: standardised intercept matrix
  write.table(O, file.path(D, "sample_overlap.txt"), quote = FALSE)
  cat("prep done:", nrow(info), "traits\n"); print(O)
} else {
  CH <- as.integer(a[2]); NCH <- as.integer(a[3]); setwd(D)
  info <- read.delim(file.path(D, "input.info.txt"))
  input <- process.input(input.info.file = file.path(D, "input.info.txt"), sample.overlap.file = file.path(D, "sample_overlap.txt"),
                         ref.prefix = BFILE, phenos = info$phenotype)
  loci <- read.loci(LOC)
  loci <- loci[!(loci$CHR == 6 & loci$STOP >= 25e6 & loci$START <= 35e6), ]
  loci <- loci[cut(seq_len(nrow(loci)), NCH, labels = FALSE) == CH, ]
  targets <- intersect(c("NONBMI", "OSA_FG"), info$phenotype); others <- intersect(EXT, info$phenotype)
  U <- list(); B <- list(); PC <- list()
  for (i in seq_len(nrow(loci))) {
    lc <- tryCatch(process.locus(loci[i, ], input), error = function(e) NULL)
    if (is.null(lc)) next
    un <- tryCatch(run.univ(lc), error = function(e) NULL); if (is.null(un)) next
    un$locus <- loci$LOC[i]; U[[i]] <- un
    ok <- un$phen[un$p < 1e-4]
    for (tg in intersect(targets, ok)) {
      oth <- intersect(others, ok)
      if (length(oth)) {
        bv <- tryCatch(run.bivar(lc, target = tg, phenos = oth), error = function(e) NULL)
        if (!is.null(bv)) { bv$locus <- loci$LOC[i]; B[[length(B) + 1]] <- bv }
      }
    }
    if (all(c("OSA_FG", "BMI_GIANT") %in% ok)) for (x in intersect(c("HTN_FG", "HTN_PANUKB", "T2D_XUE"), ok)) {
      pc <- tryCatch(run.pcor(lc, target = c("OSA_FG", x), phenos = c("OSA_FG", x, "BMI_GIANT")), error = function(e) NULL)
      if (!is.null(pc)) { pc$locus <- loci$LOC[i]; PC[[length(PC) + 1]] <- pc }
    }
  }
  w <- function(L, n) if (length(L)) data.table::fwrite(do.call(rbind, L), file.path(D, sprintf("%s_chunk%02dof%02d.tsv", n, CH, NCH)), sep = "\t")
  w(U, "univ"); w(B, "bivar"); w(PC, "pcor")
  cat("chunk", CH, "done:", nrow(loci), "loci\n")
}
