# Colocalisation of the residual-component loci (C9 addendum 04 S3), preregistered 2026-09-27 17:11.
# coloc.abf, default priors p1 = p2 = 1e-4, p12 = 1e-5; window lead +/- 500 kb; SNPs common to both GWAS (HapMap3); betas aligned to
# the allele of the first dataset. Continuous traits: type "quant", sdY = 1 (F_NonBMI is a standardized factor; GIANT BMI is a
# standardized trait). Binary traits: type "cc", s = sample case proportion, N = total cases + controls.
# Prior sensitivity (added 2026-09-27, not preregistered): calls with p12 = 5e-6 and 1e-4.
# Decision per pair: PP.H4 >= 0.80 shared causal variant; PP.H3 >= 0.80 distinct; otherwise inconclusive.
# Usage (from the repository root): Rscript shared/coloc_loci.R   Output: results/18_coloc.tsv
source("env/setup_env.R")
suppressMessages({ library(coloc); library(data.table) })
STD <- normalizePath("data/hm3_std", mustWork = FALSE)
pos <- fread(normalizePath("tools/genetics_env/ref/reference.1000G.EUR.maf.txt", mustWork = FALSE), select = c("SNP", "CHR", "BP"))
meta <- function(t) jsonlite::fromJSON(file.path(STD, paste0(t, ".meta.json")))
T <- normalizePath(".")
L <- fread(file.path(T, "figures/Table4_loci.tsv"))[, .(lead_SNP, CHR, BP = BP_hg19, gene = nearest_gene)]
SRC <- list(F_NonBMI = file.path(T, "derived/gwas/F_NonBMI_QC.tsv.gz"), BMI_GIANT = file.path(STD, "BMI_GIANT.tsv.gz"),
            SNORING_CAMPOS = file.path(STD, "SNORING_CAMPOS.tsv.gz"), OSA_FG = file.path(STD, "OSA_FG.tsv.gz"))
PAIRS <- list(c("F_NonBMI", "BMI_GIANT"), c("F_NonBMI", "SNORING_CAMPOS"), c("OSA_FG", "BMI_GIANT"))
D <- lapply(SRC, function(f) merge(fread(f, select = c("SNP", "A1", "A2", "BETA", "SE", "N", "MAF")), pos, by = "SNP"))
ds <- function(name, x) {
  base <- list(beta = x$BETA, varbeta = x$SE^2, snp = x$SNP, MAF = pmin(pmax(x$MAF, 0.001), 0.5))
  if (!grepl("^F_", name) && isTRUE(meta(name)$binary)) {
    m <- meta(name); cs <- as.numeric(sub(".*from ([0-9]+) cases / ([0-9]+) controls.*", "\\1", m$n_mode)); ct <- as.numeric(sub(".*/ ([0-9]+) controls.*", "\\1", m$n_mode))
    c(base, list(type = "cc", s = cs / (cs + ct), N = cs + ct))
  } else c(base, list(type = "quant", sdY = 1, N = median(x$N, na.rm = TRUE)))
}
rows <- list()
for (i in seq_len(nrow(L))) {
  l <- L[i]; lo <- l$BP - 5e5; hi <- l$BP + 5e5
  for (p in PAIRS) {
    a <- D[[p[1]]][CHR == l$CHR & BP >= lo & BP <= hi]; b <- D[[p[2]]][CHR == l$CHR & BP >= lo & BP <= hi]
    m <- merge(a, b, by = "SNP", suffixes = c(".a", ".b"))
    same <- m$A1.a == m$A1.b & m$A2.a == m$A2.b; flip <- m$A1.a == m$A2.b & m$A2.a == m$A1.b
    m <- m[same | flip]; m[flip[same | flip], BETA.b := -BETA.b]
    xa <- m[, .(SNP, BETA = BETA.a, SE = SE.a, N = N.a, MAF = MAF.a)]; xb <- m[, .(SNP, BETA = BETA.b, SE = SE.b, N = N.b, MAF = MAF.a)]
    r <- suppressWarnings(suppressMessages(coloc.abf(ds(p[1], xa), ds(p[2], xb))))
    # prior sensitivity (not preregistered; GWAS/PRS guideline section 5.4): p12 = 5e-6 and 1e-4
    sens <- sapply(c(5e-6, 1e-4), function(q) { z <- suppressWarnings(suppressMessages(coloc.abf(ds(p[1], xa), ds(p[2], xb), p12 = q)))$summary
      if (z["PP.H4.abf"] >= 0.8) "shared" else if (z["PP.H3.abf"] >= 0.8) "distinct" else "inconclusive" })
    s <- r$summary
    call <- if (s["PP.H4.abf"] >= 0.8) "shared causal variant" else if (s["PP.H3.abf"] >= 0.8) "distinct causal variants" else "inconclusive"
    top <- r$results[which.max(r$results$SNP.PP.H4), "snp"]
    rows[[length(rows) + 1]] <- data.frame(locus = l$gene, lead_SNP = l$lead_SNP, trait1 = p[1], trait2 = p[2], nsnps = as.integer(s["nsnps"]),
      PP.H0 = s["PP.H0.abf"], PP.H1 = s["PP.H1.abf"], PP.H2 = s["PP.H2.abf"], PP.H3 = s["PP.H3.abf"], PP.H4 = s["PP.H4.abf"],
      call = call, top_SNP_H4 = top, call_p12_5e6 = sens[1], call_p12_1e4 = sens[2])
  }
}
out <- do.call(rbind, rows); out$env <- ENV_VERSIONS; out$coloc <- as.character(packageVersion("coloc"))
write.table(out, file.path(T, "results/18_coloc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(out[, c("locus", "trait1", "trait2", "nsnps", "PP.H3", "PP.H4", "call")], digits = 3, row.names = FALSE)
