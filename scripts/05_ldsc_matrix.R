# C9 full LDSC genetic correlation matrix (Supplementary Tables 1 and 2; Figure 2a).
# Uses whatever of the listed traits have been standardized; binary traits use Neff with sample.prev 0.5
# and population prevalence = sample prevalence (preregistered convention).
source("env/setup_env.R")
TOPIC <- normalizePath(".")
STD   <- normalizePath("data/hm3_std", mustWork = FALSE)
HM3   <- normalizePath("ref/ldsc_ref/w_hm3.snplist", mustWork = FALSE)
LD    <- paste0(normalizePath("tools/genetics_env/ldsc/eur_ld", mustWork = FALSE), "/")
WLD   <- paste0(normalizePath("tools/genetics_env/ldsc/eur_wld_noMHC25_35", mustWork = FALSE), "/")
setwd(file.path(TOPIC, "derived/munged"))
want <- c("OSA_FG", "OSA_FG_INCLAVO", "OSA_CAMPOS", "OSA_CAMPOS_ADJBMI", "OSA_STRAUSZ_BMILT30", "BMI_GIANT", "WHRADJBMI_PULIT",
          "T2D_XUE", "T2D_FG", "HTN_FG", "CAD_ARAGAM", "GOUT_FG")
traits <- want[file.exists(file.path(STD, paste0(want, ".tsv.gz")))]
cat("traits available:", paste(traits, collapse = ", "), "\n")
meta <- function(t) jsonlite::fromJSON(file.path(STD, paste0(t, ".meta.json")))
for (t in traits) if (!file.exists(paste0(t, ".sumstats.gz")))
  munge(files = file.path(STD, paste0(t, ".tsv.gz")), hm3 = HM3, trait.names = t, info.filter = 0.9, maf.filter = 0.01)
prev <- sapply(traits, function(t) { m <- meta(t); if (isTRUE(m$binary)) m$sample_prev else NA })
sp <- ifelse(is.na(prev), NA, 0.5)
cov <- ldsc(sep_weights = TRUE, traits = paste0(traits, ".sumstats.gz"), sample.prev = sp, population.prev = prev, ld = LD, wld = WLD,
            trait.names = traits, ldsc.log = file.path(TOPIC, "results/05_ldsc_matrix"))
std <- ldsc(sep_weights = TRUE, traits = paste0(traits, ".sumstats.gz"), sample.prev = sp, population.prev = prev, ld = LD, wld = WLD,
            trait.names = traits, stand = TRUE, ldsc.log = file.path(TOPIC, "results/05_ldsc_matrix_std"))
k <- length(traits); idx <- matrix(0, k, k); idx[lower.tri(idx, diag = TRUE)] <- seq_len(k * (k + 1) / 2)
se_of <- function(V, i, j) sqrt(diag(V))[idx[max(i, j), min(i, j)]]
h2 <- data.frame(trait = traits, liability_or_observed_h2 = diag(cov$S),
                 h2_SE = sapply(seq_len(k), function(i) se_of(cov$V, i, i)), intercept = diag(cov$I),
                 scale = ifelse(is.na(prev), "observed (continuous)", "liability"), population_prev = prev)
rg <- do.call(rbind, lapply(seq_len(k), function(i) do.call(rbind, lapply(seq_len(k), function(j) if (j < i) {
  r <- std$S_Stand[i, j]; s <- se_of(std$V_Stand, i, j)
  data.frame(trait1 = traits[j], trait2 = traits[i], rg = r, SE = s, P = 2 * pnorm(-abs(r / s)))
}))))
write.table(h2, file.path(TOPIC, "results/05_h2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rg, file.path(TOPIC, "results/05_rg.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(cov = cov, std = std), file.path(TOPIC, "results/05_ldsc_matrix.rds"))
print(h2); print(rg[grepl("OSA_FG$", rg$trait1) | grepl("OSA_FG$", rg$trait2), ])
