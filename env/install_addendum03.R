# Additional GitHub packages for the MR and LAVA analyses (MR-PRESSO, LAVA, MRlap, CAUSE); run after install_R_packages.R
dir.create("tools/genetics_env/R_library", recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath("tools/genetics_env/R_library", mustWork = FALSE), .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = .libPaths()[1])
for (p in c("MRPRESSO", "LAVA", "MRlap", "cause")) cat(p, requireNamespace(p, quietly = TRUE), "\n")
res <- list()
for (r in c("rondolab/MR-PRESSO", "josefin-werme/LAVA", "n-mounier/MRlap", "jean997/cause")) {
  res[[r]] <- tryCatch({ remotes::install_github(r, lib = .libPaths()[1], upgrade = "never", dependencies = TRUE, build_vignettes = FALSE); "OK" },
                       error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("==>", r, res[[r]], "\n")
}
for (p in c("MRPRESSO", "LAVA", "MRlap", "cause")) cat("LOAD", p, requireNamespace(p, quietly = TRUE), "\n")
