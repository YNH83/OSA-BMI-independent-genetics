# Shared R environment (sourced at the top of every analysis script; run scripts from the repository root).
# 1. Project-local library (lavaan pinned to 0.6-19: lavaan 0.7.x rejects DWLS on continuous data as GenomicSEM 0.0.5 calls it).
# 2. Some sandboxed environments block the sysctl call behind parallel::detectCores(), which then returns NA and breaks lavaan.
#    Fix the core count explicitly; this changes no statistics, only lets packages start.
.libPaths(c(normalizePath("tools/genetics_env/R_library", mustWork = FALSE), .libPaths()))
local({
  ns <- asNamespace("parallel")
  unlockBinding("detectCores", ns)
  assign("detectCores", function(all.tests = FALSE, logical = TRUE) 8L, envir = ns)
  lockBinding("detectCores", ns)
})
options(mc.cores = 8)
suppressMessages({ library(lavaan); library(GenomicSEM) })
ENV_VERSIONS <- sprintf("R %s; GenomicSEM %s; lavaan %s", getRversion(), packageVersion("GenomicSEM"), packageVersion("lavaan"))
