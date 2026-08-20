# =============================================================================
# 44_run_robustness.R -- SI AND ROBUSTNESS SPECIFICATIONS
#
# Block A (SI, Fig. S8 counterparts): crop suitability, irrigation need,
#   overdraft conditional on irrigated area, stringent (>6 mo) scarcity,
#   double-robust (regression-augmented) intensity and Gini.
# Block B (attribution sensitivity; reviewer R5.5): groundwater-irrigation
#   Gini under the permissive GW-fraction > 0 attribution (GWShareGt0).
# Block C (river-border exclusion; reviewer R5.4): main intensity and Gini
#   contrasts re-estimated after dropping every polygon (treated or control)
#   within 50 km of a river-defined international border, complementing the
#   regression-adjusted "| Rivers" specifications.
#
# Writes derived/cache/robustness_objects.rds.
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!exists("CFG")) source(file.path(.root, "41_config.R"))
if (!exists("run_spec")) source(file.path(CFG$root, "42_core.R"))
if (is.null(.log$path)) log_open(CFG$log_file)
cfg <- CFG
t0 <- Sys.time()

si_specs  <- names(cfg$specs)[vapply(cfg$specs, function(s) s$group == "si",
                                     logical(1))]
rob_specs <- names(cfg$specs)[vapply(cfg$specs, function(s) s$group == "robust",
                                     logical(1))]

say("Block A+B (SI): ", paste(si_specs, collapse = ", "), "\n")
if (!exists("FORCE")) FORCE <- nzchar(Sys.getenv("FORCE"))
rob <- run_specs(si_specs, cfg, force = isTRUE(FORCE))

say("\nBlock C (river-border exclusion): ", paste(rob_specs, collapse = ", "), "\n")
rob <- c(rob, run_specs(rob_specs, cfg, force = isTRUE(FORCE)))

for (nm in names(rob)) stopifnot(is.finite(rob[[nm]]$best$treat_eff))

saveRDS(rob, file.path(cfg$cache_dir, "robustness_objects.rds"))
sayf("\n44_run_robustness total: %.1f min\n",
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
