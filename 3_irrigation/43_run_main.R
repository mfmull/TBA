# =============================================================================
# 43_run_main.R -- MAIN SPECIFICATIONS (Fig. 3 of the manuscript)
#
# Runs the seven preferred specifications across the 500-realization control
# ensemble and caches one result object per specification under
# derived/cache/spec_<label>.rds. Then verifies the frozen benchmark (if
# locked) and writes derived/cache/main_objects.rds for 45-47.
#
# Estimand: ATT under full optimal matching; outcome model: weighted linear
# mixed model with country random intercept (v3-identical). Nothing in this
# script changes a model relative to the legacy 4_Analyse.R -- only balance
# accounting (SMD scope, best-member rule) and diagnostics differ; those do
# not enter the ensemble model fits.
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!exists("CFG")) source(file.path(.root, "41_config.R"))
if (!exists("run_spec")) source(file.path(CFG$root, "42_core.R"))
if (is.null(.log$path)) log_open(CFG$log_file)
cfg <- CFG
t0 <- Sys.time()

main_specs <- names(cfg$specs)[vapply(cfg$specs, function(s)
  s$group %in% c("main_intensity", "main_gini"), logical(1))]
say("Main specifications: ", paste(main_specs, collapse = ", "), "\n\n")

# FORCE is resolved by 48_run_all.R before this file is sourced; the default
# here is only for running 43 on its own.
if (!exists("FORCE")) FORCE <- nzchar(Sys.getenv("FORCE"))
if (isTRUE(FORCE)) say("*** FORCE: rebuilding cached ensembles. ***\n")
main <- run_specs(main_specs, cfg, force = isTRUE(FORCE))

# ---- design sanity assertions ----------------------------------------------
stopifnot(main$IrIntens$n_treat == cfg$reference_n_treat_ir)
stopifnot(main$GWGini$n_treat  == cfg$reference_n_treat_gw_model)
for (nm in main_specs) {
  stopifnot(nrow(main[[nm]]$summary_df) >= 0.95 * cfg$n_realizations)
  stopifnot(is.finite(main[[nm]]$best$treat_eff))
}

# ---- frozen benchmark -------------------------------------------------------
chk <- function(ref, val, lab) {
  if (!is.finite(ref)) {
    if (!isTRUE(cfg$allow_unlocked_reference))
      stop("Reference for ", lab, " is unlocked and allow_unlocked_reference is FALSE.")
    sayf("  [benchmark] %s unlocked; observed %.10g (lock this in 41_config.R)\n",
         lab, val)
  } else if (abs(val - ref) > cfg$reproduce_tol * max(1, abs(ref))) {
    stop(sprintf("Frozen benchmark FAILED for %s: observed %.10g, locked %.10g",
                 lab, val, ref))
  } else sayf("  [benchmark] %s reproduced (%.10g)\n", lab, val)
}
chk(cfg$reference_IrIntens_meaneff,  main$IrIntens$counts$mean_eff,  "IrIntens mean effect")
chk(cfg$reference_IrIntens_sharesig, main$IrIntens$counts$share_sig, "IrIntens share significant")
chk(cfg$reference_IrGini_meaneff,    main$IrGini$counts$mean_eff,    "IrGini mean effect")
chk(cfg$reference_IrGini_sharesig,   main$IrGini$counts$share_sig,   "IrGini share significant")

saveRDS(main, file.path(cfg$cache_dir, "main_objects.rds"))
sayf("\n43_run_main total: %.1f min\n",
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
