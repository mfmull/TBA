# =============================================================================
# 48_run_all.R -- ONE-COMMAND PIPELINE (irrigation / land-use analysis)
#
#   Rscript 48_run_all.R
#
# Runs, in order: 43 (main + frozen benchmark), 44 (SI + robustness),
# 45 (summaries + Data_S2), 46 (figures), 47 (tables); then writes the output
# manifest and session_info.txt and runs the final verification assertions.
# Spec-level ensembles are cached under derived/cache/ with content stamps, so
# a re-run recomputes only what changed. First full run: roughly 3-4 h on
# 2 cores; cached re-run: minutes.
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
source(file.path(.root, "41_config.R"))
source(file.path(CFG$root, "42_core.R"))
log_open(CFG$log_file)
T0 <- Sys.time()

# ---- options ----------------------------------------------------------------
# FORCE rebuilds every cached ensemble (the 3-4 h path). Set it either as an
# environment variable or as a plain variable in the calling session; a session
# variable wins, so both of these work:
#
#   FORCE <- TRUE;         source("48_run_all.R")     # from RStudio
#   Sys.setenv(FORCE = 1); source("48_run_all.R")
#
# Because a session variable wins, a FORCE left over from an earlier run in the
# same RStudio session applies again -- the resolved value is echoed below so
# this stays visible. Clear it with rm(FORCE).
.opt <- function(name) {
  if (exists(name, envir = globalenv(), inherits = FALSE))
    return(isTRUE(as.logical(get(name, envir = globalenv()))))
  v <- Sys.getenv(name)
  nzchar(v) && !identical(tolower(v), "false") && !identical(v, "0")
}
FORCE <- .opt("FORCE")
sayf("Options: FORCE=%s\n", FORCE)
if (FORCE) say("*** CACHES WILL BE REBUILT -- expect roughly 3-4 hours. ***\n")

# 49 is the economic-asymmetry supplement. It reads the same data and uses the
# same fit_one() machinery, and now writes to the canonical output/, figure/
# and derived/cache/ rather than its own parallel copy of that layout.
for (s in c("43_run_main.R", "44_run_robustness.R", "45_summaries.R",
            "46_figures.R", "47_tables.R", "49_asymmetry.R")) {
  say("\n############ ", s, " ############\n")
  source(file.path(CFG$root, s))
}

# =============================================================================
# OUTPUT MANIFEST AND FINAL VERIFICATION
# =============================================================================
say("\n############ manifest and verification ############\n")
cfg <- CFG

canonical_out <- c(
  "table_main_results.csv",   "table_main_results.tex",
  "table_si_robustness.csv",  "table_si_robustness.tex",
  "table_si_ensemble.csv",    "table_si_ensemble.tex",
  "table_si_balance.csv",     "table_si_balance.tex",
  "table_si_sample_flow.csv", "table_si_sample_flow.tex",
  "table_si_summary_stats.csv", "table_si_summary_stats.tex",
  "methods_numbers.csv",      "methods_numbers.tex",
  "Dataset_S3.csv",           # SI supplementary dataset, segment level

  # economic-asymmetry supplement (49)
  "country_econ.csv",
  "asymmetry_dyads_agva_ha.csv",  "asymmetry_cells_agva_ha.csv",
  "asymmetry_dyads_gdp_pc.csv",   "asymmetry_cells_gdp_pc.csv",
  "asymmetry_quadrant_tests_agva_ha.csv",
  "asymmetry_quadrant_tests_gdp_pc.csv",
  "output_manifest.csv",      "session_info.txt")
canonical_fig <- c(
  "figure_main_irrigation.pdf",
  "figure_si_robustness.pdf",
  "figure_si_river_exclusion.pdf",
  "figure_si_balance.pdf",
  "figure_asymmetry_agva_ha.pdf",
  "figure_asymmetry_gdp_pc.pdf")

writeLines(utils::capture.output(utils::sessionInfo()),
           file.path(cfg$out_dir, "session_info.txt"))

manifest <- tibble::tibble(
  file = c(canonical_out, file.path("..", "figure", canonical_fig)),
  role = c(rep("table/data", length(canonical_out)),
           rep("figure", length(canonical_fig))))
manifest$md5 <- vapply(manifest$file, function(f) {
  p <- file.path(cfg$out_dir, f)
  if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_
}, character(1))
readr::write_csv(manifest, file.path(cfg$out_dir, "output_manifest.csv"))

# ---- assertions -------------------------------------------------------------
main <- readRDS(file.path(cfg$cache_dir, "main_objects.rds"))
rob  <- readRDS(file.path(cfg$cache_dir, "robustness_objects.rds"))
summ <- readRDS(file.path(cfg$cache_dir, "summaries.rds"))

# 1. every canonical file exists AND has real content. Existence alone is not
#    enough: a PDF device that fails writes nothing while reporting success,
#    which is exactly how the figures in this folder went stale while the
#    tables beside them kept updating. The right test differs by kind --
#      figures  a valid PDF is never under ~1 kB, so a size floor works
#      csv      check STRUCTURE: a header plus at least one data row, because
#               a correct small table can be under 200 bytes
bad <- character(0)
for (k in seq_along(c(canonical_out, canonical_fig))) {
  f <- c(canonical_out, canonical_fig)[k]
  p <- file.path(if (k <= length(canonical_out)) cfg$out_dir else cfg$fig_dir, f)
  why <- if (!file.exists(p)) "missing"
    else if (file.size(p) == 0) "empty"
    else if (grepl("[.]pdf$", f)) {
      if (file.size(p) < 1000) sprintf("only %d bytes", file.size(p)) else NA_character_
    } else if (grepl("[.]csv$", f)) {
      if (length(readLines(p, warn = FALSE)) < 2) "no data rows" else NA_character_
    } else NA_character_
  if (!is.na(why)) bad <- c(bad, sprintf("%s (%s)", f, why))
}
if (length(bad))
  stop("Canonical outputs not written correctly:\n  ",
       paste(bad, collapse = "\n  "))

# 2. figures are PDF-only in figure/
stopifnot(!any(grepl("[.](png|tif|tiff|jpe?g|svg)$",
                     list.files(cfg$fig_dir), ignore.case = TRUE)))

# 3. plotted values equal exported table values (main results)
tmr <- readr::read_csv(file.path(cfg$out_dir, "table_main_results.csv"),
                       show_col_types = FALSE)
stopifnot(abs(tmr$est_lmer[tmr$label == "IrIntens"] -
                main$IrIntens$best$treat_eff) < 1e-12)

# 4. best member selected on the STATED rule (lowest max |SMD|)
for (nm in names(main)) {
  sdf <- main[[nm]]$summary_df
  stopifnot(main[[nm]]$best_idx ==
              sdf$idx[order(sdf$max_smd, sdf$mean_smd)][1])
}

# 5. SMD statistics exclude the propensity distance row
stopifnot(isTRUE(cfg$smd_exclude_distance))

# 6. GW-subset count matches SI Dataset S3 (the 207 -> 159 fix)
stopifnot(summ$n_gw_treated == cfg$reference_n_treat_gw)

# 6b. SI Dataset S3 is exactly the declared analysis-variable manifest, one row
#     per treated segment, and still reproduces both sample filters. The column
#     assertion is the guard against silently re-admitting a descriptive column
#     (or silently losing one a new specification needs).
ds3 <- readr::read_csv(file.path(cfg$out_dir, cfg$si_dataset_name),
                       show_col_types = FALSE)
stopifnot(identical(names(ds3), cfg$si_dataset_cols),
          nrow(ds3) == cfg$reference_n_si_dataset_rows,
          nrow(dplyr::distinct(ds3, aq_id, CountryName)) == nrow(ds3),
          sum(ds3$Irrig   > 0, na.rm = TRUE) == cfg$reference_n_treat_ir,
          sum(ds3$GWIrrig > 0, na.rm = TRUE) == cfg$reference_n_treat_gw)

# 7. ensemble sizes complete
for (r in c(main, rob))
  stopifnot(nrow(r$summary_df) >= 0.95 * cfg$n_realizations)

say("\nAll verification assertions passed.\n")
sayf("TOTAL runtime: %.1f min.\n",
     as.numeric(difftime(Sys.time(), T0, units = "mins")))
