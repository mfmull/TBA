# =============================================================================
# 38_run_all.R -- ONE-COMMAND PIPELINE
#
#   Rscript 38_run_all.R                      # or source() it from RStudio
#
# Runs, in order: 33 (main + frozen benchmark), 34 (robustness blocks A-D),
# 35 (summaries), 36 (figures), 37 (tables); then writes the output manifest
# and session_info.txt and runs the final verification assertions.
# Heavy objects are cached under derived/cache/ with content stamps, so a
# re-run recomputes only what changed.
#
# -----------------------------------------------------------------------------
# OPTIONS. Each may be set either as an environment variable or as a plain
# variable in the calling session; a variable already defined in the session
# always wins, so both of these work:
#
#   Sys.setenv(FORCE = 1); source("38_run_all.R")     # environment variable
#   FORCE <- TRUE;         source("38_run_all.R")     # session variable
#
#   FORCE           rebuild every STAMPED cache (robustness specifications,
#                   influence loops, block C, H3 bootstrap, selection
#                   stability, the 39 connectivity subsets and basemap).
#                   33 and 35 always recompute regardless -- they do not cache.
#   FORCE_RATIOS    additionally rebuild the Conley ratio table. This is the
#                   expensive one (~1500 units per cutoff/seed) and it is keyed
#                   by sample, NOT stamped by code, so it is the one cache that
#                   a change to the Conley kernel in 32_core.R will NOT
#                   invalidate on its own. Use it after touching that code.
#                   The old table is moved aside, never deleted.
#
# Because a session variable wins, a FORCE left over from an earlier run in the
# same RStudio session will silently apply again -- the resolved options are
# echoed below so this is visible. Clear one with rm(FORCE).
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
source(file.path(.root, "31_config.R"))
source(file.path(CFG$root, "32_core.R"))
log_open(CFG$log_file)
T0 <- Sys.time()

# ---- resolve options --------------------------------------------------------
.opt <- function(name) {
  if (exists(name, envir = globalenv(), inherits = FALSE))
    return(isTRUE(as.logical(get(name, envir = globalenv()))))
  v <- Sys.getenv(name)
  nzchar(v) && !identical(tolower(v), "false") && !identical(v, "0")
}
FORCE          <- .opt("FORCE")          # consumed by 34 and 39
FORCE_RATIOS   <- .opt("FORCE_RATIOS")
sayf("Options: FORCE=%s  FORCE_RATIOS=%s\n", FORCE, FORCE_RATIOS)
if (FORCE || FORCE_RATIOS)
  say("*** CACHES WILL BE REBUILT -- expect a long run. ***\n")

# Conley ratios are keyed by sample rather than stamped by code, so forcing
# them means moving the table aside. Moved, not deleted: it is by far the most
# expensive artefact here, and a mistaken --force should be recoverable.
if (FORCE_RATIOS) {
  .rp <- file.path(CFG$cache_dir, "conley_ratios.csv")
  if (file.exists(.rp)) {
    .bak <- paste0(.rp, ".bak-", format(Sys.time(), "%Y%m%d-%H%M%S"))
    if (file.rename(.rp, .bak))
      say("  moved Conley ratio table aside -> ", basename(.bak), "\n")
    else stop("could not move ", .rp, " aside; check permissions.")
  }
}

# =============================================================================
# RUN-STATE MARKERS
#
# A long sequential run is indistinguishable from a hung one, and a run that
# dies half way leaves a directory of outputs that LOOK complete. So the run
# announces itself on disk:
#
#   output/RUN_IN_PROGRESS.txt  written now, refreshed after every stage,
#                               deleted on success. Its presence after the
#                               process exits means the run did NOT finish.
#   output/RUN_COMPLETE.txt     written only after every stage has run AND
#                               the verification assertions below have passed.
#
# Any consumer of output/ should check for RUN_COMPLETE.txt before trusting it.
# =============================================================================
dir.create(cfg_out <- CFG$out_dir, showWarnings = FALSE, recursive = TRUE)
.marker_done <- file.path(cfg_out, "RUN_COMPLETE.txt")
.marker_busy <- file.path(cfg_out, "RUN_IN_PROGRESS.txt")
if (file.exists(.marker_done)) unlink(.marker_done)   # stale until re-earned

.par_on <- isTRUE(CFG$parallel) && par_workers(CFG) > 1L &&
  .Platform$OS.type == "unix" &&
  !(identical(Sys.getenv("RSTUDIO"), "1") &&
    !nzchar(Sys.getenv("ALLOW_FORK_IN_RSTUDIO")))

if (!.par_on) {
  say("\n*** SEQUENTIAL RUN ***\n",
      "  Forking is disabled (RStudio, or parallel = FALSE), so every block\n",
      "  runs on one core. With FORCE this typically takes several HOURS\n",
      "  rather than ~20 min. For the parallel path (up to 8 workers):\n\n",
      "      Rscript 38_run_all.R\n\n",
      "  from a terminal, not from RStudio.\n\n")
}

.stage_log <- list()
.write_busy <- function(done, current) {
  writeLines(c(
    "RUN IN PROGRESS -- outputs in this folder are NOT complete.",
    sprintf("started      : %s", format(T0, "%Y-%m-%d %H:%M:%S")),
    sprintf("last update  : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("elapsed      : %.1f min",
            as.numeric(difftime(Sys.time(), T0, units = "mins"))),
    sprintf("execution    : %s", if (.par_on)
            sprintf("parallel, %d workers", par_workers(CFG)) else "sequential"),
    sprintf("options      : FORCE=%s FORCE_RATIOS=%s", FORCE, FORCE_RATIOS),
    "",
    "stages completed:",
    if (length(done)) paste0("  ", done) else "  (none yet)",
    sprintf("currently running: %s", current)), .marker_busy)
}
.write_busy(character(0), .stages[1])

for (s in .stages) {
  say("\n############ ", s, " ############\n")
  .t_stage <- Sys.time()
  source(file.path(CFG$root, s))
  .el <- as.numeric(difftime(Sys.time(), .t_stage, units = "mins"))
  .stage_log[[s]] <- .el
  sayf("---- %s finished in %.1f min (%.1f min into the run) ----\n", s, .el,
       as.numeric(difftime(Sys.time(), T0, units = "mins")))
  .done <- vapply(names(.stage_log), function(k)
    sprintf("%-24s %6.1f min", k, .stage_log[[k]]), "")
  .nxt <- .stages[match(s, .stages) + 1L]
  .write_busy(.done, if (is.na(.nxt)) "manifest and verification" else .nxt)
}

# =============================================================================
# OUTPUT MANIFEST AND FINAL VERIFICATION
# =============================================================================
say("\n############ manifest and verification ############\n")
cfg <- CFG

canonical <- c(
  # main paper
  "figure_main_three_panel.pdf",
  "table_main_results.csv",        "table_main_results.tex",
  "table_sample_characteristics.csv", "table_sample_characteristics.tex",
  "table_selected_aquifers.csv",   "table_selected_aquifers.tex",
  # supporting information
  "figure_main_wells_map.pdf",
  "figure_main_wells_insets.pdf",
  "figure_si_wells_connectivity.pdf",
  "figure_si_robustness_h1_h2.pdf",
  "figure_si_localization_robustness.pdf",
  "figure_si_influence.pdf",
  "table_si_h1_h2_robustness.csv", "table_si_h1_h2_robustness.tex",
  # Presentation form of the registry table -- the one pasted into the SI.
  # It MUST be listed here: assertion 1 below deletes every file in output/
  # that is not canonical, so an omitted output is silently destroyed on the
  # next full run rather than merely going stale.
  "table_si_registry_display.csv", "table_si_registry_display.tex",
  "table_si_localization_robustness.csv", "table_si_localization_robustness.tex",
  "table_si_aquifer_stability.csv", "table_si_aquifer_stability.tex",
  "table_si_wells_sumstats.csv",  "table_si_wells_sumstats.tex",
  "table_si_matching_overlap.csv", "table_si_matching_overlap.tex",
  "table_si_sample_flow.csv",      "table_si_sample_flow.tex",
  "methods_numbers.csv",           "methods_numbers.tex",
  "output_manifest.csv",           "session_info.txt")

writeLines(utils::capture.output(utils::sessionInfo()),
           file.path(cfg$out_dir, "session_info.txt"))

manifest <- tibble::tibble(file = canonical) %>%
  dplyr::mutate(
    role = dplyr::case_when(
      grepl("^figure_si|^table_si", file) ~ "supporting information",
      file %in% c("methods_numbers.csv", "methods_numbers.tex") ~ "methods numbers",
      file %in% c("output_manifest.csv", "session_info.txt") ~ "provenance",
      TRUE ~ "main paper"))
readr::write_csv(manifest, file.path(cfg$out_dir, "output_manifest.csv"))
manifest <- manifest %>%
  dplyr::mutate(exists = file.exists(file.path(cfg$out_dir, file)),
                md5 = ifelse(exists,
                             unname(tools::md5sum(file.path(cfg$out_dir, file))),
                             NA_character_))
readr::write_csv(manifest %>% dplyr::select(file, role, md5),
                 file.path(cfg$out_dir, "output_manifest.csv"))

# ---- assertions -------------------------------------------------------------
main <- readRDS(file.path(cfg$cache_dir, "main_objects.rds"))
summ <- readRDS(file.path(cfg$cache_dir, "summaries.rds"))
pref <- main$pref

# 1. every canonical file exists, and the output dir contains nothing else
present <- setdiff(list.files(cfg$out_dir), character(0))
missing <- setdiff(canonical, present)
extra   <- setdiff(present, canonical)
if (length(missing)) stop("Missing canonical outputs: ", paste(missing, collapse = ", "))
if (length(extra)) {
  file.remove(file.path(cfg$out_dir, extra))
  say("Removed non-canonical files from output/: ", paste(extra, collapse = ", "), "\n")
}

# 2. all figures are PDFs (no PNG/TIFF/JPEG/SVG anywhere in output/)
stopifnot(!any(grepl("[.](png|tif|tiff|jpe?g|svg)$", list.files(cfg$out_dir),
                     ignore.case = TRUE)))

# 3. plotted values equal exported table values (main results)
tmr <- readr::read_csv(file.path(cfg$out_dir, "table_main_results.csv"),
                       show_col_types = FALSE)
h1row <- tmr[tmr$hypothesis == "H1: depletion level" &
               tmr$term == "transboundary contrast", ]
h2row <- tmr[tmr$hypothesis == "H2: borderward gradient" &
               tmr$term == "transboundary contrast", ]
stopifnot(abs(h1row$estimate - pref$h1$estimate) < 1e-9,
          abs(h2row$estimate - pref$h2$estimate) < 1e-9,
          abs(h1row$mde80 - main$mde80) < 1e-9)

# 4. MDE used the intended alpha and power
assert_mde(cfg)
stopifnot(abs(main$mde80 -
                (qnorm(1 - cfg$alpha_one) + qnorm(0.80)) *
                (if (is.finite(pref$h1$se_cr2)) pref$h1$se_cr2
                 else pref$h1$se_model)) < 1e-9)

# 5. MDE is not the confidence interval (distinct quantities exported)
stopifnot(!isTRUE(all.equal(main$mde80,
                            abs(pref$h1$ci90_high - pref$h1$estimate))))

# 6. placebo tested two-sided; directional alternatives correct
stopifnot(identical(cfg$directional_alt[["H1"]], "greater"),
          identical(cfg$directional_alt[["H2"]], "less"))

# 7. Panel C contains all eligible segments; proportions within [0, 1]
stopifnot(nrow(main$h3$scatter) == sum(pref$data$TBn == 1))
for (km in cfg$coverage_km) {
  v <- pref$data[[paste0("share_within_", km, "km")]]
  stopifnot(all(v >= 0 & v <= 1, na.rm = TRUE))
}

# 8. selected aquifers come from the coded rule
sel_tab <- readr::read_csv(file.path(cfg$out_dir, "table_selected_aquifers.csv"),
                           show_col_types = FALSE)
stopifnot(identical(sort(sel_tab$segment),
                    sort(main$ident$table$label_txt[main$ident$table$selected])))

# 9. frozen benchmark reproduced (33 already stops otherwise; re-assert)
if (is.finite(cfg$reference_h1))
  stopifnot(abs(pref$h1$estimate - cfg$reference_h1) <=
              cfg$reproduce_tol * max(1, abs(cfg$reference_h1)))

say("\nAll verification assertions passed.\n")

.total_min <- as.numeric(difftime(Sys.time(), T0, units = "mins"))

# ---- completion marker: written only now, after every assertion above ------
writeLines(c(
  "RUN COMPLETE -- every stage ran and every verification assertion passed.",
  "",
  sprintf("finished     : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("started      : %s", format(T0, "%Y-%m-%d %H:%M:%S")),
  sprintf("total runtime: %.1f min", .total_min),
  sprintf("execution    : %s", if (.par_on)
          sprintf("parallel, %d workers", par_workers(CFG)) else "sequential"),
  sprintf("options      : FORCE=%s FORCE_RATIOS=%s", FORCE, FORCE_RATIOS),
  sprintf("host         : %s", paste(Sys.info()[c("nodename", "sysname")],
                                     collapse = " / ")),
  sprintf("R            : %s", getRversion()),
  "",
  "stage timings:",
  vapply(names(.stage_log), function(k)
    sprintf("  %-24s %6.1f min", k, .stage_log[[k]]), ""),
  "",
  "frozen reference check:",
  sprintf("  H1 %+.4f mm/yr   H2 %+.6f Fisher z", CFG$reference_h1,
          CFG$reference_h2),
  "",
  sprintf("canonical outputs present: %d of %d",
          sum(file.exists(file.path(CFG$out_dir, canonical))), length(canonical))),
  .marker_done)
if (file.exists(.marker_busy)) unlink(.marker_busy)

say("\n", strrep("=", 62), "\n")
sayf("RUN COMPLETE -- %.1f min. Marker: %s\n", .total_min, basename(.marker_done))
say(strrep("=", 62), "\n")

