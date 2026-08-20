# =============================================================================
# 58_run_all.R -- ONE-COMMAND PIPELINE (dyad typology)
#
#   Rscript 58_run_all.R                      # or source() it from RStudio
#
# Runs, in order: 53 (preferred classification + frozen benchmark),
# 54 (robustness variants), 55 (summaries), 56 (figures), 57 (tables); then
# writes the output manifest and session_info.txt and runs the final
# verification assertions.
#
# The whole pipeline is a deterministic classification over a few hundred rows:
# a full run takes seconds, and nothing here needs parallelism. Caching exists
# so figures can be re-rendered without recomputing, and is content-stamped.
#
# -----------------------------------------------------------------------------
# OPTIONS. Set either as an environment variable or as a plain variable in the
# calling session; a session variable wins, so both of these work:
#
#   Sys.setenv(FORCE = 1); source("58_run_all.R")
#   FORCE <- TRUE;         source("58_run_all.R")
#
#   FORCE   rebuild the cached robustness variants rather than reusing them.
#
# Because a session variable wins, a FORCE left over from an earlier run in the
# same RStudio session applies again -- the resolved value is echoed below.
# Clear it with rm(FORCE).
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
source(file.path(.root, "51_config.R"))
source(file.path(CFG$root, "52_core.R"))
log_open(CFG$log_file)
T0 <- Sys.time()

.opt <- function(name) {
  if (exists(name, envir = globalenv(), inherits = FALSE))
    return(isTRUE(as.logical(get(name, envir = globalenv()))))
  v <- Sys.getenv(name)
  nzchar(v) && !identical(tolower(v), "false") && !identical(v, "0")
}
FORCE <- .opt("FORCE")
sayf("Options: FORCE=%s\n", FORCE)

# 59 runs last and reports only: it diffs the published classification against
# the dyad-specific border correction and writes its own tables. It skips
# itself when the GEE exports are absent, and it never alters the preferred
# specification -- so its outputs are deliberately NOT in the canonical
# manifest below, which asserts the published deliverables.
for (s in c("53_run_main.R", "54_run_robustness.R", "55_summaries.R",
            "56_figures.R", "57_tables.R", "59_dyad_correction.R")) {
  say("\n############ ", s, " ############\n")
  source(file.path(CFG$root, s))
}

# =============================================================================
# OUTPUT MANIFEST AND FINAL VERIFICATION
# =============================================================================
say("\n############ manifest and verification ############\n")
cfg <- CFG

canonical_out <- c(
  "Dataset_S4.csv",                  # SI supplementary dataset, dyad level
  "table_main_typology.csv",         "table_main_typology.tex",
  "table_si_group_stats.csv",        "table_si_group_stats.tex",
  "table_si_contingency.csv",        "table_si_contingency.tex",
  "table_si_contingency_tests.csv",  "table_si_contingency_tests.tex",
  "table_si_robustness.csv",         "table_si_robustness.tex",
  "table_si_bl_profile.csv",         "table_si_bl_profile.tex",
  "table_si_sample_flow.csv",        "table_si_sample_flow.tex",
  "methods_numbers.csv",             "methods_numbers.tex",
  "output_manifest.csv",             "session_info.txt")
canonical_fig <- c("figure_main_typology.pdf",
                   "figure_si_typology_bars.pdf",
                   "figure_si_landuse_scale.pdf",
                   "figure_si_robustness.pdf")
if (have_vcd)
  canonical_fig <- c(canonical_fig,
                     paste0("figure_si_mosaic_",
                            names(readRDS(file.path(cfg$cache_dir,
                              "main_objects.rds"))$mosaics), ".pdf"))

writeLines(utils::capture.output(utils::sessionInfo()),
           file.path(cfg$out_dir, "session_info.txt"))

manifest <- tibble::tibble(
  file = c(canonical_out, canonical_fig),
  role = c(rep("table/data", length(canonical_out)),
           rep("figure", length(canonical_fig))),
  dir  = c(rep("output", length(canonical_out)),
           rep("figure", length(canonical_fig))))
manifest$md5 <- vapply(seq_len(nrow(manifest)), function(i) {
  p <- file.path(if (manifest$dir[i] == "output") cfg$out_dir else cfg$fig_dir,
                 manifest$file[i])
  if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_
}, character(1))
readr::write_csv(manifest, file.path(cfg$out_dir, "output_manifest.csv"))

# ---- assertions --------------------------------------------------------------
main <- readRDS(file.path(cfg$cache_dir, "main_objects.rds"))
summ <- readRDS(file.path(cfg$cache_dir, "summaries.rds"))
robj <- readRDS(file.path(cfg$cache_dir, "robustness_objects.rds"))

# 1. every canonical file exists AND has real content. Existence alone is not
#    enough -- a PDF device that fails writes nothing while reporting success --
#    but the right test differs by kind, and a single byte floor gets it wrong:
#    a correct 5-row table is under 200 bytes.
#      figures  a valid PDF is never under ~1 kB, so a floor works
#      csv/tex  check STRUCTURE instead: a header plus at least one data row
bad <- character(0)
for (i in seq_len(nrow(manifest))) {
  f <- manifest$file[i]
  p <- file.path(if (manifest$dir[i] == "output") cfg$out_dir else cfg$fig_dir, f)
  why <- if (!file.exists(p)) "missing"
    else if (file.size(p) == 0) "empty"
    else if (grepl("[.]pdf$", f)) {
      if (file.size(p) < 1000) sprintf("only %d bytes", file.size(p)) else NA_character_
    } else if (grepl("[.]csv$", f)) {
      nl <- length(readLines(p, warn = FALSE))
      if (nl < 2) "no data rows" else NA_character_
    } else NA_character_
  if (!is.na(why)) bad <- c(bad, sprintf("%s (%s)", f, why))
}
if (length(bad))
  stop("Canonical outputs not written correctly:\n  ",
       paste(bad, collapse = "\n  "))

# 2. figures are PDF only
stopifnot(!any(grepl("[.](png|tif|tiff|jpe?g|svg)$", list.files(cfg$fig_dir),
                     ignore.case = TRUE)))

# 3. the exported table matches the cached classification
tmt <- readr::read_csv(file.path(cfg$out_dir, "table_main_typology.csv"),
                       show_col_types = FALSE)
stopifnot(identical(as.integer(tmt$n_current),
                    as.integer(main$counts$current[cfg$class_levels])),
          identical(as.integer(tmt$n_future),
                    as.integer(main$counts$future[cfg$class_levels])))

# 4. Dataset S4 covers every dyad exactly once, with no unclassified row, is
#    exactly the declared column manifest, and its stage-3 columns were folded
#    into class/classF rather than shipped twice.
ds4 <- readr::read_csv(file.path(cfg$out_dir, cfg$si_dataset_name),
                       show_col_types = FALSE)
stopifnot(nrow(ds4) == cfg$reference_n_dyads,
          nrow(dplyr::distinct(ds4, code, CC_1, CC_2)) == nrow(ds4),
          !any(is.na(ds4$class)), !any(is.na(ds4$classF)),
          identical(names(ds4), cfg$si_dataset_cols),
          !any(c("class_3", "classF_3") %in% names(ds4)))
# the stage ledger must still agree with the classification the run produced
stopifnot(identical(as.integer(table(ds4$class)[cfg$class_levels]),
                    as.integer(cfg$reference_class[cfg$class_levels])))

# 5. the frozen benchmark (53 already stops otherwise; re-assert here)
stopifnot(identical(as.integer(main$counts$current[cfg$class_levels]),
                    as.integer(cfg$reference_class[cfg$class_levels])),
          identical(as.integer(main$counts$future[cfg$class_levels]),
                    as.integer(cfg$reference_classF[cfg$class_levels])))

# 6. every robustness variant reproduces its frozen counts
for (nm in names(robj$variants)) {
  ref <- cfg$reference_variants[[nm]]
  if (is.null(ref)) next
  stopifnot(identical(robj$variants[[nm]]$counts$current, as.integer(ref$class)),
            identical(robj$variants[[nm]]$counts$future,  as.integer(ref$classF)))
}

# 7. the claim the SI makes about robustness: the three-way partition is
#    stable even where the BL-DA split is not. Asserted, not asserted-in-prose.
p <- robj$partition
stopifnot(max(p$urban) - min(p$urban) <= 60,
          max(p$interior) - min(p$interior) <= 60)

# 8. Does BL keep its defining profile -- highest overdraft exposure -- in
#    every variant? This is a scientific finding, not a code invariant, so it
#    is REPORTED rather than asserted: a variant in which it stopped holding
#    would be a result worth reading, not a reason to halt the pipeline.
if (all(robj$bl_profile$BL_is_max)) {
  say("BL has the highest overdraft exposure in every robustness variant.\n")
} else {
  say("\n*** NOTE: BL is NOT the highest-overdraft class in: ",
      paste(robj$bl_profile$variant[!robj$bl_profile$BL_is_max],
            collapse = ", "),
      "\n    The robustness paragraph in the SI should say so explicitly.\n\n")
}

say("\nAll verification assertions passed.\n")
sayf("TOTAL runtime: %.2f min.\n",
     as.numeric(difftime(Sys.time(), T0, units = "mins")))
