# =============================================================================
# assemble_SI_datasets.R -- BUILD THE SUBMITTED SUPPLEMENTARY DATASETS
#
# The SI_datasets/ folder used to be assembled by hand, which is how the
# numbering drifted from the pipelines that produce it. This script is the
# single place where a pipeline output becomes a numbered SI dataset. Run it
# from the repository root AFTER the three pipelines have run:
#
#   Rscript 2_wells/38_run_all.R
#   Rscript 3_irrigation/48_run_all.R
#   Rscript 4_typology/58_run_all.R
#   Rscript assemble_SI_datasets.R
#
# NUMBERING (2026-08-19). Datasets are numbered in SI section order --
# wells, then irrigation, then typology -- which is also the order of the
# \dataset{} blocks at the end of tex/_SI_vFINAL.tex. Because \dataset uses a
# positional LaTeX counter, the numbers follow the order of those blocks and
# every \ref{data:...} in the MS and SI resolves automatically.
#
#   S1  wells segments (44)          2_wells/output/table_si_aquifer_stability.csv
#   S2  IGRAC inventory (990)        0_preprocess_igrac/igrac_segment_inventory.csv
#   S3  irrigation segments (929)    3_irrigation/output/Dataset_S3.csv          [NEW]
#   S4  dyad classification (568)    4_typology/output/Dataset_S4.csv       [MERGED]
#   S5  audit census (568)           4_typology/output/Data_S3_classification_audit.csv
#
# WHAT CHANGED relative to the 2026-08-14 assembly:
#   * the old S3 (dyad classification) and the old S5 (stage-by-stage build)
#     are now ONE dataset -- the new S4. They were 568 rows on the same
#     (code, CC_1, CC_2) key, so the merge is 1:1 and lossless. class_3 and
#     classF_3 are not carried: they are identical to class and classF.
#   * a NEW segment-level dataset (S3) publishes the irrigation analysis
#     sample -- every outcome, sample filter and covariate the 41-series
#     specifications use, and nothing else.
#   * the audit census moves from S4 to S5. Its contents are unchanged.
# =============================================================================

suppressPackageStartupMessages({library(readr); library(dplyr)})

ROOT <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})()
OUT <- file.path(ROOT, "supplementary_datasets")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

say <- function(...) cat(paste0(...))

# name, source path, expected rows, key columns -------------------------------
SPEC <- list(
  list(n = "Dataset_S1.csv",
       src = file.path(ROOT, "2_wells", "output", "table_si_aquifer_stability.csv"),
       rows = 44L,  key = c("segment")),
  list(n = "Dataset_S2.csv",
       src = file.path(ROOT, "0_preprocess_igrac", "igrac_segment_inventory.csv"),
       rows = 990L, key = c("aq_id", "CC")),
  list(n = "Dataset_S3.csv",
       src = file.path(ROOT, "3_irrigation", "output", "Dataset_S3.csv"),
       rows = 929L, key = c("aq_id", "CountryName")),
  list(n = "Dataset_S4.csv",
       src = file.path(ROOT, "4_typology", "output", "Dataset_S4.csv"),
       rows = 568L, key = c("code", "CC_1", "CC_2")),
  list(n = "Dataset_S5.csv",
       src = file.path(ROOT, "4_typology", "output",
                       "Data_S3_classification_audit.csv"),
       rows = 568L, key = c("aquifer_code", "country_1", "country_2")))

say("\n############ assembling SI_datasets/ ############\n")
for (s in SPEC) {
  if (!file.exists(s$src))
    stop("missing pipeline output for ", s$n, ": ", s$src,
         "\n  Run the producing pipeline first (see the header of this file).")
  d <- readr::read_csv(s$src, show_col_types = FALSE, progress = FALSE)

  # A supplementary dataset that silently changes length is the failure mode
  # this script exists to prevent, so length and key uniqueness are asserted,
  # not reported.
  if (nrow(d) != s$rows)
    stop(s$n, ": expected ", s$rows, " rows, got ", nrow(d))
  if (!all(s$key %in% names(d)))
    stop(s$n, ": missing key column(s) ",
         paste(setdiff(s$key, names(d)), collapse = ", "))
  if (any(duplicated(d[, s$key])))
    stop(s$n, ": key ", paste(s$key, collapse = " x "), " is not unique")

  file.copy(s$src, file.path(OUT, s$n), overwrite = TRUE)
  say(sprintf("  %-16s %4d rows x %2d cols   <- %s\n", s$n, nrow(d), ncol(d),
              sub(paste0("^", ROOT, "/"), "", s$src)))
}

# The two datasets that used to be separate must still be recoverable from the
# merged one; check the columns of both ancestors survived.
ds4 <- readr::read_csv(file.path(OUT, "Dataset_S4.csv"), show_col_types = FALSE)
stopifnot(all(c("class", "classF", "cooperation_level", "has_river_border",
                "dist_class") %in% names(ds4)),          # ex-Dataset S3
          all(c("class_stage1", "class_stage2", "measured_stage2",
                "changed_overall", "audit_error_class", "audit_tier",
                "near_share_stage1", "near_share_stage2",
                "comment_stage2", "comment_stage3") %in% names(ds4)))

manifest <- bind_rows(lapply(SPEC, function(s) {
  d <- readr::read_csv(file.path(OUT, s$n), show_col_types = FALSE)
  tibble(dataset = s$n, rows = nrow(d), cols = ncol(d),
         md5 = unname(tools::md5sum(file.path(OUT, s$n))))
}))
readr::write_csv(manifest, file.path(OUT, "MANIFEST.csv"))
say("\nWrote SI_datasets/MANIFEST.csv\n")
print(as.data.frame(manifest), row.names = FALSE)
say("\nAll assembly assertions passed.\n")
