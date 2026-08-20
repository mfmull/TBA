# =============================================================================
# 59_dyad_correction.R -- CORRECTION STEP 2: DYAD-SPECIFIC BORDER ZONES
#
# Emits the class / classF correction for each dyad when the near (10 km) and
# far (100 km) zones are measured against that dyad's OWN border instead of
# against `lb`, the global land-border layer. This is step 2 of the four-stage
# schema:
#
#   1  segment-based zones, distance sensitivity      (near5 / far200 variants)
#   2  dyad-border correction                          <- THIS SCRIPT
#   3  external audit adjudication                     (AUDIT_T1 / AUDIT_T12)
#   4  inference on the finalised dataset
#
# The ledger is written as class_in -> class_out with a `step` column, so
# steps 2 and 3 stack in either order. See the note below on the four dyads
# where the order is not immaterial.
#
# Changes nothing: CFG$dyad_zones stays FALSE, 53_run_main.R still reproduces
# the frozen benchmark, and every published number stands. Skips itself
# cleanly when the GEE exports are absent.
#
# Writes:
#   output/table_dyad_border_ledger.csv    one row per dyad, in -> out
#   output/table_dyad_border_matrix.csv    correction matrices, long form
#   output/table_dyad_border_changed.csv   the reclassified dyads only
#   output/table_dyad_border_excess.csv    written by apply_dyad_zones()
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!exists("CFG")) source(file.path(.root, "51_config.R"))
if (!exists("classify_dyads")) source(file.path(CFG$root, "52_core.R"))
if (is.null(.log$path)) log_open(CFG$log_file)
cfg <- CFG
t0 <- Sys.time()

say("\n======== CORRECTION STEP 2: DYAD-SPECIFIC BORDER ZONES =========\n")

if (!dyad_zone_ready(cfg)) {
  say("Correction exports not present; nothing to compute.\n")
  for (s in names(cfg$dyad_zone_files))
    if (!length(dyad_zone_parts(cfg, s)))
      sayf("  no file matching %s  (zone %s)\n", cfg$dyad_zone_files[[s]], s)
  say("Run gee/dyad_border_corrections.js, then gee/finish_dyad_10k.js for\n",
      "any aquifer it could not finish, and drop the CSVs in 1_data/geeOut/.\n")
} else {

sayf("Corrected zones: %s   (near = %s, far = %s)\n",
     paste(names(cfg$dyad_zone_files), collapse = ", "),
     cfg$near_zone, cfg$far_zone)

res <- dyad_border_correction(cfg)
led <- res$ledger
lv  <- cfg$class_levels

readr::write_csv(led,        file.path(cfg$out_dir, "table_dyad_border_ledger.csv"))
readr::write_csv(res$matrix, file.path(cfg$out_dir, "table_dyad_border_matrix.csv"))

# ---- coverage ---------------------------------------------------------------
sayf("\n%d of %d dyads measured against their own border; %d keep the\n",
     sum(led$measured), nrow(led), sum(!led$measured))
say("segment-wide value because the export covers 3+-riparian aquifers only.\n")
say("  An unmeasured dyad is NOT a verified-unchanged dyad -- it is one that\n",
    "  was not measured. dyadBilateralRisk.csv is the evidence on those.\n")

# ---- counts ------------------------------------------------------------------
say("\n---- class counts (IU0, IU1, UB, DA, BL) ----\n")
sayf("%-14s %-24s %s\n", "", "current", "future")
sayf("%-14s %-24s %s\n", "in  (published)",
     paste(res$counts$current_in,  collapse = ", "),
     paste(res$counts$future_in,   collapse = ", "))
sayf("%-14s %-24s %s\n", "out (corrected)",
     paste(res$counts$current_out, collapse = ", "),
     paste(res$counts$future_out,  collapse = ", "))
sayf("%-14s %-24s %s\n", "delta",
     paste(sprintf("%+d", res$counts$current_out - res$counts$current_in),
           collapse = ", "),
     paste(sprintf("%+d", res$counts$future_out - res$counts$future_in),
           collapse = ", "))

# ---- correction matrices -----------------------------------------------------
for (scen in c("current", "future")) {
  sayf("\n---- correction matrix, %s (rows in, cols out) ----\n", scen)
  m <- res$matrix %>% filter(scenario == scen)
  print(stats::xtabs(n ~ published + corrected, m))
}

# The rule orders the classes from least to most bilateral interaction, and
# the dyad zone is a subset of the lb zone, so a corrected dyad can only move
# DOWN that order. Any upward move would mean the correction invented land
# use, and is worth stopping to look at rather than reporting.
ord <- setNames(seq_along(lv), lv)
up  <- led %>% filter(ord[class_out]  > ord[class_in] |
                      ord[classF_out] > ord[classF_in])
sayf("\nMonotonicity: %d dyad(s) move UP the IU0<IU1<UB<DA<BL order.\n", nrow(up))
if (nrow(up))
  say("  *** Unexpected. A smaller zone cannot add land use; inspect these.\n")
else
  say("  As expected -- every change is toward the interior end, because a\n",
      "  dyad-specific zone is a subset of the all-borders zone.\n")

# ---- what moved --------------------------------------------------------------
sayf("\n%d dyad(s) change class, %d change classF, %d change either.\n",
     sum(led$changed), sum(led$changedF), sum(led$changed | led$changedF))
sayf("Of the %d measured, %d (%.0f%%) keep both classes.\n", sum(led$measured),
     sum(led$measured & !led$changed & !led$changedF),
     100 * sum(led$measured & !led$changed & !led$changedF) /
       max(sum(led$measured), 1))

changed <- led %>% filter(changed | changedF) %>%
  arrange(desc(abs(near_share_in - near_share_out)))
readr::write_csv(changed,
                 file.path(cfg$out_dir, "table_dyad_border_changed.csv"))
if (nrow(changed)) {
  say("\n---- reclassified dyads (largest near-zone shift first) ----\n")
  print(as.data.frame(changed %>%
    select(code, CC_1, CC_2, class_in, class_out, classF_in, classF_out,
           near_share_in, near_share_out) %>% utils::head(50)),
    row.names = FALSE, digits = 3)
}

# ---- collision with step 3 ---------------------------------------------------
say("\n---- overlap with the external audit (step 3) ----\n")
both <- led %>% filter((changed | changedF) & audit_t1)
sayf("%d of the %d dyads this step moves are also AUDIT_T1 rows.\n",
     nrow(both), sum(led$changed | led$changedF))
if (nrow(both)) {
  print(as.data.frame(both %>%
    select(code, CC_1, CC_2, class_in, class_out, audit_from, audit_to,
           audit_from_stale)), row.names = FALSE)
  say("\nThese are the only rows where the order of steps 2 and 3 matters.\n",
      "  audit_from_stale = TRUE means the audit recorded a starting class\n",
      "  that this step has already changed, so applying step 3 afterwards\n",
      "  would rewrite a class the audit never actually examined. Where\n",
      "  audit_to == class_out the two agree independently, which is\n",
      "  corroboration rather than a conflict.\n")
}

# ---- freeze the variant ------------------------------------------------------
say("\n---- to freeze dyad_border in reference_variants (51_config.R) ----\n")
sayf("    dyad_border = list(class  = c(%s),\n",
     paste(sprintf("%dL", res$counts$current_out), collapse = ", "))
sayf("                       classF = c(%s)),\n",
     paste(sprintf("%dL", res$counts$future_out), collapse = ", "))

sayf("\n59_dyad_correction.R done in %.2f min.\n",
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
