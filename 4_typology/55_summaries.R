# =============================================================================
# 55_summaries.R -- METHODS NUMBERS, SI2 CLASSIFICATION TABLE, SAMPLE FLOW
#
# Every count quoted in the typology section of the manuscript is generated
# here, not typed. Writes derived/cache/summaries.rds; 57_tables.R exports.
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!exists("CFG")) source(file.path(.root, "51_config.R"))
if (!exists("classify_dyads")) source(file.path(CFG$root, "52_core.R"))
if (is.null(.log$path)) log_open(CFG$log_file)
cfg <- CFG
t0 <- Sys.time()

say("\n================ SUMMARIES =====================\n")
main <- readRDS(file.path(cfg$cache_dir, "main_objects.rds"))
rob  <- readRDS(file.path(cfg$cache_dir, "robustness_objects.rds"))
tp <- main$typology
cur <- main$counts$current; fut <- main$counts$future

# ---- SI Dataset S4: the consolidated dyad-level table ------------------------
# The manuscript cites this dataset for the final cooperation scores, so the
# dyad characteristics used in the mosaics travel with the classification
# rather than living only in the figure.
#
# 2026-08-19: the stage-by-stage build ledger (formerly shipped separately as
# Dataset S5) is merged in here, so the SI carries ONE dyad-level dataset. Key
# is (code, CC_1, CC_2); both sides are 568 rows on that key, so the join is
# 1:1 and nothing is dropped or duplicated -- asserted below rather than
# assumed.
si_class <- tp %>%
  mutate(cooperation_level = .data$ArrangementScore,
         dist_class = distance_class(load_dyads(cfg), cfg)) %>%
  select(code, name, CC_1, CC_2, class, classF,
         cooperation_level, has_river_border, dist_class)

stages <- readr::read_csv(cfg$stages_file, show_col_types = FALSE) %>%
  rename(class_stage1  = class_1,  class_stage2  = class_2,
         classF_stage1 = classF_1, classF_stage2 = classF_2,
         changed_stage2 = changed_2, changed_stage3 = changed_3,
         near_share_stage1 = near_share_1, near_share_stage2 = near_share_2)

key <- c("code", "CC_1", "CC_2")
stopifnot(
  nrow(si_class) == cfg$reference_n_dyads,
  nrow(stages)   == cfg$reference_n_dyads,
  !any(duplicated(si_class[key])), !any(duplicated(stages[key])),
  nrow(dplyr::anti_join(si_class, stages, by = key)) == 0,
  nrow(dplyr::anti_join(stages, si_class, by = key)) == 0)

si_dyads <- si_class %>%
  left_join(stages %>% select(-name), by = key) %>%
  arrange(code, CC_1, CC_2)

# Stage 3 IS the finalised classification. If the ledger and the pipeline ever
# disagree, the ledger is stale and the merge must not silently ship both.
# (as.character: `class` may arrive as a factor from classify_dyads() while the
# ledger is read as plain character, and identical() would fail on the type.)
stopifnot(identical(as.character(si_dyads$class),
                    as.character(si_dyads$class_3)),
          identical(as.character(si_dyads$classF),
                    as.character(si_dyads$classF_3)))
si_dyads <- si_dyads[, cfg$si_dataset_cols, drop = FALSE]
sayf("Dataset S4: %d dyads x %d columns (classification + stage ledger).\n",
     nrow(si_dyads), ncol(si_dyads))

si2 <- si_dyads   # legacy name, kept so older drivers keep running

# ---- methods numbers ---------------------------------------------------------
mn <- list()
put <- function(key, description, value, unit, source, where) {
  mn[[length(mn) + 1]] <<- tibble(key = key, description = description,
                                  value = value, unit = unit,
                                  source = source, location = where)
}
n <- main$counts$n_dyads
put("nDyads", "transboundary aquifer dyads", n, "dyads",
    "dyad table", "Methods, typology")
put("nCountries", "countries appearing in at least one dyad",
    main$counts$n_countries, "countries", "dyad table", "Methods, typology")
put("nAquifers", "aquifers with at least one dyad",
    main$counts$n_aquifers, "aquifers", "dyad table", "Methods, typology")
for (k in cfg$class_levels) {
  put(paste0("n", k), paste("dyads classified", k), cur[[k]], "dyads",
      "classification", "Results, typology")
  put(paste0("share", k), paste("share of dyads classified", k),
      round(100 * cur[[k]] / n, 1), "%", "classification", "Results, typology")
}
put("nInterior", "interior dyads (IU0 + IU1)", cur[["IU0"]] + cur[["IU1"]],
    "dyads", "classification", "Results, typology")
put("shareInterior", "share of dyads that are interior",
    round(100 * (cur[["IU0"]] + cur[["IU1"]]) / n, 1), "%",
    "classification", "Results, typology")
put("shareInteriorUrban", "share interior or urban-border (IU + UB)",
    round(100 * (cur[["IU0"]] + cur[["IU1"]] + cur[["UB"]]) / n, 1), "%",
    "classification", "Results, typology")
for (k in cfg$class_levels)
  put(paste0("nF", k), paste("dyads classified", k, "under the future scenario"),
      fut[[k]], "dyads", "future classification", "Results, typology")
put("blGrowth", "factor increase in BL under the future scenario",
    round(fut[["BL"]] / cur[["BL"]], 2), "x",
    "future classification", "Results, typology")
put("daGrowth", "factor increase in DA under the future scenario",
    round(fut[["DA"]] / cur[["DA"]], 2), "x",
    "future classification", "Results, typology")
# UNITS, verified empirically (see SI S5.6). The cropland bands (CR, IR, GW,
# CR3, GW3, IR3) carry 1e8 m2 = 10,000 ha per unit; the URBAN band carries
# 1e7 m2 = 1,000 ha per unit. They differ by ten because Crop1k_0307 is not a
# 0-1 fraction while the urban band is binary. So the same eps is 10 ha of
# groundwater-irrigated cropland and 1 ha of urban land, and both are reported.
put("epsHa", "detection threshold, groundwater-irrigated cropland",
    cfg$eps * 1e4, "ha",
    "51_config.R", "Methods, typology")
put("epsHaUrban", "detection threshold, urban land", cfg$eps * 1e3, "ha",
    "51_config.R", "Methods, typology")
put("threshPct", "concentration threshold for BL", 100 * cfg$thresh, "%",
    "51_config.R", "Methods, typology")
put("tukeyConf", "Tukey HSD confidence level", 100 * cfg$tukey_conf, "%",
    "51_config.R", "Fig. 4 caption")

# ---- transition matrix -------------------------------------------------------
transitions <- tp %>% count(class, classF, name = "n")
stayed <- sum(transitions$n[transitions$class == transitions$classF])
sayf("Class transitions: %d of %d dyads (%.0f%%) keep their class under the future scenario.\n",
     stayed, n, 100 * stayed / n)
put("nUnchanged", "dyads keeping their class in the future scenario", stayed,
    "dyads", "future classification", "Results, typology")
methods_numbers <- bind_rows(mn)

# ---- sample flow -------------------------------------------------------------
raw_n <- nrow(read.csv(cfg$dyad_file, stringsAsFactors = FALSE))
sample_flow <- tibble(
  step = c("dyads in FinalDiads.csv",
           "one row per (aquifer, country pair)",
           "classified (rule is exhaustive)",
           "entering Fig. 4B/C statistics (land use > 0)"),
  n = c(raw_n, n, sum(!is.na(tp$class)), sum(tp$land_total > 0, na.rm = TRUE)))
sample_flow$dropped <- c(NA_integer_, diff(sample_flow$n))
print(as.data.frame(sample_flow), row.names = FALSE)

# ---- robustness summary for the SI text -------------------------------------
rob_summary <- rob$partition %>%
  mutate(interior_share = round(100 * interior / n, 1),
         urban_share    = round(100 * urban / n, 1),
         bilateral_share = round(100 * bilateral / n, 1))

summ <- list(si_dyads = si_dyads, si2 = si2, methods_numbers = methods_numbers,
             transitions = transitions, sample_flow = sample_flow,
             rob_summary = rob_summary, bl_profile = rob$bl_profile)
saveRDS(summ, file.path(cfg$cache_dir, "summaries.rds"))
sayf("\n55_summaries.R done in %.2f min.\n",
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
