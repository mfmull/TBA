# =============================================================================
# 45_summaries.R -- METHODS NUMBERS, SAMPLE FLOW, SI DATASET S3
#
#   output/methods_numbers.csv / .tex    every count quoted in the Methods
#   output/table_si_sample_flow.csv/.tex sample construction, treated+controls
#   output/Dataset_S3.csv                SI supplementary dataset, segment
#                                        level (FIX vii + viii; was Data_S2)
#
# Writes derived/cache/summaries.rds for 46-48.
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!exists("CFG")) source(file.path(.root, "41_config.R"))
if (!exists("run_spec")) source(file.path(CFG$root, "42_core.R"))
if (is.null(.log$path)) log_open(CFG$log_file)
cfg <- CFG
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

d    <- load_main_data(cfg)
sets <- load_control_sets(cfg)
tr   <- d %>% filter(type == "treat")
cn   <- d %>% filter(type != "treat")

# ---- SI Dataset S3 (regenerated from the current build; FIX vii + viii) -----
s2 <- build_si_irrigation_dataset(cfg)
write.csv(s2, file.path(cfg$out_dir, cfg$si_dataset_name), row.names = FALSE)
say("  wrote ", cfg$si_dataset_name, " (", nrow(s2), " rows, ",
    ncol(s2), " columns)\n")

# ---- methods numbers --------------------------------------------------------
mn <- tribble(
  ~quantity, ~value, ~where_quoted,
  "TBA-country segments (treated)",            nrow(tr),                       "Methods: Transboundary aquifers",
  "Control basins in pool",                    nrow(cn),                       "Methods: controls",
  "Treated segments with irrigated cropland",  sum(tr$Ir > 0, na.rm = TRUE),   "Methods: 'restricted to segments overlapping irrigated cropland'",
  "Treated segments with nonzero GW irrigation", sum(tr$GW > 0, na.rm = TRUE), "Methods (UPDATED: was 207 from the pre-rebuild Data_S2)",
  "Control basins with irrigated cropland",    sum(cn$Ir > 0, na.rm = TRUE),   "Methods: controls",
  "Control basins with nonzero GW irrigation", sum(cn$GW > 0, na.rm = TRUE),   "Methods: controls",
  "Non-overlapping control realizations",      length(sets),                   "Methods: ensemble",
  "Median control-set size",                   median(lengths(sets)),          "Methods: ensemble",
  "Treated segments within 50 km of a river border", sum(tr$near_river_border, na.rm = TRUE), "Methods: '| Rivers' adjustment",
  "Control basins within 50 km of a river border",   sum(cn$near_river_border, na.rm = TRUE), "Methods: '| Rivers' adjustment",
  "Countries represented (treated segments)",  n_distinct(tr$CntrName),        "Methods: random intercepts")
write.csv(mn, file.path(cfg$out_dir, "methods_numbers.csv"), row.names = FALSE)

tex <- c("\\begin{tabular}{lrl}", "\\toprule",
         "Quantity & Value & Quoted in \\\\", "\\midrule",
         sprintf("%s & %s & %s \\\\",
                 gsub("&", "\\\\&", mn$quantity),
                 format(mn$value, big.mark = ","),
                 gsub("[_&]", " ", mn$where_quoted)),
         "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(cfg$out_dir, "methods_numbers.tex"))
say("  wrote methods_numbers.csv/.tex\n")

# ---- sample flow ------------------------------------------------------------
flow <- tribble(
  ~step, ~treated, ~controls,
  "Master table (_dataMain.csv)",            nrow(tr), nrow(cn),
  "Irrigated (Ir > 0): intensity + Ir Gini", sum(tr$Ir > 0, na.rm = TRUE), sum(cn$Ir > 0, na.rm = TRUE),
  "GW-irrigated (GW > 0): GW Gini",          sum(tr$GW > 0, na.rm = TRUE), sum(cn$GW > 0, na.rm = TRUE),
  "Ir > 0, not near river border (excl. robustness)",
    sum(tr$Ir > 0 & !tr$near_river_border, na.rm = TRUE),
    sum(cn$Ir > 0 & !cn$near_river_border, na.rm = TRUE),
  "GW > 0, not near river border (excl. robustness)",
    sum(tr$GW > 0 & !tr$near_river_border, na.rm = TRUE),
    sum(cn$GW > 0 & !cn$near_river_border, na.rm = TRUE))
write.csv(flow, file.path(cfg$out_dir, "table_si_sample_flow.csv"), row.names = FALSE)
tex <- c("\\begin{tabular}{lrr}", "\\toprule",
         "Sample restriction & Treated segments & Control basins \\\\", "\\midrule",
         sprintf("%s & %d & %d \\\\", gsub("[_&>]", " ", flow$step),
                 flow$treated, flow$controls),
         "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(cfg$out_dir, "table_si_sample_flow.tex"))
say("  wrote table_si_sample_flow.csv/.tex\n")

summ <- list(methods_numbers = mn, sample_flow = flow,
             n_gw_treated = sum(tr$GW > 0, na.rm = TRUE),
             si_dataset_rows = nrow(s2),
             si_dataset_cols = names(s2),
             data_s2_rows = nrow(s2))   # legacy name, kept for old drivers
saveRDS(summ, file.path(cfg$cache_dir, "summaries.rds"))
