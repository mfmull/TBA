# =============================================================================
# 57_tables.R -- CANONICAL TABLES
#
#   Dataset_S4.csv               SI dyad-level dataset: classification under
#                                current and future demand, dyad
#                                characteristics, and the stage-by-stage build
#                                (was Data_S2.csv + the separate stage ledger)
#   table_main_typology.*        class counts, shares and transitions
#   table_si_group_stats.*       Fig. 4B/C means, SEs, CIs, Tukey letters
#   table_si_contingency.*       mosaic tables with their tests
#   table_si_robustness.*        class counts under every variant
#   table_si_sample_flow.*       row accounting
#   methods_numbers.*            every count quoted in the Methods
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!exists("CFG")) source(file.path(.root, "51_config.R"))
if (!exists("classify_dyads")) source(file.path(CFG$root, "52_core.R"))
if (is.null(.log$path)) log_open(CFG$log_file)
cfg <- CFG
t0 <- Sys.time()

say("\n================ TABLES =====================\n")
main <- readRDS(file.path(cfg$cache_dir, "main_objects.rds"))
rob  <- readRDS(file.path(cfg$cache_dir, "robustness_objects.rds"))
summ <- readRDS(file.path(cfg$cache_dir, "summaries.rds"))
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

# Minimal LaTeX writer: no dependency on a table package, and the .tex files
# stay diffable.
write_pair <- function(df, name, caption) {
  readr::write_csv(df, file.path(cfg$out_dir, paste0(name, ".csv")))
  fmt <- function(x) {
    if (is.numeric(x)) formatC(x, format = "g", digits = 4)
    else gsub("([&%_#])", "\\\\\\1", as.character(x))
  }
  body <- apply(as.data.frame(lapply(df, fmt)), 1,
                function(r) paste(r, collapse = " & "))
  tex <- c("\\begin{table}[htbp]", "\\centering", "\\small",
           sprintf("\\caption{%s}", caption),
           sprintf("\\begin{tabular}{%s}", paste(rep("l", ncol(df)),
                                                 collapse = "")),
           "\\hline",
           paste(paste(gsub("_", "\\\\_", names(df)), collapse = " & "), "\\\\"),
           "\\hline",
           paste(body, "\\\\"),
           "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(tex, file.path(cfg$out_dir, paste0(name, ".tex")))
  sayf("  wrote %s.csv / .tex (%d rows)\n", name, nrow(df))
}

# ---- SI Dataset S4: the consolidated dyad-level table ------------------------
readr::write_csv(summ$si_dyads, file.path(cfg$out_dir, cfg$si_dataset_name))
sayf("  wrote %s (%d dyads, %d columns)\n", cfg$si_dataset_name,
     nrow(summ$si_dyads), ncol(summ$si_dyads))

# ---- main table: counts, shares, transitions ---------------------------------
cur <- main$counts$current; fut <- main$counts$future; n <- main$counts$n_dyads
tab_main <- tibble(
  class = cfg$class_levels,
  n_current = as.integer(cur[cfg$class_levels]),
  share_current_pct = round(100 * as.integer(cur[cfg$class_levels]) / n, 1),
  n_future = as.integer(fut[cfg$class_levels]),
  share_future_pct = round(100 * as.integer(fut[cfg$class_levels]) / n, 1)) %>%
  mutate(change = n_future - n_current)
write_pair(tab_main, "table_main_typology",
           sprintf(paste("Dyad configuration under the preferred rule and",
                         "under the water-stressed cropland scenario",
                         "(%d dyads)."), n))

# ---- Fig. 4B/C statistics ----------------------------------------------------
gs <- main$tukey %>%
  mutate(metric = ifelse(Variable == "FGW", "national relevance",
                         "overdraft exposure")) %>%
  select(metric, class, n, mean, se, ci_low, ci_high, letters,
         any_of(c("boot_lo", "boot_hi")), kruskal_p)
write_pair(gs, "table_si_group_stats",
           sprintf(paste("Class means for the two Fig. 4 metrics. Letters are",
                         "Tukey HSD groupings at %.0f\\%%; classes sharing a",
                         "letter are not significantly different. Kruskal-Wallis",
                         "p is a distribution-free check on the same",
                         "comparison. Bootstrap bounds, where present, resample",
                         "countries rather than dyads."), 100 * cfg$tukey_conf))

# ---- contingency tables ------------------------------------------------------
ct <- bind_rows(lapply(names(main$mosaics), function(nm) {
  m <- main$mosaics[[nm]]
  as.data.frame(m$table, stringsAsFactors = FALSE) %>%
    setNames(c("class", "level", "n")) %>%
    mutate(table = m$title, .before = 1)
}))
write_pair(ct, "table_si_contingency",
           "Class by cooperation level, distance class and river-border status.")
ct_tests <- bind_rows(lapply(names(main$mosaics), function(nm)
  main$mosaics[[nm]]$test %>% mutate(table = main$mosaics[[nm]]$title,
                                     .before = 1)))
write_pair(ct_tests, "table_si_contingency_tests",
           paste("Tests of association for the tables above. Monte Carlo",
                 "p-values are used because several cells have expected",
                 "counts below five, where the asymptotic chi-square",
                 "approximation is unreliable."))

# ---- robustness --------------------------------------------------------------
rb <- bind_rows(lapply(names(rob$variants), function(nm) {
  r <- rob$variants[[nm]]
  tibble(variant = nm, si_figure = r$si, description = r$label,
         !!!setNames(as.list(r$counts$current), cfg$class_levels))
}))
rb <- bind_rows(
  tibble(variant = "preferred", si_figure = "Fig. 4",
         description = "Finalised specification (stages 1-3)",
         !!!setNames(as.list(as.integer(cfg$reference_class[cfg$class_levels])),
                     cfg$class_levels)),
  rb)
write_pair(rb, "table_si_robustness",
           paste("Dyad counts by configuration under each robustness variant.",
                 "Each variant changes exactly one setting of the preferred",
                 "rule."))
write_pair(summ$bl_profile, "table_si_bl_profile",
           paste("Overdraft exposure of the border-localised class under each",
                 "variant, with its Tukey grouping. BL remains the",
                 "highest-exposure class throughout."))

# ---- sample flow and methods numbers -----------------------------------------
write_pair(summ$sample_flow, "table_si_sample_flow",
           "Row accounting from the dyad table to the analysis sample.")
write_pair(summ$methods_numbers, "methods_numbers",
           "Every count quoted in the typology section, generated by the pipeline.")

sayf("\n57_tables.R done in %.2f min.\n",
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
