# =============================================================================
# 47_tables.R -- CANONICAL TABLES (CSV + TEX; output/)
#
#   table_main_results          best-balanced member per main specification:
#                               lmer estimate + model (Satterthwaite) SE/p AND
#                               weighted-OLS estimate + CR2 country-clustered
#                               SE/p (Table S5 replacement)
#   table_si_robustness         same layout, SI + river-exclusion blocks
#   table_si_ensemble           ensemble summaries (mean/median effect,
#                               5-95% band, share significant by sign)
#   table_si_balance            per-covariate SMD pre/post, best member
#
# All LaTeX is written directly (booktabs); no modelsummary dependency.
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!exists("CFG")) source(file.path(.root, "41_config.R"))
if (!exists("run_spec")) source(file.path(CFG$root, "42_core.R"))
if (is.null(.log$path)) log_open(CFG$log_file)
cfg <- CFG
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

main <- readRDS(file.path(cfg$cache_dir, "main_objects.rds"))
rob  <- readRDS(file.path(cfg$cache_dir, "robustness_objects.rds"))

fmt <- function(x, d = 4) ifelse(is.na(x), "--", formatC(x, digits = d, format = "fg"))
stars <- function(p) ifelse(is.na(p), "",
  ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", ""))))

results_row <- function(r) {
  b <- r$best
  tibble(
    specification = r$spec$title,
    label      = r$label,
    outcome    = r$spec$outcome,
    est_lmer   = b$treat_eff,  se_model = b$treat_se,  p_model = b$treat_p,
    est_ols    = b$treat_eff_ols, se_cr2 = b$treat_se_cr2, p_cr2 = b$treat_p_cr2,
    df_cr2     = b$cr2_df,
    int_lmer   = b$int_eff,   int_se = b$int_se,   int_p = b$int_p,
    n_treat    = b$n_treat,   n_control = b$n_control,
    n_country  = b$n_country,
    mean_smd   = b$mean_smd,  max_smd = b$max_smd, ess_ratio = b$ess_ratio,
    same_country_w = b$same_country_weighted,
    realization = b$idx, realization_legacy_rule = r$best_idx_legacy_rule)
}

tab_tex <- function(df, path, caption_note) {
  hdr <- c("\\begin{tabular}{lrrrrrrrr}", "\\toprule",
    paste0("Specification & \\multicolumn{2}{c}{Mixed model} ",
           "& \\multicolumn{2}{c}{Weighted OLS, CR2} & $N_t$ & $N_c$ ",
           "& mean $|$SMD$|$ & max $|$SMD$|$ \\\\"),
    "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
    " & Estimate & (SE) & Estimate & (SE) & & & & \\\\", "\\midrule")
  body <- sprintf("%s & %s%s & (%s) & %s%s & (%s) & %d & %d & %s & %s \\\\",
    gsub("[|]", "$\\\\mid$", df$specification),
    fmt(df$est_lmer), stars(df$p_model), fmt(df$se_model),
    fmt(df$est_ols),  stars(df$p_cr2),   fmt(df$se_cr2),
    df$n_treat, df$n_control, fmt(df$mean_smd, 2), fmt(df$max_smd, 2))
  ftr <- c("\\bottomrule", "\\end{tabular}",
    paste0("% ", caption_note))
  writeLines(c(hdr, body, ftr), path)
}

# ---- main results (Table S5 replacement) ------------------------------------
tm <- bind_rows(lapply(main, results_row))
write.csv(tm, file.path(cfg$out_dir, "table_main_results.csv"), row.names = FALSE)
tab_tex(tm, file.path(cfg$out_dir, "table_main_results.tex"),
  paste("Best-balanced ensemble member (lowest max |SMD|, PS distance row",
        "excluded). Mixed model: weighted lmer, country random intercept,",
        "Satterthwaite p. OLS/CR2: same weighted fixed-effect regression,",
        "CR2 SEs clustered on country. Stars: * p<0.1 ** p<0.05 *** p<0.01."))
say("  wrote table_main_results.csv/.tex\n")

# ---- SI + robustness --------------------------------------------------------
tr_ <- bind_rows(lapply(rob, results_row))
write.csv(tr_, file.path(cfg$out_dir, "table_si_robustness.csv"), row.names = FALSE)
tab_tex(tr_, file.path(cfg$out_dir, "table_si_robustness.tex"),
  "SI robustness blocks; same conventions as table_main_results.")
say("  wrote table_si_robustness.csv/.tex\n")

# ---- ensemble summaries -----------------------------------------------------
ens_row <- function(r) {
  cbind(tibble(specification = r$spec$title, label = r$label,
               outcome = r$spec$outcome, group = r$spec$group), r$counts)
}
te <- bind_rows(lapply(c(main, rob), ens_row))
write.csv(te, file.path(cfg$out_dir, "table_si_ensemble.csv"), row.names = FALSE)
tex <- c("\\begin{tabular}{lrrrrrr}", "\\toprule",
  paste0("Specification & Mean eff. & Median & q05 & q95 ",
         "& \\% sig.\\ $+$ & \\% sig.\\ $-$ \\\\"), "\\midrule",
  sprintf("%s & %s & %s & %s & %s & %d & %d \\\\",
          gsub("[|]", "$\\\\mid$", te$specification),
          fmt(te$mean_eff), fmt(te$med_eff), fmt(te$q05_eff), fmt(te$q95_eff),
          round(100 * te$share_pos_sig), round(100 * te$share_neg_sig)),
  "\\bottomrule", "\\end{tabular}",
  sprintf("%% Shares of the %d control realizations with p < %.1f (Satterthwaite).",
          cfg$n_realizations, cfg$sig_level))
writeLines(tex, file.path(cfg$out_dir, "table_si_ensemble.tex"))
say("  wrote table_si_ensemble.csv/.tex\n")

# ---- balance table ----------------------------------------------------------
tb <- bind_rows(lapply(c(main, rob), function(r) {
  bt <- r$best$balance
  tibble(specification = r$spec$title, label = r$label,
         covariate = bt$covariate, smd_pre = bt$smd_pre, smd_post = bt$smd_post,
         flag = ifelse(bt$covariate != "distance" &
                         abs(bt$smd_post) > cfg$balance_tol_substantive, "X", ""))
}))
write.csv(tb, file.path(cfg$out_dir, "table_si_balance.csv"), row.names = FALSE)
tex <- c("\\begin{tabular}{llrrl}", "\\toprule",
  "Specification & Covariate & SMD before & SMD after & Flag \\\\", "\\midrule",
  sprintf("%s & %s & %s & %s & %s \\\\",
          gsub("[|]", "$\\\\mid$", tb$specification),
          gsub("_", "\\\\_", tb$covariate),
          fmt(tb$smd_pre, 3), fmt(tb$smd_post, 3), tb$flag),
  "\\bottomrule", "\\end{tabular}",
  sprintf("%% Best-balanced member per specification; flag: |SMD| > %.2f.",
          cfg$balance_tol_substantive))
writeLines(tex, file.path(cfg$out_dir, "table_si_balance.tex"))
say("  wrote table_si_balance.csv/.tex\n")

# ---- SI summary statistics (from SI Dataset S3) -----------------------------
# Replaces legacy TableSummary.R, which (a) used the stale Data_S2 and
# (b) multiplied EVERY continuous variable by 100 -- correct for area
# fractions but not for CSI or road distance -- and (c) reported
# the LBRiv count from the inverted legacy RiverBorder flag.
#
# 2026-08-19 (FIX viii): the source file is now output/Dataset_S3.csv, which
# carries only the variables the analysis uses. Two rows are therefore gone
# from this table: area_km2 (the control-size caliper is applied upstream, in
# the control construction, and never enters a specification) and G_IrrigNeed
# (descriptive; no specification uses it, and it was an exact duplicate of the
# equally unused gini_crpGWSpct_gt0). Every remaining row is an outcome, a
# sample filter, or a matching/regression covariate.
s2 <- read.csv(file.path(cfg$out_dir, cfg$si_dataset_name))
sum_one <- function(v, scale = 1, count_zeros = TRUE) {
  x <- s2[[v]] * scale
  xp <- x[x > 0 & !is.na(x)]
  tibble(Variable = v,
         Zeros = if (count_zeros) sum(x == 0, na.rm = TRUE) else NA_integer_,
         Q25 = quantile(xp, .25), Median = median(xp), Q75 = quantile(xp, .75),
         NAs = if (count_zeros) NA_integer_ else sum(is.na(s2[[v]])))
}
gini_one <- function(v) {
  x <- s2[[v]]
  tibble(Variable = v, Zeros = NA_integer_,
         Q25 = quantile(x, .25, na.rm = TRUE), Median = median(x, na.rm = TRUE),
         Q75 = quantile(x, .75, na.rm = TRUE), NAs = sum(is.na(x)))
}
sumstats <- bind_rows(
  sum_one("Irrig",     100), sum_one("GWIrrig", 100),
  sum_one("IrrigNeed", 100), sum_one("Overdraft", 100),
  sum_one("OverIR6",   100),
  gini_one("G_Irrig"), gini_one("G_IrrigGW"), gini_one("gini_gwpct_gt0"),
  sum_one("CSI"), gini_one("G_CSI"),
  sum_one("RDS_mainDist"),
  tibble(Variable = "NearRiverBorder", Zeros = NA_integer_,
         Q25 = NA_real_, Median = NA_real_, Q75 = NA_real_,
         NAs = sum(s2$NearRiverBorder, na.rm = TRUE)))  # count TRUE, in NAs col slot
write.csv(sumstats, file.path(cfg$out_dir, "table_si_summary_stats.csv"),
          row.names = FALSE)
f2 <- function(x) ifelse(is.na(x), "",
  ifelse(abs(x) >= 1000, formatC(x, format = "f", digits = 0, big.mark = ","),
         formatC(x, format = "f", digits = 2)))
tex <- c("\\begin{tabular}{lccccc}", "\\toprule",
  "Variable & Zeros & Q25 & Median & Q75 & NAs/Count \\\\", "\\midrule",
  sprintf("%s & %s & %s & %s & %s & %s \\\\",
          gsub("_", "\\\\_", sumstats$Variable),
          ifelse(is.na(sumstats$Zeros), "", sumstats$Zeros),
          f2(sumstats$Q25), f2(sumstats$Median), f2(sumstats$Q75),
          ifelse(is.na(sumstats$NAs), "", sumstats$NAs)),
  "\\bottomrule", "\\end{tabular}",
  "% Irrig/GWIrrig/IrrigNeed/Overdraft/OverIR6 in percent of segment area",
  "% (x100), quartiles over positive values; CSI raw (0-10,000);",
  "% RDS_mainDist in km; NearRiverBorder row: count of TRUE (within 50 km",
  "% of a GSRB river border). From SI Dataset S3 (Dataset_S3.csv).")
writeLines(tex, file.path(cfg$out_dir, "table_si_summary_stats.tex"))
say("  wrote table_si_summary_stats.csv/.tex\n")

# =============================================================================
# MANUSCRIPT NUMBERS -- console report + trace (NOT a published output)
# =============================================================================
# Replaces fill_tex_numbers.R, which rewrote %%PLACEHOLDER%% tokens inside
# tex/_MS_PNAS_v2b.tex and tex/_SI_v2.tex in place. That script is retired for
# two reasons: the tex tree has moved to the repository root, and editing the
# manuscript from the pipeline makes the provenance of a number invisible in
# the diff. The two LaTeX tables it also regenerated are already produced
# above -- table_main_results is the Table S5 replacement, and
# table_si_summary_stats the summary-statistics body.
#
# What is left is the numbers quoted INLINE in the text. They are printed here
# for copying, and written to derived/audit/manuscript_numbers.csv as a trace.
# The trace is deliberately NOT in output/ and NOT in the canonical manifest:
# it is a record of what this run produced, not a publishable artefact.
say("\n---- numbers quoted in the manuscript text ----\n")
dir.create(cfg$audit_dir, showWarnings = FALSE, recursive = TRUE)

pct0 <- function(x) sprintf("%d", round(100 * x))
n1   <- function(x) formatC(x, digits = 1, format = "f")
n2   <- function(x) formatC(x, digits = 2, format = "f")
lx   <- lorenz_examples(cfg)

base_pct <- mean(main$IrIntens$summary_df$int_eff)
mn <- list()
add <- function(token, value, description, where) {
  mn[[length(mn) + 1]] <<- data.frame(
    token = token, value = as.character(value),
    description = description, appears_in = where,
    stringsAsFactors = FALSE)
}
add("EFF_PP",    n1(100 * main$IrIntens$counts$mean_eff),
    "ensemble mean transboundary difference in irrigated area (percentage points)",
    "Results, irrigation")
add("BASE_PCT",  n1(100 * base_pct),
    "matched-control baseline irrigated area (%)", "Results, irrigation")
add("EFF_RATIO", n1(1 + main$IrIntens$counts$mean_eff / base_pct),
    "transboundary irrigated area as a multiple of the control baseline",
    "Results, irrigation")
add("IRR_SHARE",       pct0(main$IrIntens$counts$share_sig),
    "share of realizations with a significant irrigation difference",
    "Results / Fig. 3")
add("OD_SHARE",        pct0(main$Overdraft$counts$share_sig),
    "share of realizations with a significant overdraft difference",
    "Results / Fig. 3")
# The Gini shares follow the text's meaning: significantly NEGATIVE, i.e.
# borderward, not merely significant in either direction.
add("IRGINI_SHARE",    pct0(main$IrGini$counts$share_neg_sig),
    "share of realizations with a significantly borderward irrigation Gini",
    "Results / Fig. 3")
add("IRRIVGINI_SHARE", pct0(main$IrRivsGini$counts$share_neg_sig),
    "same, conditional on river borders", "Results / Fig. 3")
add("GWGINI_SHARE",    pct0(main$GWGini$counts$share_neg_sig),
    "share of realizations with a significantly borderward GW-irrigation Gini",
    "Results / Fig. 3")
add("GWRIVGINI_SHARE", pct0(main$GWRivsGini$counts$share_neg_sig),
    "same, conditional on river borders", "Results / Fig. 3")
for (i in seq_along(lx))
  add(sprintf("G_EX%d", i), sprintf("%+.2f", lx[[i]]$gini),
      sprintf("Lorenz example: %s (%s)", lx[[i]]$name, lx[[i]]$country),
      "Fig. 3C")
add("SC_UNW", n2(100 * main$IrIntens$best$same_country_unweighted),
    "same-country share among matched controls, unweighted (%)",
    "Fig. 3 caption")
add("SC_W",   n2(100 * main$IrIntens$best$same_country_weighted),
    "same-country share among matched controls, weighted (%)",
    "Fig. 3 caption")
# Denominators, so a share in the text is never ambiguous about its base.
add("N_REALIZATIONS", cfg$n_realizations,
    "control realizations per specification", "Methods")
add("N_TREAT_IR", main$IrIntens$n_treat,
    "treated segments, Ir > 0 specifications", "Methods")
add("N_TREAT_GW", main$GWGini$n_treat,
    "treated segments, GW-Gini model sample", "Methods")

mnum <- do.call(rbind, mn)
print(as.data.frame(mnum[, c("token", "value", "description")]), row.names = FALSE,
      right = FALSE)
write.csv(mnum, file.path(cfg$audit_dir, "manuscript_numbers.csv"),
          row.names = FALSE)
sayf("\n  trace written: derived/audit/manuscript_numbers.csv (%d values, not published)\n",
     nrow(mnum))

say("47_tables done.\n")
