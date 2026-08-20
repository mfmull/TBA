# =============================================================================
# 35_summaries.R -- METHODS NUMBERS, SAMPLE FLOW, DESCRIPTIVES
#
# Builds every quantity quoted in the Methods and the descriptive statistics
# needed to write the paper, all from the fitted objects -- no number is typed
# by hand. Results are cached to derived/cache/summaries.rds; 36/37 render.
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!exists("CFG")) source(file.path(.root, "31_config.R"))
if (!exists("run_specification")) source(file.path(CFG$root, "32_core.R"))
if (is.null(.log$path)) log_open(CFG$log_file)
cfg <- CFG
t0 <- Sys.time()

main <- readRDS(file.path(cfg$cache_dir, "main_objects.rds"))
rob  <- readRDS(file.path(cfg$cache_dir, "robustness_objects.rds"))
if (isTRUE(cfg$cache_check)) {
  st_now <- cache_stamp(cfg)$data
  for (obj in list(main, rob))
    if (!identical(attr(obj, "stamp")$data, st_now))
      stop("Cached objects were computed on different input data; re-run 33/34.")
}
pref <- main$pref; d <- pref$data; h3 <- main$h3; ident <- main$ident
bsum <- main$balance_summary

# =============================================================================
# SAMPLE FLOW (wells)
# =============================================================================
raw <- read.csv(cfg$wells_file, stringsAsFactors = FALSE) %>% as_tibble()
n_raw    <- nrow(raw)
w1 <- raw %>% filter(is.finite(lon), is.finite(lat))
w2 <- w1 %>% filter(is.finite(GWSlp))
w3 <- w2 %>% filter(!is.na(Aquifer), nzchar(Aquifer), !is.na(CC),
                    !is.na(TB), is.finite(dist_LB_km))
w4 <- w3 %>% filter(!CC %in% cfg$drop_cc)
w4 <- w4 %>% mutate(unit_id = paste0(Aquifer, "_", CC))
seg_n <- w4 %>% count(unit_id)
w5 <- w4 %>% inner_join(seg_n %>% filter(n >= cfg$min_wells) %>% select(unit_id),
                        by = "unit_id")
w6 <- w5 %>% filter(unit_id %in% d$unit_id)   # matching + outcome eligibility
sample_flow <- tibble(
  step = c("raw wells in source file",
           "valid coordinates",
           "valid depletion outcome (GWSlp)",
           "valid aquifer-country assignment and border distance",
           "country with a land border (islands dropped)",
           sprintf("segment has >= %d wells", cfg$min_wells),
           "segment passes outcome and matching-covariate eligibility [final H1 = H2 sample]"),
  wells = c(n_raw, nrow(w1), nrow(w2), nrow(w3), nrow(w4), nrow(w5), nrow(w6))) %>%
  mutate(dropped = c(NA, -diff(wells)))
# H1 and H2 use the same final well sample by construction (assertion):
stopifnot(sum(d$n) == nrow(w6))

# =============================================================================
# DESCRIPTIVES
# =============================================================================
wf <- w6 %>% left_join(d %>% select(unit_id, TBn), by = "unit_id")
miss_tab <- raw %>%
  summarise(across(c(GWSlp, dist_LB_km, lon, lat, urbkHaKm2, CS_max,
                     LB_river, prec_mm),
                   ~ sum(!is.finite(.x)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(share_missing = n_missing / n_raw)

dist_desc <- wf %>%
  mutate(group = if_else(TBn == 1, "transboundary", "domestic")) %>%
  group_by(group) %>%
  summarise(n_wells = n(),
            median_dist_km = median(dist_LB_km),
            q25_dist_km = quantile(dist_LB_km, .25),
            q75_dist_km = quantile(dist_LB_km, .75),
            share_within_10km = mean(dist_LB_km <= 10),
            share_within_50km = mean(dist_LB_km <= 50),
            share_within_100km = mean(dist_LB_km <= 100),
            .groups = "drop")

seg_desc <- d %>%
  mutate(group = if_else(TBn == 1, "transboundary", "domestic")) %>%
  group_by(group) %>%
  summarise(n_segments = n(),
            n_countries = n_distinct(CC),
            n_aquifers = n_distinct(aq_id),
            wells_total = sum(n),
            wells_median = median(n),
            wells_q25 = quantile(n, .25), wells_q75 = quantile(n, .75),
            .groups = "drop")

cc_conc <- d %>% group_by(CC) %>%
  summarise(segments = n(), treated = sum(TBn == 1), wells = sum(n),
            .groups = "drop") %>%
  arrange(desc(wells)) %>%
  mutate(share_wells = wells / sum(wells))

cov_desc <- d %>%
  mutate(group = if_else(TBn == 1, "transboundary", "domestic")) %>%
  group_by(group) %>%
  summarise(across(all_of(attr(d, "ps_cov")),
                   list(mean = ~ mean(.x), sd = ~ sd(.x)),
                   .names = "{.col}__{.fn}"),
            .groups = "drop") %>%
  pivot_longer(-group, names_to = c("covariate", "stat"), names_sep = "__") %>%
  pivot_wider(names_from = c(group, stat), values_from = value)

tba_cov <- d %>% filter(TBn == 1) %>%
  summarise(across(starts_with("share_within_"),
                   list(seg_any = ~ mean(.x > 0), mean_share = ~ mean(.x))))

overlap_cc <- d %>% group_by(CC) %>%
  summarise(k = n_distinct(TBn), wells = sum(n), .groups = "drop")
n_overlap_cc <- sum(overlap_cc$k > 1)
wells_overlap <- sum(overlap_cc$wells[overlap_cc$k > 1])

# =============================================================================
# METHODS NUMBERS (canonical machine-readable file)
# =============================================================================
NUM <- list()
put <- function(key, label, value, unit, source_object, manuscript_location,
                pr = "preferred") {
  NUM[[length(NUM) + 1]] <<- tibble(
    key = key, label = label, value = value, unit = unit,
    source_object = source_object, manuscript_location = manuscript_location,
    preferred_or_robustness = pr)
  invisible(value)
}
tb_aq <- unique(d$aq_id[d$TBn == 1])
usa_wells <- cc_conc$wells[cc_conc$CC == "USA"]

put("nWellsRaw", "raw wells in source file", n_raw, "wells", "sample_flow", "Methods, sample")
put("nWellsCoord", "wells with valid coordinates", nrow(w1), "wells", "sample_flow", "Methods, sample")
put("nWellsOutcome", "wells with valid depletion outcome", nrow(w2), "wells", "sample_flow", "Methods, sample")
put("nWellsAssigned", "wells with valid aquifer-country assignment", nrow(w3), "wells", "sample_flow", "Methods, sample")
put("nWellsFinal", "final analysis wells (H1 and H2)", nrow(w6), "wells", "sample_flow", "Methods, sample")
put("nWellsTB", "transboundary wells", sum(d$n[d$TBn == 1]), "wells", "segment table", "Methods, sample")
put("nWellsCtrl", "domestic-control wells", sum(d$n[d$TBn == 0]), "wells", "segment table", "Methods, sample")
put("nSegments", "aquifer-country segments", nrow(d), "segments", "segment table", "Methods, sample")
put("nSegmentsTB", "transboundary segments", sum(d$TBn == 1), "segments", "segment table", "Methods, sample")
put("nSegmentsCtrl", "domestic segments", sum(d$TBn == 0), "segments", "segment table", "Methods, sample")
put("nAquifers", "physical aquifers", n_distinct(d$aq_id), "aquifers", "segment table", "Methods, sample")
put("nAquifersTB", "transboundary physical aquifers", length(tb_aq), "aquifers", "segment table", "Methods, sample")
put("nCountries", "countries", n_distinct(d$CC), "countries", "segment table", "Methods, sample")
put("nCountriesOverlap", "countries with both treated and control segments",
    n_overlap_cc, "countries", "overlap_cc", "Methods/SI, identification")
put("nWellsOverlapCC", "wells in overlap countries", wells_overlap, "wells", "overlap_cc", "SI, identification")
put("shareWellsOverlapCC", "share of wells in overlap countries",
    round(wells_overlap / sum(d$n), 3), "share", "overlap_cc", "SI, identification")
put("nWellsUSA", "wells in the United States", usa_wells %||% NA_real_, "wells",
    "cc_conc", "SI, concentration")
put("shareWellsUSA", "share of wells in the United States",
    round((usa_wells %||% NA_real_) / sum(d$n), 3), "share", "cc_conc", "SI, concentration")
for (km in cfg$coverage_km) {
  v <- wf %>% filter(TBn == 1)
  put(sprintf("nTBwells%dkm", km), sprintf("TBA wells within %d km of border", km),
      sum(v$dist_LB_km <= km), "wells", "well table", "Results, localization")
  put(sprintf("shareTBwells%dkm", km),
      sprintf("share of TBA wells within %d km", km),
      round(mean(v$dist_LB_km <= km), 3), "share", "well table", "Results, localization")
  s <- d %>% filter(TBn == 1)
  put(sprintf("nTBseg%dkm", km), sprintf("TBA segments with wells within %d km", km),
      sum(s[[paste0("share_within_", km, "km")]] > 0), "segments",
      "segment table", "Results, localization")
  put(sprintf("shareTBseg%dkm", km),
      sprintf("share of TBA segments with wells within %d km", km),
      round(mean(s[[paste0("share_within_", km, "km")]] > 0), 3), "share",
      "segment table", "Results, localization")
}
put("medianWellsPerSeg", "median wells per segment", median(d$n), "wells", "segment table", "Methods, sample")
put("iqrWellsPerSegLo", "wells per segment, 25th pct", quantile(d$n, .25), "wells", "segment table", "Methods, sample")
put("iqrWellsPerSegHi", "wells per segment, 75th pct", quantile(d$n, .75), "wells", "segment table", "Methods, sample")
put("minWellsThreshold", "minimum wells per segment", cfg$min_wells, "wells", "CFG", "Methods, sample")
put("essControlATO", "effective control sample (ATO)", round(bsum$ess_control_ato, 1),
    "segments", "balance_summary", "Methods, design")
put("essTreatedATO", "effective treated sample (ATO)", round(bsum$ess_treated_ato, 1),
    "segments", "balance_summary", "Methods, design")
put("matchingEstimand", "matching estimand", NA_real_, cfg$estimand, "CFG", "Methods, design")
put("maxSMDATO", "max |SMD| under overlap weights", signif(bsum$max_abs_smd_ato, 2),
    "SMD", "balance_summary", "Methods, design")
put("maxSMDunweighted", "max |SMD| unweighted", signif(bsum$max_abs_smd_unweighted, 2),
    "SMD", "balance_summary", "Methods, design")
put("nCountryClusters", "country clusters for inference", n_distinct(d$CC),
    "clusters", "segment table", "Methods, inference")
put("nWCB", "wild cluster bootstrap repetitions", cfg$n_wcb, "reps", "CFG", "Methods, inference")
put("nProfileBoot", "country bootstrap repetitions (profiles)", cfg$profile_boot,
    "reps", "CFG", "Methods, inference")
put("nH3Boot", "H3 rank-stability bootstrap repetitions", cfg$h3_boot_reps,
    "reps", "CFG", "SI, H3")
put("conleyCutoff", "Conley cutoff", cfg$conley_km, "km", "CFG", "Methods, first stage")
put("firstStageEstimator", "first-stage estimator", NA_real_,
    "segment mean depletion (HC3) and Fisher z of depletion-distance correlation",
    "CFG", "Methods, first stage")
put("secondStageEstimator", "second-stage estimator", NA_real_,
    "weighted mixed-effects meta-regression (rma.mv, REML), random = country + segment",
    "CFG", "Methods, second stage")
put("weightDefinition", "second-stage weight", NA_real_,
    "w_ATO / (sampling variance + tau2), winsorised", "CFG", "Methods, second stage")
put("tauH1", "between-segment SD, H1 (level)", round(sqrt(pref$h1$tau2), 1),
    "mm/yr", "pref", "Results, H3")
put("tauH2", "between-segment SD, H2 (gradient)", round(sqrt(pref$h2$tau2), 3),
    "Fisher z", "pref", "Results, H3")
put("H1estimate", "H1 preferred estimate", round(pref$h1$estimate, 1), "mm/yr",
    "pref", "Results, H1")
put("H1ciLow", "H1 90% CI lower", round(pref$h1$ci90_low, 1), "mm/yr", "pref", "Results, H1")
put("H1ciHigh", "H1 90% CI upper", round(pref$h1$ci90_high, 1), "mm/yr", "pref", "Results, H1")
put("H1pOneSided", "H1 one-sided p (predicted direction)", round(pref$h1$p_one, 3),
    "p", "pref", "Results, H1")
put("H1mde80", "H1 80%-power MDE (one-sided alpha 0.05)", round(main$mde80, 0),
    "mm/yr", "pref", "Results, H1 detection limit")
put("H1mde90", "H1 90%-power MDE (one-sided alpha 0.05)", round(main$mde90, 0),
    "mm/yr", "pref", "SI, H1 detection limit")
put("H1mdeSource", "MDE uncertainty source", NA_real_, pref$h1$se_source,
    "pref", "Methods, inference")
put("H2estimate", "H2 preferred estimate", round(pref$h2$estimate, 4), "Fisher z",
    "pref", "Results, H2")
put("H2ciLow", "H2 90% CI lower", round(pref$h2$ci90_low, 4), "Fisher z", "pref", "Results, H2")
put("H2ciHigh", "H2 90% CI upper", round(pref$h2$ci90_high, 4), "Fisher z", "pref", "Results, H2")
put("H2pOneSided", "H2 one-sided p (predicted direction)", round(pref$h2$p_one, 3),
    "p", "pref", "Results, H2")
put("H2placebo", "H2 intercept (placebo-type diagnostic)",
    round(pref$h2$baseline, 4), "Fisher z", "pref", "Results, H2")
put("H2placeboP", "H2 placebo two-sided p", round(pref$h2$baseline_p, 3), "p",
    "pref", "Results, H2")
put("nH3Segments", "H3 transboundary segments", nrow(h3$scatter), "segments",
    "h3", "Results, H3")
put("nSelectedAquifers", "aquifers selected by the coded rule",
    nrow(ident$selected), "segments", "ident", "Results, H3")
methods_numbers <- bind_rows(NUM)

# =============================================================================
# BLOCK-B TABLE: consolidated H1 robustness (one row per family member)
# =============================================================================
tab <- rob$tab
h1_rob <- tab %>%
  filter(ok, !(h1_same %in% TRUE)) %>%
  transmute(family, specification, estimand,
            h1_estimate,
            h1_ci90_low, h1_ci90_high,
            h1_p_one, h1_p_one_source,
            # Bootstrap interval and the one-sided bound that inverts the same
            # test as h1_p_one. The CR2 interval above is kept for transparency,
            # but it is not the one to print beside h1_p_one.
            h1_wcb_lo, h1_wcb_hi, h1_wcb_bound_one, h1_wcb_bound_one_alt,
            uncertainty_source = h1_se_source,
            h1_mde80, h1_mde90,
            n_treated, n_control, n_countries, n_wells,
            status = if_else(converged, "converged", "check convergence"))

# =============================================================================
# BLOCK-A TABLE: H2 robustness
# =============================================================================
h2_rob <- tab %>%
  filter(ok) %>%
  transmute(family, specification, estimand,
            h2_unit, h2_scale, table_only,
            h2_estimate, h2_ci90_low, h2_ci90_high,
            h2_p_one, h2_p_cr2, h2_p_wild,
            h2_wcb_lo, h2_wcb_hi, h2_wcb_bound_one, h2_wcb_bound_one_alt,
            h2_baseline, h2_baseline_p,
            n_treated, n_control, n_countries,
            uncertainty_estimator_matched, ratio_source,
            status = if_else(converged, "converged", "check convergence"))

# family-level summaries
fam_sum <- function(t, est) {
  t %>% group_by(family) %>%
    summarise(n_specs = n(),
              n_negative = sum(.data[[est]] < 0),
              n_positive = sum(.data[[est]] > 0),
              min = min(.data[[est]]), max = max(.data[[est]]),
              median = median(.data[[est]]), .groups = "drop")
}
h2_fam <- fam_sum(h2_rob %>% filter(h2_scale == "z"), "h2_estimate")
h1_fam <- fam_sum(h1_rob, "h1_estimate")

# =============================================================================
# persist
# =============================================================================
summaries <- list(
  sample_flow = sample_flow, missing = miss_tab, dist_desc = dist_desc,
  seg_desc = seg_desc, cc_conc = cc_conc, cov_desc = cov_desc,
  tba_cov = tba_cov, methods_numbers = methods_numbers,
  h1_rob = h1_rob, h2_rob = h2_rob, h1_fam = h1_fam, h2_fam = h2_fam,
  n_overlap_cc = n_overlap_cc, wells_overlap = wells_overlap)
cache_write(summaries, file.path(cfg$cache_dir, "summaries.rds"),
            cache_stamp(cfg, list(object = "summaries")))
sayf("35_summaries.R done in %.1f min.\n",
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
