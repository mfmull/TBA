# =============================================================================
# 34_run_robustness.R -- ROBUSTNESS, ORGANISED AROUND THE NARRATIVE
#
# Four inferential blocks:
#   A  robustness of H2 (the main positive result)
#   B  robustness of H1 and its detection limit (consolidated by family)
#   C  robustness of the descriptive localization pattern (no bin tests)
#   D  robustness of H3 and the selected aquifers
#
# Every check is a call to run_specification() with overrides; nothing
# re-implements fitting code. Each specification is cached to
# derived/cache/specs/<id>.rds with a stamp covering the configuration, the
# override list, the code and the input data; stale caches are refused.
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!exists("CFG")) source(file.path(.root, "31_config.R"))
if (!exists("run_specification")) source(file.path(CFG$root, "32_core.R"))
if (is.null(.log$path)) log_open(CFG$log_file)
# FORCE=TRUE rebuilds every stamped cache below. 38_run_all.R resolves it from
# the environment or the calling session and defines it before sourcing this
# file, so only supply a default when running this script on its own.
if (!exists("FORCE")) FORCE <- nzchar(Sys.getenv("FORCE"))
if (isTRUE(FORCE)) say("*** FORCE: rebuilding cached specifications. ***\n")
# SMOKE mode (SMOKE=1 environment variable): truncates the registry and the
# leave-one-out loops to a few entries so the whole pipeline can be exercised
# end-to-end in minutes. Smoke outputs are NOT quotable; the canonical run has
# SMOKE unset. (Successor of the former 25_smoke_test.R.)
SMOKE <- nzchar(Sys.getenv("SMOKE"))
if (SMOKE) say("*** SMOKE MODE: truncated registry and loops; not quotable. ***\n")
t0 <- Sys.time()

cfg <- CFG
spec_dir <- file.path(cfg$cache_dir, "specs")
dir.create(spec_dir, showWarnings = FALSE, recursive = TRUE)

ref_path <- file.path(cfg$cache_dir, "main_reference.rds")
if (!file.exists(ref_path))
  stop("Run 33_run_main.R first: no frozen reference to compare against.")
ref <- readRDS(ref_path)
ref_key <- .md5_of(list(h1 = ref$h1, h2 = ref$h2, n = ref$n_segments))

say("\n================ ROBUSTNESS (blocks A-D) =====================\n")
sayf("Frozen preferred: H1 %+.4f mm/yr, H2 %+.6f Fisher z\n", ref$h1, ref$h2)

# ---- registry ---------------------------------------------------------------
# `family` groups conceptual specifications; parameter variants (Conley
# cutoffs, subsample seeds, winsorisation, minimum-well cutoffs) live INSIDE a
# family and are never counted as independent confirmations.
# `table_only`: H2 outcome not on the Fisher-z scale (excluded from H2 forest).
# `h1_same`: override touches only the H2 outcome; H1 identical to preferred.
REG <- list()
add <- function(id, family, label, ..., table_only = FALSE, h1_same = FALSE) {
  REG[[length(REG) + 1]] <<- list(spec_id = id, family = family,
                                  specification = label, overrides = list(...),
                                  table_only = table_only, h1_same = h1_same)
}

# Empirical common distance support (0.01/0.99 quantile overlap of treated and
# control well-to-border distance distributions).
common_support_window <- function(cfg, probs = c(0.01, 0.99)) {
  d <- build_segments(cfg)
  wells <- fs_wells(cfg) %>% inner_join(d %>% select(unit_id, TBn), by = "unit_id")
  qs <- wells %>% group_by(TBn) %>%
    summarise(lo = quantile(dist_LB_km, probs[1]),
              hi = quantile(dist_LB_km, probs[2]), .groups = "drop")
  c(max(qs$lo), min(qs$hi))
}
sw <- tryCatch(common_support_window(cfg), error = function(e) NULL)
if (!is.null(sw)) sayf("Empirical common distance support: %.1f-%.1f km\n", sw[1], sw[2])

## Design and target population
add("pref",   "Design", "Preferred specification (ATO)")
add("cfe",    "Design", "Country fixed effects (within-country identification)",
    country_effect = "fixed")
add("att",    "Design", "ATT estimand, nearest-neighbour matching (different estimand)",
    estimand = "ATT")
add("noriver","Design", "Exclude river-border segments",
    exclude_river_borders = TRUE)
## Border-crossing rivers (Marc; R5.4). Distinct from `noriver`: that row
## removes segments whose border IS a river, this one removes aquifers whose
## two national parts are connected BY a river. Classification comes from
## classify_river_crossings.R (GloRiC v1.0); see SI Section S3. The dropped-id
## vector is sorted (deterministic ov_md5); if the classification file is
## absent the row is skipped so the pipeline still runs without the river data.
.rc_file <- file.path(cfg$root, "1_data", "river_crossings_by_unit.csv")
if (file.exists(.rc_file)) {
  .rc  <- read.csv(.rc_file, stringsAsFactors = FALSE)
  .q   <- unique(.rc$q_min_cms)
  stopifnot(length(.q) == 1L)
  # aquifer-level promotion: a crossing river connects both national parts by
  # construction, so any flagged segment flags its whole aquifer.
  .aq  <- unique(.rc$Aquifer[.rc$has_crossing_river])
  .ids <- sort(unique(.rc$unit_id[.rc$Aquifer %in% .aq]))
  add("nocrossriver", "Design",
      sprintf("Exclude aquifers with a border-crossing river (Q >= %g m3/s)", .q),
      drop_treated_segments = .ids)
} else {
  say("  [skipped] nocrossriver: 1_data/river_crossings_by_unit.csv not found\n")
}
## Support diagnostics -- NOT causal confirmations of a near-border H1.
## Reported with their control counts; if near-border control support is thin,
## these rows are read as support/sensitivity diagnostics only (block C notes).
add("supp_fix", "Support diagnostics", "Wells within 0-100 km only (support diagnostic)",
    distance_window_km = c(0, 100))
if (!is.null(sw))
  add("supp_emp", "Support diagnostics",
      sprintf("Empirical common distance support (%.0f-%.0f km)", sw[1], sw[2]),
      distance_window_km = sw)
## First stage
add("spear",  "First stage", "Spearman rank correlation for H2",
    first_stage_h2 = "spearman", h1_same = TRUE)
add("robfs",  "First stage", "Robust (Huber M) first-stage level (approx. spatial SE)",
    first_stage_h1 = "robust_mean")
add("trimfs", "First stage", "Trimmed-mean first-stage level (approx. spatial SE)",
    first_stage_h1 = "trimmed_mean")
add("slope",  "First stage", "H2 as physical slope, mm/yr/km (table only)",
    first_stage_h2 = "physical_slope", table_only = TRUE, h1_same = TRUE)
add("bins",   "First stage", "H2 as near (0-10 km) vs interior (10-100 km) contrast (table only)",
    first_stage_h2 = "distance_bins", table_only = TRUE, h1_same = TRUE)
add("sdscreen", "First stage", "Stricter distance-variation screen (SD >= 10 km)",
    min_sd_dist_km = 10)
## Spatial uncertainty (one family)
for (k in c(5, 25, 50))
  add(sprintf("conley%d", k), "Spatial uncertainty",
      sprintf("Conley cutoff %d km", k), conley_km = k)
for (p in c(0, 1))
  add(sprintf("neff%d", p), "Spatial uncertainty",
      sprintf("Effective sample size n/rho^%d", p), n_eff_power = p,
      h1_same = TRUE)
for (s in 2:6)
  add(sprintf("seed%02d", s), "Spatial uncertainty",
      sprintf("Dense-segment Conley subsample, seed %d", s), conley_seed = s)
## Second stage and parameters (one family)
add("nowins",  "Second stage and parameters", "No variance winsorisation",
    variance_winsor = NULL)
add("wins595", "Second stage and parameters", "Variance winsorisation at 5th/95th",
    variance_winsor = c(0.05, 0.95))
add("nosegre", "Second stage and parameters", "No segment random effect",
    include_segment_re = FALSE)
for (m in c(50, 100))
  add(sprintf("minw%d", m), "Second stage and parameters",
      sprintf("Minimum %d wells per segment", m), min_wells = m)
add("unadj",   "Second stage and parameters", "Unadjusted outcome model (~ 1 + TBn)",
    outcome_adjustment = FALSE)
add("augcov",  "Second stage and parameters",
    "Augmented covariate set (+ precipitation)",
    ps_covariates = c(cfg$cov_match, cfg$cov_extra))
add("equalinfo", "Second stage and parameters",
    "ATO-weighted, equal segment information",
    precision_weight = FALSE)

if (SMOKE) REG <- REG[vapply(REG, function(e)
  e$spec_id %in% c("pref", "att", "supp_fix", "spear", "slope", "nosegre"),
  logical(1))]
sayf("%d registry specifications, %d families.\n", length(REG),
     n_distinct(vapply(REG, function(e) e$family, "")))

# ---- pre-build every Conley combination the registry needs ------------------
# (sequential, so the parallel spec loop below never writes the ratio cache)
prebuild <- unique(c(
  lapply(c(5, 10, 25, 50), function(k) list(conley_km = k, conley_seed = 1)),
  lapply(2:6, function(s) list(conley_km = 10, conley_seed = s))))
if (SMOKE) prebuild <- prebuild[2]
for (pb in prebuild) {
  cfgx <- cfg_with(cfg, conley_km = pb$conley_km, conley_seed = pb$conley_seed)
  invisible(get_ratios(cfgx$conley_km, cfgx, wells = fs_wells(cfgx)))
}
for (win in Filter(Negate(is.null), list(c(0, 100), sw))) {
  cfgx <- cfg_with(cfg, distance_window_km = win)
  invisible(get_ratios(cfgx$conley_km, cfgx, wells = fs_wells(cfgx)))
}
# Warm the first-stage memoisation before forking.
invisible(first_stage_from_wells(cfg))

# ---- run the registry -------------------------------------------------------
# ALL cache I/O happens in this (parent) process; the forked workers compute
# and return, and touch no files. This matters: the previous version stamped,
# read and wrote inside each worker, which meant eight forked children were
# concurrently md5-summing wellsData.csv (37 MB) and writing .rds files on a
# cloud-synced volume. That is the one par_map call site in the pipeline whose
# workers did file I/O, and the only one whose workers died en masse.
spec_file <- function(e) file.path(spec_dir, paste0(e$spec_id, ".rds"))
# cache_stamp() hashes the configuration, the source files and wellsData.csv.
# Those three are IDENTICAL for all 30 specifications -- only the override and
# the extra vary -- so hash the files once and vary the two cheap fields. The
# resulting stamps are byte-for-byte what the previous code produced, so
# caches written by earlier runs still validate.
.base_stamp <- cache_stamp(cfg, list(), extra = NULL)

# ---- test-inversion bootstrap intervals for the registry --------------------
# The registry inherits cfg$wcb_ci = FALSE -- 33_run_main.R switches the
# bootstrap interval on for the preferred specification only -- which is why
# every registry row cached before 2026-08-13 carries NA in wcb_lo, wcb_hi and
# wcb_bound_one. Table S2 needs them for every row: with 19 country clusters the
# reported p_1 comes from the bootstrap, and pairing it with a CR2 interval lets
# the two contradict each other. They do in the preferred row, where H2 has
# p_1 = 0.045 beside a CR2 90% interval of [-0.160, 0.002] that contains zero;
# the bootstrap interval for that row is [-0.124, -0.002].
#
# This is applied as a per-entry OVERRIDE rather than by flipping CFG$wcb_ci,
# and the distinction is the point. An override reaches the cache key through
# ov_md5 below, so exactly the 30 registry specifications refit. Changing CFG
# would move cfg_md5, which .base_stamp carries and every other stamp in this
# file inherits, dragging the 103 influence refits, locC, the H3 bootstrap and
# the selection-stability loops along with it -- none of which appear in
# Table S2, and none of which would gain anything.
REG_CI_OV <- list(wcb_ci = TRUE)

spec_stamp <- function(e) {
  st <- .base_stamp
  st$ov_md5    <- .md5_of(.canon_list(modifyList(e$overrides, REG_CI_OV)))
  st$extra_md5 <- .md5_of(list(ref = ref_key, spec = e$specification))
  st
}
spec_stamps <- lapply(REG, spec_stamp)
spec_files  <- vapply(REG, spec_file, character(1))

# Compute only: safe to fork.
run_spec <- function(e) {
  r <- do.call(run_specification,
               c(list(spec_id = e$spec_id, specification = e$specification,
                      family = e$family), modifyList(e$overrides, REG_CI_OV)))
  r$table_only <- isTRUE(e$table_only)
  r$data <- NULL; r$s1 <- NULL; r$s2 <- NULL; r$h3 <- NULL
  r
}
results <- lapply(seq_along(REG), function(i)
  cache_read(spec_files[i], spec_stamps[[i]], cfg, force = FORCE))
todo <- which(vapply(results, is.null, logical(1)))
sayf("  registry: %d of %d cached, %d to compute.\n",
     length(REG) - length(todo), length(REG), length(todo))
if (length(todo)) {
  fresh <- par_map(REG[todo], run_spec, label = "registry", cfg = cfg)
  for (k in seq_along(todo)) {
    i <- todo[k]
    results[[i]] <- fresh[[k]]
    # A job that still failed after the sequential retry must not be cached as
    # if it were a result; it will surface as a failed row below instead.
    if (!is.null(fresh[[k]]))
      cache_write(fresh[[k]], spec_files[i], spec_stamps[[i]])
  }
}
# The mutate() below matches results to REG by position, so alignment is an
# assumption worth asserting rather than discovering three tables later.
stopifnot(length(results) == length(REG))

# spec_row() is defined in 32_core.R and is deliberately NOT extended to carry
# the bootstrap bounds. That file is hashed into cache_stamp()$code, so editing
# it would invalidate every cache in the pipeline rather than just the registry.
# The bounds are read straight off the packed h1/h2 lists, which run_spec()
# keeps (it drops only data, s1, s2 and h3).
wcb_row <- function(r) {
  if (!is.list(r)) r <- list()
  g  <- function(x, f) { v <- x[[f]]
    if (is.null(v) || !length(v)) NA_real_ else as.numeric(v[1]) }
  gc <- function(x, f) { v <- x[[f]]
    if (is.null(v) || !length(v)) NA_character_ else as.character(v[1]) }
  tibble(
    h1_wcb_lo = g(r$h1, "wcb_lo"), h1_wcb_hi = g(r$h1, "wcb_hi"),
    h1_wcb_bound_one     = g(r$h1, "wcb_bound_one"),
    h1_wcb_bound_one_alt = gc(r$h1, "wcb_bound_one_alt"),
    h2_wcb_lo = g(r$h2, "wcb_lo"), h2_wcb_hi = g(r$h2, "wcb_hi"),
    h2_wcb_bound_one     = g(r$h2, "wcb_bound_one"),
    h2_wcb_bound_one_alt = gc(r$h2, "wcb_bound_one_alt"))
}
tab <- bind_cols(bind_rows(lapply(results, spec_row)),
                 bind_rows(lapply(results, wcb_row))) %>%
  mutate(table_only = vapply(REG, function(e) isTRUE(e$table_only), logical(1)),
         h1_same    = vapply(REG, function(e) isTRUE(e$h1_same),    logical(1)),
         spec_id    = coalesce(spec_id,
                               vapply(REG, function(e) e$spec_id, "")))
stopifnot(nrow(tab) == length(REG),
          identical(tab$spec_id, vapply(REG, function(e) e$spec_id, "")))

# A stale cache written before REG_CI_OV existed would return rows whose
# bootstrap columns are entirely NA, and an all-NA column reaches Table S2 as a
# column of blanks rather than as an error. Fail loudly instead: the preferred
# row must carry a finite bootstrap interval for both hypotheses.
.pref_wcb <- tab %>% filter(spec_id == "pref")
if (!nrow(.pref_wcb) || !all(is.finite(c(.pref_wcb$h1_wcb_lo, .pref_wcb$h1_wcb_hi,
                                         .pref_wcb$h2_wcb_lo, .pref_wcb$h2_wcb_hi))))
  stop("Registry ran without test-inversion bootstrap intervals. ",
       "Check REG_CI_OV, and delete derived/cache/specs/ if the cache predates it.")
sayf("  bootstrap intervals present for %d of %d registry rows (H1), %d (H2).\n",
     sum(is.finite(tab$h1_wcb_lo)), nrow(tab), sum(is.finite(tab$h2_wcb_lo)))

failed <- tab %>% filter(!ok)
if (nrow(failed)) {
  say("\n*** FAILED SPECIFICATIONS (reported, not dropped) ***\n")
  print(as.data.frame(failed %>% select(spec_id, specification, warnings)),
        row.names = FALSE)
}

# preferred row must reproduce the frozen reference
pref_row <- tab %>% filter(spec_id == "pref")
stopifnot(nrow(pref_row) == 1)
if (abs(pref_row$h1_estimate - ref$h1) > cfg$reproduce_tol * max(1, abs(ref$h1)) ||
    abs(pref_row$h2_estimate - ref$h2) > cfg$reproduce_tol * max(1, abs(ref$h2)))
  stop("Registry 'pref' does not reproduce the frozen main result.")
say("\nPreferred specification reproduces the frozen main result.\n")

# ---- influence: leave-one-out loops (parallel; ratios cached) ---------------
base_d <- add_design_weights(build_segments(cfg), cfg)
influence_loop <- function(name, ids, override_name, extra_of) {
  f  <- file.path(cfg$cache_dir, paste0("influence_", name, ".rds"))
  st <- cache_stamp(cfg, list(override = override_name),
                    extra = list(ref = ref_key, ids = sort(as.character(ids))))
  hit <- cache_read(f, st, cfg, force = FORCE)
  if (!is.null(hit)) { say("  [cached] influence: ", name, "\n"); return(hit) }
  one <- function(i) {
    id <- ids[i]
    ov <- list(id); names(ov) <- override_name
    r <- do.call(run_specification,
                 c(list(spec_id = paste0(name, "_", i),
                        specification = paste("drop", id),
                        family = "Influence"), ov))
    bind_cols(tibble(dropped = id), extra_of(id), spec_row(r) %>%
                select(-spec_id, -family, -specification))
  }
  out <- par_map_dfr(seq_along(ids), one,
                     label = paste0("leave-one-", name, "-out"), cfg = cfg) %>%
    mutate(d_h1 = h1_estimate - ref$h1,
           d_h2 = h2_estimate - ref$h2,
           sign_reversal_h2 = is.finite(h2_estimate) &
             sign(h2_estimate) != sign(ref$h2)) %>%
    arrange(desc(abs(d_h2)))
  cache_write(out, f, st)
  out
}
.trim <- function(ids) if (SMOKE) head(ids, 3) else ids
inf_cc <- influence_loop("country", .trim(sort(unique(base_d$CC))), "drop_countries",
  function(id) tibble(
    wells_removed = sum(base_d$n[base_d$CC == id], na.rm = TRUE),
    treated_removed = sum(base_d$CC == id & base_d$TBn == 1)))
inf_seg <- influence_loop("treated_segment",
  .trim(sort(base_d$unit_id[base_d$TBn == 1])), "drop_treated_segments",
  function(id) tibble(
    wells_removed = base_d$n[match(id, base_d$unit_id)],
    treated_removed = 1L))
inf_aq <- influence_loop("physical_aquifer",
  .trim(sort(unique(base_d$aq_id[base_d$TBn == 1]))), "drop_aquifers",
  function(id) tibble(
    segments_removed = sum(base_d$aq_id == id),
    wells_removed = sum(base_d$n[base_d$aq_id == id], na.rm = TRUE),
    treated_removed = sum(base_d$aq_id == id & base_d$TBn == 1)))

for (x in list(list(inf_cc, "country"), list(inf_seg, "treated segment"),
               list(inf_aq, "physical aquifer"))) {
  nrev <- sum(x[[1]]$sign_reversal_h2, na.rm = TRUE)
  sayf("Leave-one-%s-out: %d H2 sign reversal(s) of %d.\n",
       x[[2]], nrev, nrow(x[[1]]))
}

# ---- block C: robustness of the descriptive localization pattern ------------
say("\n---- Block C: localization robustness (descriptive) ----\n")
locC_path <- file.path(cfg$cache_dir, "localization_robustness.rds")
locC_st <- cache_stamp(cfg, list(analysis = "locC"), extra = list(ref = ref_key))
locC <- cache_read(locC_path, locC_st, cfg, force = FORCE)
if (is.null(locC)) {
  d <- base_d
  # influential countries/aquifers for exclusion checks: largest treated well counts
  big_cc <- base_d %>% filter(TBn == 1) %>% count(CC, wt = n, sort = TRUE) %>%
    slice_head(n = 2) %>% pull(CC)
  big_aq <- base_d %>% filter(TBn == 1) %>% arrange(desc(n)) %>%
    slice_head(n = 2) %>% pull(aq_id)
  specs <- list(
    list(lab = "preferred (weighted, centred trend)", args = list()),
    list(lab = "unweighted", args = list(weighted = FALSE)),
    list(lab = "uncentred trend", args = list(centred_trend = FALSE)),
    list(lab = "min 50 wells per segment",
         args = list(), cfgov = list(min_wells = 50)),
    list(lab = sprintf("drop %s", big_cc[1]),
         args = list(drop_cc_extra = big_cc[1])),
    list(lab = sprintf("drop %s", big_cc[2]),
         args = list(drop_cc_extra = big_cc[2])),
    list(lab = sprintf("drop largest aquifer (%s)", big_aq[1]),
         args = list(drop_aq_extra = big_aq[1])),
    list(lab = sprintf("drop 2nd-largest aquifer (%s)", big_aq[2]),
         args = list(drop_aq_extra = big_aq[2])))
  loc_tab <- map_dfr(specs, function(s) {
    cfgx <- if (!is.null(s$cfgov)) do.call(cfg_with, c(list(cfg), s$cfgov)) else cfg
    dx <- if (!is.null(s$cfgov)) add_design_weights(build_segments(cfgx), cfgx) else d
    do.call(localization_summary,
            c(list(d = dx, cfg = cfgx, label = s$lab), s$args))
  })
  # leave-one-country-out profile contrasts (compact: contrast only)
  ccs_tb <- sort(unique(base_d$CC[base_d$TBn == 1]))
  loo_prof <- map_dfr(ccs_tb, function(cc)
    localization_summary(d, cfg, drop_cc_extra = cc, boot = 199,
                         label = paste("drop", cc)))
  # alternative bins for the Panel B profile
  alt_bins <- list(
    default = cfg$profile_bins,
    coarse  = c(0, 25, 100, 250, Inf),
    fine    = c(0, 5, 10, 25, 50, 75, 100, 150, 200, 300, Inf))
  profiles <- imap(alt_bins, function(b, nm)
    profile_by_distance(d, cfg, bins = b) %>% mutate(bins = nm))
  # unweighted and centred variants of the default profile
  prof_unw <- profile_by_distance(d, cfg, weighted = FALSE) %>%
    mutate(bins = "default (unweighted)")
  prof_ctr <- profile_by_distance(d, cfg, centred = TRUE) %>%
    mutate(bins = "default (within-segment centred)")
  locC <- list(table = loc_tab, loo = loo_prof,
               profiles = bind_rows(profiles), prof_unw = prof_unw,
               prof_ctr = prof_ctr, windows = window_summary(d, cfg))
  cache_write(locC, locC_path, locC_st)
}

# ---- block D: robustness of H3 and the selected aquifers --------------------
say("\n---- Block D: H3 stability ----\n")
boot_path <- file.path(cfg$cache_dir, "h3_bootstrap.rds")
boot_st <- cache_stamp(cfg, list(analysis = "h3boot"), extra = list(ref = ref_key))
h3b <- cache_read(boot_path, boot_st, cfg, force = FORCE)
if (is.null(h3b)) {
  h3b <- h3_bootstrap(cfg)
  cache_write(h3b, boot_path, boot_st)
}
sayf("H3 bootstrap: %d of %d replicates valid.\n", h3b$n_valid, h3b$n_requested)

# Selection stability across robustness families: does the coded selection
# reproduce under each defensible full-specification variant?
selstab_path <- file.path(cfg$cache_dir, "h3_selection_stability.rds")
selstab_st <- cache_stamp(cfg, list(analysis = "selstab"),
                          extra = list(ref = ref_key))
selstab <- cache_read(selstab_path, selstab_st, cfg, force = FORCE)
if (is.null(selstab)) {
  fam_specs <- list(
    list(id = "pref",    lab = "preferred",              ov = list()),
    list(id = "unadj",   lab = "unadjusted outcome",     ov = list(outcome_adjustment = FALSE)),
    list(id = "nowins",  lab = "no winsorisation",       ov = list(variance_winsor = NULL)),
    list(id = "minw50",  lab = "min 50 wells",           ov = list(min_wells = 50)),
    list(id = "conley25",lab = "Conley 25 km",           ov = list(conley_km = 25)),
    list(id = "spear",   lab = "Spearman H2",            ov = list(first_stage_h2 = "spearman")),
    list(id = "noriver", lab = "no river borders",       ov = list(exclude_river_borders = TRUE)))
  if (SMOKE) fam_specs <- fam_specs[c(1, 2)]
  one_sel <- function(s) {
    cfgx <- do.call(cfg_with, c(list(cfg), s$ov))
    dx <- add_design_weights(build_segments(cfgx), cfgx)
    if (is.null(dx)) return(NULL)
    id <- tryCatch(identify_segments(dx, cfgx), error = function(e) NULL)
    if (is.null(id)) return(NULL)
    sc <- tryCatch(h3_analysis(dx, cfgx)$scatter, error = function(e) NULL)
    out <- id$table %>%
      transmute(family_spec = s$lab, .gid, label_txt, selected,
                rank_raw, consistent)
    if (!is.null(sc))
      out <- out %>% left_join(sc %>% select(.gid, pattern), by = ".gid")
    out
  }
  selstab <- par_map_dfr(fam_specs, one_sel, label = "H3 selection stability",
                         cfg = cfg)
  cache_write(selstab, selstab_path, selstab_st)
}



# =============================================================================
# AUGMENTED SPECIFICATIONS AND HYDRAULIC-CONNECTIVITY SUBSETS
#
# Folded in from the reviewer-response scripts of the revision round, which
# no longer exist as separate files. Everything below extends the registry
# built above and is written into the same `rob` object, so 36 and 37 read it
# exactly as they read any other robustness result.
#
#   * IGRAC hydraulic properties and the Fig. 2C crosswalk
#   * connectivity subsets            -> SI Fig. S2
#   * well-level regressions          -> SI Table S4 (two rows)
#   * governance-proxy augmentation   -> SI Table S4 (one row) and Table S5
#   * natural water-table depth       -> SI Table S4 (one row)
#   * crossing-river complement       -> SI Table S4, note c
#
# Each block self-skips when its optional input is absent, so the canonical
# paper outputs never depend on a file a user has not downloaded.
# =============================================================================


d0    <- d                                   # 39b called the preferred frame d0
wells <- fs_wells(cfg)
have_lme4      <- requireNamespace("lme4",        quietly = TRUE)
have_club      <- requireNamespace("clubSandwich", quietly = TRUE)
have_sandwich  <- requireNamespace("sandwich",    quietly = TRUE)
have_terra     <- requireNamespace("terra",       quietly = TRUE)
augmented_rows <- list()                     # collected -> rob$augmented


# ---- input resolver (1_data/ then map/) --------------------------------------

find_input <- function(name, dirs = c("1_data", "map"), cfg = CFG) {
  cand <- file.path(cfg$root, dirs, name)
  hit  <- cand[file.exists(cand)]
  if (!length(hit))
    stop("Required input '", name, "' not found in any of: ",
         paste(file.path(dirs), collapse = ", "),
         " (searched under ", cfg$root, ").")
  if (length(hit) > 1)
    notes_add(sprintf("input '%s' found in %d locations; using %s",
                      name, length(hit), dirname(hit[1])))
  sayf("  input %-24s <- %s/\n", name, basename(dirname(hit[1])))
  hit[1]
}

# ---- IGRAC hydraulic properties ---------------------------------------------

aqf <- read.csv(find_input("aqfData.csv"),
                stringsAsFactors = FALSE) %>% as_tibble()
aqf_need <- c("Aquifer", "CC", "area_km2", "GWIrKHaKm2", "crpkHaKm2",
              "urbkHaKm2", "prec_mm", "pop")
if (length(setdiff(aqf_need, names(aqf))))
  stop("aqfData.csv is missing required column(s): ",
       paste(setdiff(aqf_need, names(aqf)), collapse = ", "))
# ROI (100-yr radius of influence) is COMPUTED here from Tr and S, never read
# from the source file. The distributed IGRAC table carried its own ROI column
# that ran backwards (falling as Tr rises, e.g. Mimbres Tr = 0.0022 m^2/s ->
# "ROI" = 6310 m, vs. Tijuana-San Diego Tr = 0.0186 m^2/s -> "ROI" = 631 m),
# which is not physically sensible for a radius of influence and indicates the
# source column is not Eqn S22. That column has since been deleted from
# 1_data/IGRAC_Properties.csv; the select(-any_of(...)) below keeps the
# deletion enforced, so a re-downloaded copy cannot silently reintroduce it.
# Eqn S22 (SI): ROI = 1.5 * sqrt(Tr * t / S), with Tr in m^2/s, t = 100 yr in
# seconds, S dimensionless storativity.
# No other script in the 3x series reads this file or references ROI, Tr, S,
# Lithology or Confinement (verified); the recomputation is local to 39.
igrac_roi_years <- 100
ig_raw <- read.csv(find_input("IGRAC_Properties.csv"),
                   stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM") %>%
  as_tibble()
# The first column is the IGRAC aquifer code. Older copies of this file left it
# unnamed (or BOM-mangled), hence the positional rename -- but only rename when
# it is needed, so that a file which already names it aq_code, or which leads
# with some other column entirely, fails loudly instead of creating duplicates.
if (!identical(names(ig_raw)[1], "aq_code")) {
  if ("aq_code" %in% names(ig_raw))
    stop("IGRAC_Properties.csv: column 1 is '", names(ig_raw)[1],
         "' but an 'aq_code' column already exists at position ",
         match("aq_code", names(ig_raw)),
         ". This is not the raw IGRAC property table -- check the file.")
  ig_raw <- ig_raw %>% rename_with(~ "aq_code", 1)
}
ig_need <- c("aq_code", "aquifer", "Lithology", "Confinement.description.",
             "Tr", "S")
if (length(setdiff(ig_need, names(ig_raw))))
  stop("IGRAC_Properties.csv is missing required column(s): ",
       paste(setdiff(ig_need, names(ig_raw)), collapse = ", "))
ig <- ig_raw %>%
  select(-any_of(c("ROI", "diffusivity"))) %>%   # never inherit a tabulated ROI
  mutate(Tr = suppressWarnings(as.numeric(Tr)),
         S = suppressWarnings(as.numeric(S)),
         # Tr == 0 encodes "no data" in this file.
         Tr = if_else(Tr <= 0, NA_real_, Tr),
         S  = if_else(S <= 0, NA_real_, S),
         ROI = 1.5 * sqrt(Tr * (igrac_roi_years * 365.25 * 24 * 3600) / S),
         diffusivity = Tr / S)
sayf("  IGRAC properties: %d records, %d with usable Tr and S.\n",
     nrow(ig), sum(is.finite(ig$Tr) & is.finite(ig$S)))


# =============================================================================
# 2. IGRAC CROSSWALK (curated, deterministic)
# =============================================================================
# Jasechko aquifer units do not systematically coincide with IGRAC polygons
# (Methods), so the match below is by hydrogeologic correspondence and is
# graded: "direct" (same named system), "approximate" (nearest enclosing or
# adjacent IGRAC system), "none" (no plausible IGRAC counterpart in the
# 202-record property table). Only matched segments enter the connectivity
# bins; coverage is reported.
XW <- tribble(
  ~unit_id,                                                ~aq_code, ~quality,
  "Tijuana-San Diego Basin_USA",                           "8N",  "direct",
  "Mexicali Valley_USA",                                   "9N",  "direct",
  "Yuma Basin_USA",                                        "9N",  "approximate",
  "Eastern Imperial and Amos and Oligby Valleys_USA",      "9N",  "approximate",
  "Avra Valley_USA",                                       "12N", "approximate",
  "San Pedro Basin_USA",                                   "13N", "direct",
  "Mimbres Basin_USA",                                     "18N", "direct",
  "Valle de Juarez and Hueco Bolson_USA",                  "15N", "direct",
  "Rio Grande Delta_USA",                                  "17N", "direct",
  "Rio Grande Delta_MEX",                                  "17N", "direct",
  "Edwards Plateau_USA",                                   "16N", "direct",
  "Stockton Plateau_USA",                                  "16N", "approximate",
  "Western Carrizo-Wilcox_USA",                            NA,    "none",
  "Williston Basin_USA",                                   "6N",  "approximate",
  "Devonian and Overlying Unconsolidated Aquifers_CAN",    "6N",  "approximate",
  "Sandlands Aquifer_CAN",                                 "6N",  "approximate",
  "Southern Odanah Shale_CAN",                             "6N",  "approximate",
  "Transboundary Odanah Shale_CAN",                        "6N",  "approximate",
  "Southern St. Lawrence Lowlands_CAN",                    "7N",  "approximate",
  "Bari Doab_PAK",                                         "AS78", "direct",
  "Chaj Doab_PAK",                                         "AS78", "direct",
  "Rechna Doab_PAK",                                       "AS78", "direct",
  "Eastside Middle Indus Plain_PAK",                       "AS78", "direct",
  "Eastside Middle Indus Plain_IND",                       "AS78", "direct",
  "Trans Gangetic Plain_IND",                              "AS78", "approximate",
  "Upper Gangetic Plain_IND",                              "AS80", "approximate",
  "Middle Gangetic Plain_IND",                             "AS80", "approximate",
  "Western Bengal Basin_IND",                              "AS80", "approximate",
  "Barind Tract and Central Floodplains_IND",              "AS80", "approximate",
  "Barind Tract and Central Floodplains_BGD",              "AS80", "approximate",
  "Northern Piedmont and Tista Fan_IND",                   "AS80", "approximate",
  "Northern Piedmont and Tista Fan_BGD",                   "AS80", "approximate",
  "Madhupur Tract and Eastern Floodplains_BGD",            "AS80", "approximate",
  "Southcentral Floodplains_BGD",                          "AS80", "approximate",
  "Southeastern Bengal Basin_BGD",                         "AS80", "approximate",
  "Sylhet Basin_BGD",                                      "AS80", "approximate",
  "Eastern Mekong Delta_VNM",                              "AS89", "approximate",
  "Mioplioquatenario de Elvas-Campo Maior and Vegas Baja Aquifers_ESP",
                                                           "EU10", "approximate",
  "Pomfret-Vergelegen Dolomitic Aquifer_ZAF",              "AF6",  "approximate",
  "Mons Sedimentary Basin_BEL",                            NA,    "none",
  "Tournaisis Region Carbonate Aquifers_BEL",              NA,    "none",
  "Alsatian Aquifer and Upper Rhine Valley_FRA",           NA,    "none",
  "West Jutland Quaternary and Tertiary Sand Aquifer System_DNK", NA, "none",
  "Lom- Pleven Depression_BGR",                            NA,    "none")
stopifnot(all(d$unit_id[d$TBn == 1] %in% XW$unit_id))
xw <- XW %>%
  left_join(ig %>% select(aq_code, igrac_name = aquifer, Lithology,
                          Confinement = Confinement.description., Tr, S, ROI,
                          diffusivity),
            by = "aq_code") %>%
  mutate(lith_class = case_when(
    grepl("^Sediment -", Lithology) ~ "Unconsolidated sediment",
    grepl("^Sedimentary rocks", Lithology) ~ "Sedimentary rock",
    grepl("^Crystalline|^Metamorphic", Lithology) ~ "Crystalline / metamorphic",
    TRUE ~ NA_character_))
readr::write_csv(xw, file.path(cfg$out_dir, "table_si_igrac_crosswalk.csv"))
sayf("Crosswalk: %d direct, %d approximate, %d unmatched; %d with Tr/S data.\n",
     sum(xw$quality == "direct"), sum(xw$quality == "approximate"),
     sum(xw$quality == "none"), sum(is.finite(xw$ROI)))


# =============================================================================
# 3. H1/H2 BY HYDRAULIC CONNECTIVITY AND LITHOLOGY (R5.4, R3.3, R5.3, R1.2)
# =============================================================================
say("\n---- connectivity subsets ----\n")
tb_x <- d %>% filter(TBn == 1) %>% select(unit_id, n) %>%
  left_join(xw, by = "unit_id")
roi_ok <- tb_x %>% filter(is.finite(ROI))
roi_cut <- quantile(roi_ok$ROI, c(1/3, 2/3))
tb_x <- tb_x %>%
  mutate(roi_bin = case_when(
    !is.finite(ROI) ~ NA_character_,
    ROI <= roi_cut[1] ~ "short ROI (highest connectivity attenuation)",
    ROI <= roi_cut[2] ~ "intermediate ROI",
    TRUE ~ "long ROI"))
# NOTE ON DIRECTION. ROI (Eqn S22) increases with transmissivity and
# decreases with storativity, as physically expected for a drawdown radius:
# short ROI means the 100-yr cone of influence stays close to the well
# (highest attenuation, weakest far-field connectivity), long ROI means it
# reaches furthest. Bins are labelled by ROI itself, with Tr, S and
# diffusivity carried alongside for interpretation, rather than asserting a
# single "connectivity" ordering from Tr or diffusivity alone.
subsets <- c(
  split(tb_x$unit_id, tb_x$roi_bin),
  split(tb_x$unit_id, tb_x$lith_class))
subsets <- subsets[lengths(subsets) >= 3]  # need enough treated segments

# cache_stamp() md5-sums wellsData.csv (37 MB) on every call, and the file
# hashes are identical across subsets -- only the override varies. Hash once.
.base_stamp <- cache_stamp(cfg, list(), extra = NULL)
run_subset <- function(keep_ids, lab, sid) {
  drop_ids <- setdiff(d$unit_id[d$TBn == 1], keep_ids)
  f <- file.path(file.path(cfg$cache_dir, "specs_augmented"), paste0(sid, ".rds"))
  st <- .base_stamp
  st$ov_md5 <- .md5_of(.canon_list(list(subset = sort(keep_ids),
                                        wcb_ci = TRUE)))
  st$extra_md5 <- .md5_of(list(what = "connectivity_wcb"))
  hit <- cache_read(f, st, cfg, force = FORCE)
  if (!is.null(hit)) return(hit)
  # wcb_ci = TRUE: the wild-cluster-bootstrap CI (not just its p-value) is
  # needed here, because the figure plots it -- see note by conn_tab below.
  r <- run_specification(sid, lab, "Connectivity subsets",
                         drop_treated_segments = drop_ids, wcb_ci = TRUE)
  r$keep_ids <- keep_ids
  # keep r$data for the profile, drop the heavy fit objects
  r$s1 <- NULL; r$s2 <- NULL
  cache_write(r, f, st)
  r
}
subruns <- imap(subsets, function(ids, lab)
  run_subset(ids, lab, paste0("sub_", substr(gsub("[^a-z]", "", tolower(lab)), 1, 14))))

conn_tab <- imap_dfr(subruns, function(r, lab) {
  props <- tb_x %>% filter(unit_id %in% r$keep_ids)
  tibble(subset = lab,
         n_treated = r$sample$n_treated,
         n_wells_treated = sum(props$n),
         Tr_median = median(props$Tr, na.rm = TRUE),
         S_median = median(props$S, na.rm = TRUE),
         ROI_median = median(props$ROI, na.rm = TRUE),
         # 90% two-sided WILD CLUSTER BOOTSTRAP interval (test inversion) --
         # the SAME inference route as h*_p_one below, not the CR2/
         # Satterthwaite interval. Plotting a CR2 interval next to a WCB
         # p-value let the bar and p1 visibly disagree about significance
         # (e.g. a bar not crossing zero next to a non-significant p1, and
         # vice versa); this keeps bar and p1 consistent by construction.
         h1_estimate = r$h1$estimate, h1_lo = r$h1$wcb_lo,
         h1_hi = r$h1$wcb_hi, h1_p_one = r$h1$p_one,
         h2_estimate = r$h2$estimate, h2_lo = r$h2$wcb_lo,
         h2_hi = r$h2$wcb_hi, h2_p_one = r$h2$p_one,
         h2_placebo = r$h2$baseline, h2_placebo_p = r$h2$baseline_p,
         segments = paste(sub("_[A-Z]+$", "", r$keep_ids), collapse = "; "))
})
readr::write_csv(conn_tab, file.path(cfg$out_dir, "table_si_wells_connectivity.csv"))


# =============================================================================
# 5. WELL-LEVEL UNPOOLED REGRESSION (R1.4)
# =============================================================================
# GWSlp_i = a + g TBn + covariates + (1|CC) + (1|unit_id), each well carrying
# its segment's overlap weight divided by the segment's well count (so a
# segment's total influence equals its design weight). This is a qualitative
# sensitivity, not an ATO estimator: lme4 treats these as precision weights.
# Clustered (CR2) SEs are unavailable for weighted lmer fits (clubSandwich
# declines prior weights), so the wild cluster bootstrap on country is the
# clustered inference; the segment-clustered objection is addressed by the
# segment random intercept plus the country bootstrap.
say("\n---- well-level regression (R1.4) ----\n")
wl_tab <- NULL
if (have_lme4) {
  ps_cov <- attr(d, "ps_cov")
  keep <- d %>% select(unit_id, CC, TBn, w_ato, n_seg = n,
                       all_of(paste0(ps_cov, "_ctr")))
  wv <- wells %>% select(unit_id, GWSlp, dist_LB_km) %>%
    inner_join(keep, by = "unit_id") %>%
    group_by(unit_id) %>%
    mutate(dist_z = if (sd(dist_LB_km) > 0)
             (dist_LB_km - mean(dist_LB_km)) / sd(dist_LB_km) else NA_real_,
           gw_z = if (sd(GWSlp) > 0)
             (GWSlp - mean(GWSlp)) / sd(GWSlp) else NA_real_) %>%
    ungroup() %>%
    mutate(w_well = w_ato / n_seg)
  # The adjustment covariates span degrees of latitude (tens), CS_max
  # (thousands) and urban density (thousandths), which is what lme4's
  # "predictor variables are on very different scales" warning is about: it is
  # a conditioning complaint, not a validity one. Scaling each to unit SD fixes
  # the conditioning and leaves the reported terms untouched -- rescaling a
  # covariate rescales only its own coefficient, and TBn, dist_z and their
  # interaction are deliberately left on their original scales.
  cov_ctr <- paste0(ps_cov, "_ctr")
  wv <- wv %>% mutate(across(all_of(cov_ctr), function(v) {
    s <- stats::sd(v, na.rm = TRUE)
    if (is.finite(s) && s > 0) v / s else v
  }))
  fit_wl <- function(form, dat, term, lab, unit) {
    dat$.w <- dat$w_well
    environment(form) <- environment()
    f <- tryCatch(suppressMessages(lme4::lmer(
      form, data = dat, weights = .w,
      control = lme4::lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))),
      error = function(e) NULL)
    if (is.null(f)) return(NULL)
    cf <- lme4::fixef(f)
    if (!term %in% names(cf)) return(NULL)
    est <- unname(cf[term]); se <- sqrt(diag(as.matrix(vcov(f))))[term]
    X <- model.matrix(lme4::nobars(form), data = dat)
    wb <- wcb_p(dat[[all.vars(form)[1]]], X, dat$w_well, dat$CC,
                jname = term, cfg = cfg)
    tibble(outcome = lab, term = term, unit = unit, estimate = est,
           se_model = unname(se),
           ci90_low = est - 1.645 * se, ci90_high = est + 1.645 * se,
           p_model_two_sided = 2 * pnorm(abs(est / se), lower.tail = FALSE),
           p_wild_country_two_sided = if (is.null(wb)) NA_real_ else wb$p,
           n_wells = nrow(dat), n_segments = n_distinct(dat$unit_id),
           n_countries = n_distinct(dat$CC),
           note = paste("Precision-weighted mixed model; segment and country",
                        "random intercepts; CR2 unavailable for weighted lmer",
                        "fits, wild cluster bootstrap on country reported."))
  }
  cc_terms <- paste0(ps_cov, "_ctr")
  f1 <- reformulate(c("TBn", cc_terms, "(1 | CC)", "(1 | unit_id)"),
                    response = "GWSlp")
  f2 <- reformulate(c("TBn * dist_z", cc_terms, "(1 | CC)", "(1 | unit_id)"),
                    response = "gw_z")
  wl_tab <- bind_rows(
    fit_wl(f1, wv %>% filter(is.finite(GWSlp), w_well > 0), "TBn",
           "H1 analogue: well-level depletion", "mm/yr"),
    fit_wl(f2, wv %>% filter(is.finite(gw_z), is.finite(dist_z), w_well > 0),
           "TBn:dist_z",
           "H2 analogue: within-segment standardised gradient",
           "SD depletion per SD distance"))
  readr::write_csv(wl_tab, file.path(cfg$out_dir, "table_si_well_level.csv"))
  print(as.data.frame(wl_tab[, 1:8]), digits = 3, row.names = FALSE)
}


# =============================================================================
# SHARED: refit the two-stage design on an augmented segment frame
# =============================================================================
# Reproduces the run_specification() inner flow (recentre covariates at the
# transboundary mean of the estimation sample, refit propensity + overlap
# weights, second stage for H1 and H2) but takes the segment frame and the
# propensity covariate set as arguments, so a country- or segment-level
# variable that does not live in wellsData.csv can enter the design. The
# preferred covariate set is ALWAYS refit on the same complete-case subsample
# so that augmented-vs-preferred differences are never composition effects.

.aug_row <- function(s, hyp, unit, spec_label, dd, cfg) {
  if (is.null(s)) return(NULL)
  po <- one_sided_p(s, hyp, cfg)
  sp <- se_of(s)
  tibble::tibble(
    specification = spec_label, hypothesis = hyp, unit = unit,
    estimate  = s$est,
    ci90_low  = if (is.finite(s$cc_lo %||% NA_real_)) s$cc_lo else s$ci_lo,
    ci90_high = if (is.finite(s$cc_hi %||% NA_real_)) s$cc_hi else s$ci_hi,
    p_one_sided = po$p, p_one_source = po$source,
    p_two_cr2 = s$p_cc %||% NA_real_, p_two_wild = s$p_wcb %||% NA_real_,
    mde80 = if (hyp == "H1") mde_of(sp$se, 0.80, cfg) else NA_real_,
    baseline = s$b0 %||% NA_real_, tau2 = s$tau2 %||% NA_real_,
    n_segments = nrow(dd), n_treated = sum(dd$TBn == 1),
    n_countries = dplyr::n_distinct(dd$CC))
}

augmented_refits <- function(d0, cfg, extra, by, newvars, tag,
                             extra_note = "") {
  d <- d0 %>%
    inner_join(extra, by = by) %>%
    filter(if_all(all_of(newvars), is.finite))
  n_lost <- nrow(d0) - nrow(d)
  if (nrow(d) < 30 || length(unique(d$TBn)) < 2) {
    say("  [", tag, "] too few complete-case segments (", nrow(d), "); skipped.\n")
    return(NULL)
  }
  sayf("  [%s] complete-case subsample: %d of %d segments (%d treated; %d dropped).\n",
       tag, nrow(d), nrow(d0), sum(d$TBn == 1), n_lost)
  
  run_one <- function(ps_cov, spec_label) {
    dd <- d
    cov_all <- intersect(unique(c(attr(d0, "cov_all"), newvars)), names(dd))
    for (v in cov_all) {
      mu <- wmean(dd[[v]][dd$TBn == 1], rep(1, sum(dd$TBn == 1)))
      dd[[paste0(v, "_ctr")]] <- dd[[v]] - mu
    }
    attr(dd, "ps_cov")  <- intersect(ps_cov, names(dd))
    attr(dd, "cov_all") <- cov_all
    for (a in c("ratio_source", "h2_source")) attr(dd, a) <- attr(d0, a)
    dd <- add_design_weights(dd, cfg)
    if (is.null(dd)) { say("    propensity model failed: ", spec_label, "\n"); return(NULL) }
    o1 <- fit_stage2(dd, "beta_0", "bv_0", cfg, tag = paste(tag, spec_label, "H1"))
    o2 <- fit_stage2(dd, "z",      "bv_z", cfg, tag = paste(tag, spec_label, "H2"))
    s1 <- summarise_fit(o1, cfg, paste(tag, spec_label, "H1"), hyp = "H1")
    s2 <- summarise_fit(o2, cfg, paste(tag, spec_label, "H2"), hyp = "H2")
    out <- dplyr::bind_rows(
      .aug_row(s1, "H1", "mm/yr",    spec_label, dd, cfg),
      .aug_row(s2, "H2", "Fisher z", spec_label, dd, cfg))
    if (is.null(out) || !nrow(out)) return(NULL)
    # balance of the added variable(s) in this design
    for (v in newvars) {
      out[[paste0("smd_", v, "_unweighted")]] <-
        smd(dd[[v]], dd$TBn, rep(1, nrow(dd)))
      out[[paste0("smd_", v, "_ato")]] <- smd(dd[[v]], dd$TBn, dd$w_ato)
    }
    out
  }
  
  ps0 <- attr(d0, "ps_cov")
  res <- dplyr::bind_rows(
    run_one(ps0, "preferred covariates, complete-case subsample"),
    run_one(union(ps0, newvars),
            paste0("augmented propensity + outcome (+ ",
                   paste(newvars, collapse = " + "), ")")))
  if (!is.null(res) && nrow(res)) {
    res$segments_dropped_missing <- n_lost
    res$note <- extra_note
  }
  res
}


# =============================================================================
# PART 1 (A1 / R1.4). WELL-LEVEL FIXED-EFFECTS OLS, CLUSTERED INFERENCE
# =============================================================================
# The mixed-model version (segment + country random intercepts, wild cluster
# bootstrap) already ships in table_si_well_level.csv. This is the reviewer's
# literal alternative: pooled OLS at the well level with country dummies,
# each well carrying its segment design weight divided by the segment's well
# count (so a segment's total influence equals its weight in the preferred
# design), inference clustered on country (CR2, Satterthwaite df) and, as a
# second tier, two-way on country and segment.
say("\n---- PART 1: well-level OLS with country FE (R1.4 literal) ----\n")
wl_fe <- NULL
{
  ps_cov  <- attr(d0, "ps_cov")
  cov_ctr <- paste0(ps_cov, "_ctr")
  keep <- d0 %>% select(unit_id, CC, TBn, w_ato, n_seg = n, all_of(cov_ctr))
  wv <- wells %>% select(unit_id, GWSlp, dist_LB_km) %>%
    inner_join(keep, by = "unit_id") %>%
    group_by(unit_id) %>%
    mutate(dist_z = if (sd(dist_LB_km) > 0)
      (dist_LB_km - mean(dist_LB_km)) / sd(dist_LB_km) else NA_real_,
      gw_z = if (sd(GWSlp) > 0)
        (GWSlp - mean(GWSlp)) / sd(GWSlp) else NA_real_) %>%
    ungroup() %>%
    mutate(w_well = w_ato / n_seg) %>%
    # unit-SD scaling of adjustment covariates (conditioning only; the
    # reported terms TBn, dist_z and TBn:dist_z are left untouched)
    mutate(across(all_of(cov_ctr), function(v) {
      s <- stats::sd(v, na.rm = TRUE); if (is.finite(s) && s > 0) v / s else v
    }))
  
  fit_fe <- function(form, dat, term, lab, unit) {
    yv  <- all.vars(form)[1]
    dat <- dat %>% filter(is.finite(.data[[yv]]), w_well > 0)
    f <- tryCatch(lm(form, data = dat, weights = w_well), error = function(e) NULL)
    if (is.null(f) || !term %in% names(coef(f))) return(NULL)
    est <- unname(coef(f)[term])
    # tier 1: CR2 on country (Satterthwaite)
    se_cr2 <- df_cr2 <- p_cr2 <- NA_real_
    cr2_source <- NA_character_
    if (have_club) {
      ct <- tryCatch(clubSandwich::coef_test(
        f, vcov = "CR2", cluster = dat$CC, coefs = term, test = "Satterthwaite"),
        error = function(e) { say("    [CR2 failed] ", conditionMessage(e), "\n"); NULL })
      if (!is.null(ct)) {
        se_cr2 <- ct$SE[1]; df_cr2 <- ct$df[1]; p_cr2 <- ct$p_Satt[1]
        cr2_source <- "clubSandwich CR2 (Satterthwaite)"
      }
    }
    if (!is.finite(se_cr2) && have_sandwich) {
      # fallback: one-way country clustering (CR1/HC1), t with G-1 df
      V1 <- tryCatch(sandwich::vcovCL(f, cluster = dat$CC, type = "HC1"),
                     error = function(e) { say("    [CL1 failed] ", conditionMessage(e), "\n"); NULL })
      if (!is.null(V1) && term %in% rownames(V1)) {
        se_cr2 <- sqrt(V1[term, term])
        df_cr2 <- dplyr::n_distinct(dat$CC) - 1
        p_cr2  <- 2 * pt(abs(est / se_cr2), df = df_cr2, lower.tail = FALSE)
        cr2_source <- "sandwich CL1 fallback (t, G-1 df)"
      }
    }
    # tier 2: two-way country x segment (CGM), HC1 small-sample
    se_2w <- p_2w <- NA_real_
    if (have_sandwich) {
      V <- tryCatch(sandwich::vcovCL(
        f, cluster = dat[, c("CC", "unit_id")], type = "HC1"),
        error = function(e) NULL)
      if (!is.null(V) && term %in% rownames(V)) {
        se_2w <- sqrt(V[term, term])
        p_2w  <- 2 * pt(abs(est / se_2w), df = dplyr::n_distinct(dat$CC) - 1,
                        lower.tail = FALSE)
      }
    }
    tibble::tibble(
      outcome = lab, term = term, unit = unit, estimate = est,
      se_ols = sqrt(diag(vcov(f)))[term],
      se_cr2_country = se_cr2, df_cr2 = df_cr2, p_cr2_country = p_cr2,
      cr2_source = cr2_source,
      se_twoway_cc_segment = se_2w, p_twoway_cc_segment = p_2w,
      n_wells = nrow(dat), n_segments = dplyr::n_distinct(dat$unit_id),
      n_countries = dplyr::n_distinct(dat$CC),
      note = paste("OLS, country fixed effects, segment-total weights",
                   "(w_ATO/n_k). CR2 country clustering (Satterthwaite df);",
                   "two-way CGM country x segment where available. The",
                   "reviewer's specification; complements the mixed-model",
                   "rows of table_si_well_level.csv."))
  }
  
  f1 <- reformulate(c("TBn", cov_ctr, "factor(CC)"), response = "GWSlp")
  f2 <- reformulate(c("TBn * dist_z", cov_ctr, "factor(CC)"), response = "gw_z")
  wl_fe <- dplyr::bind_rows(
    fit_fe(f1, wv, "TBn", "H1 analogue: well-level depletion", "mm/yr"),
    fit_fe(f2, wv %>% filter(is.finite(dist_z)), "TBn:dist_z",
           "H2 analogue: within-segment standardised gradient",
           "SD depletion per SD distance"))
  if (!is.null(wl_fe) && nrow(wl_fe)) {
    readr::write_csv(wl_fe, file.path(cfg$out_dir, "table_si_well_level_fe.csv"))
    print(as.data.frame(wl_fe[, 1:9]), digits = 3, row.names = FALSE)
  } else say("  fixed-effects fits failed; nothing written.\n")
  if (!have_club)     say("  NOTE: 'clubSandwich' missing -> no CR2 tier.\n")
  if (!have_sandwich) say("  NOTE: 'sandwich' missing -> no two-way tier.\n")
}


# =============================================================================
# PART 2 (A3 / R5.3). GOVERNANCE-PROXY ROBUSTNESS
# =============================================================================
say("\n---- PART 2: governance proxy (R5.3) ----\n")
gov_path <- file.path(cfg$root, "1_data", "governance_proxy.csv")
gov_tab <- NULL
if (!file.exists(gov_path)) {
  tmpl <- d0 %>% group_by(CC) %>%
    summarise(n_segments = dplyr::n(), n_treated = sum(TBn == 1),
              n_wells = sum(n), .groups = "drop") %>%
    arrange(desc(n_wells)) %>%
    mutate(gov_score = NA_real_, source = "", notes = "")
  readr::write_csv(tmpl, gov_path)
  say("  TEMPLATE WRITTEN: ", gov_path, "\n",
      "  Fill 'gov_score' for the ", nrow(tmpl), " analysis countries",
      " (e.g. 0/1 statutory permit regime, or an ordinal coding after",
      " Edwards 2016 / Edwards & Guilfoos 2021), cite each value in",
      " 'source', then rerun this script. Part 2 skipped.\n")
} else {
  gov <- readr::read_csv(gov_path, show_col_types = FALSE)
  if (!all(c("CC", "gov_score") %in% names(gov)))
    stop("governance_proxy.csv must carry columns CC and gov_score.")
  gov <- gov %>% select(CC, gov_score) %>% distinct()
  n_scored <- sum(is.finite(gov$gov_score[gov$CC %in% d0$CC]))
  if (n_scored == 0) {
    say("  governance_proxy.csv present but all gov_score are NA; skipped.\n")
  } else {
    if (n_scored < dplyr::n_distinct(d0$CC))
      say("  NOTE: ", dplyr::n_distinct(d0$CC) - n_scored,
          " of ", dplyr::n_distinct(d0$CC),
          " countries unscored; their segments drop from this check.\n")
    gov_tab <- augmented_refits(
      d0, cfg, gov, by = "CC", newvars = "gov_score", tag = "governance",
      extra_note = paste("Country-level groundwater-governance proxy added to",
                         "propensity and outcome models (R5.3). Preferred row",
                         "is refit on the identical scored subsample."))
    if (!is.null(gov_tab)) {
      readr::write_csv(gov_tab, file.path(cfg$out_dir, "table_si_governance_robustness.csv"))
      show_cols <- intersect(c("specification", "hypothesis", "estimate",
                               "p_one_sided", "mde80", "n_segments",
                               "n_treated", "n_countries"), names(gov_tab))
      print(as.data.frame(gov_tab[, show_cols]), digits = 3, row.names = FALSE)
    }
  }
}


# =============================================================================
# PART 3 (A4 / R1.1). DEPTH-TO-WATER ROBUSTNESS (post-treatment; check only)
# =============================================================================
say("\n---- PART 3: depth to water (A4; post-treatment robustness) ----\n")
wtd_well_path <- file.path(cfg$root, "1_data", "wtd_by_well.csv")
wtd_seg_path <- file.path(cfg$root, "1_data", "wtd_by_segment.csv")
wtd_tif_path <- file.path(cfg$root, "1_data", "wtd_global.tif")
wtd_cache    <- file.path(cfg$cache_dir, "wtd_by_segment_sampled.csv")
wtd_tab <- NULL
wtd <- NULL

if (file.exists(wtd_well_path)) {
  ww <- readr::read_csv(wtd_well_path, show_col_types = FALSE)
  if (!"id" %in% names(ww))
    stop("wtd_by_well.csv must carry the well 'id' column of wellsData.csv.")
  wcols <- grep("^wtd", names(ww), value = TRUE)
  if (!length(wcols))
    stop("wtd_by_well.csv carries no column starting with 'wtd'.")
  if (length(wcols) > 1)
    say("  NOTE: several wtd_* columns found (",
        paste(wcols, collapse = ", "), "); using ", wcols[1],
        ". Reorder columns to choose another.\n")
  ww <- ww %>% transmute(id, .wtd = .data[[wcols[1]]])
  n_matched <- sum(ww$id %in% wells$id)
  wtd <- wells %>% select(id, unit_id) %>%
    inner_join(ww, by = "id") %>%
    filter(is.finite(.wtd)) %>%
    group_by(unit_id) %>%
    summarise(wtd_m = mean(.wtd), n_wells_wtd = dplyr::n(), .groups = "drop")
  sayf("  using well-level depths from %s (column %s): %d wells matched, %d segments covered.\n",
       basename(wtd_well_path), wcols[1], n_matched, nrow(wtd))
} else if (file.exists(wtd_seg_path)) {
  wtd <- readr::read_csv(wtd_seg_path, show_col_types = FALSE)
  if (!all(c("unit_id", "wtd_m") %in% names(wtd)))
    stop("wtd_by_segment.csv must carry columns unit_id and wtd_m.")
  wtd <- wtd %>% select(unit_id, wtd_m) %>% distinct()
  say("  using segment-level depths from ", basename(wtd_seg_path), ".\n")
} else if (file.exists(wtd_cache)) {
  wtd <- readr::read_csv(wtd_cache, show_col_types = FALSE)
  say("  using cached raster sample ", basename(wtd_cache), ".\n")
} else if (file.exists(wtd_tif_path)) {
  if (!have_terra) {
    say("  wtd_global.tif found but package 'terra' is missing; skipped.\n")
  } else {
    say("  sampling ", basename(wtd_tif_path), " at ", nrow(wells), " wells...\n")
    r  <- terra::rast(wtd_tif_path)
    xy <- terra::vect(as.data.frame(wells[, c("lon", "lat")]),
                      geom = c("lon", "lat"), crs = "EPSG:4326")
    v  <- terra::extract(r, xy)[, 2]
    wtd <- wells %>% mutate(.wtd = v) %>%
      filter(is.finite(.wtd)) %>%
      group_by(unit_id) %>%
      summarise(wtd_m = mean(.wtd), n_wells_sampled = dplyr::n(),
                .groups = "drop")
    readr::write_csv(wtd, wtd_cache)
    say("  segment means cached to ", wtd_cache, ".\n")
  }
} else {
  tmpl <- d0 %>% transmute(unit_id, CC, TBn, n_wells = n,
                           lat_centroid = lat_c, lon_centroid = lon_c,
                           wtd_m = NA_real_)
  tp <- file.path(cfg$root, "1_data", "wtd_by_segment_TEMPLATE.csv")
  readr::write_csv(tmpl, tp)
  say("  No depth layer found. Provide either:\n",
      "    1_data/wtd_global.tif  (global water-table depth, m; e.g. Fan et",
      " al. 2013 or Bierkens et al. 2022) -- sampled automatically; or\n",
      "    1_data/wtd_by_segment.csv (unit_id, wtd_m) -- template written to ",
      tp, ".\n  Part 3 skipped.\n")
}

if (!is.null(wtd)) {
  wtd_tab <- augmented_refits(
    d0, cfg, wtd %>% select(unit_id, wtd_m), by = "unit_id",
    newvars = "wtd_m", tag = "wtd",
    extra_note = paste("Depth to water is an outcome of past pumping",
                       "(post-treatment); reported as a robustness check",
                       "only, per the response to R1.1. Preferred row is",
                       "refit on the identical subsample with depth data."))
  if (!is.null(wtd_tab)) {
    readr::write_csv(wtd_tab, file.path(cfg$out_dir, "table_si_wtd_robustness.csv"))
    show_cols <- intersect(c("specification", "hypothesis", "estimate",
                             "p_one_sided", "mde80", "n_segments",
                             "n_treated", "n_countries"), names(wtd_tab))
    print(as.data.frame(wtd_tab[, show_cols]), digits = 3, row.names = FALSE)
  }
}


# ---- consolidate the augmented rows -----------------------------------------
augmented <- if (length(augmented_rows))
  dplyr::bind_rows(augmented_rows) else NULL
if (!is.null(augmented))
  sayf("  augmented specifications: %d rows\n", nrow(augmented))


# ---- persist the consolidated robustness objects ----------------------------
rob <- list(tab = tab, inf_cc = inf_cc, inf_seg = inf_seg, inf_aq = inf_aq,
            locC = locC, h3_bootstrap = h3b, selection_stability = selstab,
            support_window = sw,
            augmented = if (exists("augmented")) augmented else NULL,
            conn_tab  = if (exists("conn_tab"))  conn_tab  else NULL,
            xw        = if (exists("xw"))        xw        else NULL)
cache_write(rob, file.path(cfg$cache_dir, "robustness_objects.rds"),
            cache_stamp(cfg, list(object = "robustness"),
                        extra = list(ref = ref_key)))
sayf("\n34_run_robustness.R done in %.1f min.\n",
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
