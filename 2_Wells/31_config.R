# =============================================================================
# 31_config.R -- CONFIGURATION
#
# Every substantive analytical choice lives here and nowhere else. The defaults
# ARE the preferred specification: run_specification() with no overrides must
# reproduce the frozen benchmark below, and 33_run_main.R asserts that it does.
#
# The pipeline is self-contained on ONE input file, 1_data/wellsData.csv.
# Every first-stage quantity (segment depletion means, depletion-distance
# correlations, segment covariates) is recomputed from the wells; there is no
# second input file to drift out of sync with.
#
# Sourcing this file defines CFG and nothing else. It fits nothing and writes
# nothing.
# =============================================================================

# ---- project root -----------------------------------------------------------
# All paths are resolved against the directory containing this file, so the
# pipeline runs from any working directory ("here"-style without the package).
.this_file <- (function() {
  # When source()d, the calling frame records the file path:
  for (fr in rev(sys.frames())) {
    f <- tryCatch(get("ofile", envir = fr), error = function(e) NULL)
    if (is.character(f) && grepl("31_config[.]R$", f))
      return(normalizePath(f))
  }
  # When run directly via Rscript:
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) && grepl("31_config[.]R$", f[1])) return(normalizePath(f[1]))
  # Fallback: the working directory must contain this file.
  if (file.exists("31_config.R")) return(normalizePath("31_config.R"))
  stop("Cannot locate 31_config.R; run from the repository root.")
})()
PROJECT_ROOT <- dirname(.this_file)

CFG <- list(

  ## ---- paths (all relative to PROJECT_ROOT) --------------------------------
  root         = PROJECT_ROOT,
  wells_file   = file.path(PROJECT_ROOT, "1_data", "wellsData.csv"),
  out_dir      = file.path(PROJECT_ROOT, "output"),          # canonical outputs ONLY
  cache_dir    = file.path(PROJECT_ROOT, "derived", "cache"),
  audit_dir    = file.path(PROJECT_ROOT, "derived", "audit"),
  log_file     = file.path(PROJECT_ROOT, "derived", "audit", "run_log.txt"),

  ## ---- sample --------------------------------------------------------------
  drop_cc      = c("AUS", "NZL", "TWN", "JPN", "IRL", "GBR"),  # no land border
  min_wells    = 20,
  cov_match    = c("lat_c", "lon_c", "urbkHaKm2", "CS_max", "LB_river"),
  cov_extra    = "prec_mm",         # augmented propensity set = cov_match + this
  drop_countries        = NULL,     # leave-one-country-out
  drop_treated_segments = NULL,     # leave-one-treated-segment-out
  drop_aquifers         = NULL,     # leave-one-physical-aquifer-out
  exclude_river_borders = FALSE,

  ## ---- distance thresholds (km) -- centralised, used everywhere ------------
  # 100 km is the theory- and hydrogeology-motivated limit of plausible
  # cross-border hydraulic influence (SI); it is NOT estimated from the data.
  influence_zone_km = 100,
  coverage_km       = c(10, 50, 100),        # within-distance coverage summaries
  profile_bins      = c(0, 10, 25, 50, 100, 200, Inf),  # Panel B bins
  profile_boot      = 999,                   # country-cluster bootstrap reps
  gradient_windows  = c(25, 50, 100, 200, 500, Inf),    # cumulative windows (SI)
  gradient_min_wells = 10,

  ## ---- first stage ---------------------------------------------------------
  first_stage_h1     = "mean",       # mean | robust_mean | trimmed_mean
  first_stage_h2     = "pearson",    # pearson | spearman | physical_slope | distance_bins
  trim_frac          = 0.10,
  distance_window_km = NULL,         # c(lo, hi) restriction on dist_LB_km
  distance_bins      = c(10, 100),   # near (<=10) vs interior (10-100) contrast
  min_wells_per_bin  = 5,
  min_sd_dist_km     = 0,            # stricter distance-variation screen

  ## ---- first-stage spatial (Conley) uncertainty ----------------------------
  conley_km    = 10,
  conley_seed  = 1,
  n_max_wells  = 2500,
  n_eff_power  = 2,                  # n_eff = n_w / ratio_1^power
  ratios_follow_sample = TRUE,

  ## ---- design ---------------------------------------------------------------
  estimand      = "ATO",             # ATO (overlap) | ATT (full matching)
  ps_covariates = NULL,              # NULL -> cov_match

  ## ---- second stage ----------------------------------------------------------
  outcome_adjustment = TRUE,
  country_effect     = "random",     # random | fixed
  include_segment_re = TRUE,
  precision_weight   = TRUE,
  variance_winsor    = c(0.01, 0.99),
  w_winsor           = c(0.01, 0.99),
  max_iter           = 25,
  tol                = 1e-6,

  ## ---- inference: significance levels and directional alternatives ----------
  # PRE-STATED DIRECTIONS (committed in writing before this analysis):
  #   H1: transboundary depletion is MORE POSITIVE than the matched domestic
  #       baseline (alternative "greater").
  #   H2: the transboundary gradient is MORE NEGATIVE (more borderward) than
  #       the matched domestic baseline (alternative "less").
  #   H2 intercept (placebo-type diagnostic): TWO-SIDED, no direction.
  # The primary directional test is one-sided at alpha_one = 0.05, whose
  # critical value equals the near-zero end of the two-sided 90% interval.
  directional_alt = c(H1 = "greater", H2 = "less"),
  alpha_one    = 0.05,               # one-sided directional test level
  ci_level     = 0.90,               # two-sided transparency interval
  n_wcb        = 999,
  wcb_seed     = 7,
  wcb_weights  = "rademacher",
  wcb_ci       = FALSE,              # test-inversion bootstrap CI (main spec only)
  # MDE: minimum detectable effect for the H1 one-sided test at alpha_one,
  # at the powers below. mde = (qnorm(1 - alpha_one) + qnorm(power)) * SE.
  # Verified by an assertion in 32_core.R (assert_mde()).
  power_mde    = c(0.80, 0.90),

  ## ---- H3 -------------------------------------------------------------------
  h3_by            = "unit_id",
  h3_centre        = "domestic",
  h3_variance      = "preferred",    # preferred | precision | legacy
  h3_deviations    = "total",
  q_cut            = 0.10,
  top_k            = 10,
  h3_boot_reps     = 200,            # rank-stability bootstrap replicates
  h3_boot_seed     = 20240601,
  h3_boot_level    = "country",
  h3_boot_resample_wells = TRUE,
  # The bootstrap refits the shrinkage under the (non-iterated) precision
  # variance path -- one of the two voting paths -- for computational
  # tractability; the SI states this scope. BLUP validation is skipped inside
  # replicates (it is validated on every non-bootstrap fit).
  h3_boot_variance = "precision",
  h3_blup_validate = TRUE,
  h3_blup_tol      = 1e-6,
  # H3 SELECTION RULE (coded; the phrase "selected under every analysis
  # considered defensible" refers to exactly this):
  #   selected <=> rank_raw <= top_k on the RAW (unmodelled) combined deviation
  #            AND consistent under EVERY voting variance path, where
  #   consistent <=> top-k on level OR gradient OR combined deviation, or in
  #                  the fast+borderward quadrant, under every voting path.
  h3_paths_voting = c("preferred", "precision"),
  h3_paths_shown  = c("preferred", "precision", "legacy"),
  # Panel C fill classes for the share of wells within influence_zone_km:
  share100_breaks = c(0, 0.25, 0.75, 1),

  ## ---- balance thresholds ----------------------------------------------------
  balance_tol_numeric     = 1e-6,
  balance_tol_substantive = 0.10,

  ## ---- parallel execution ----------------------------------------------------
  # Reproducibility never depends on these: all random loops draw from
  # pre-generated L'Ecuyer-CMRG streams, one per job.
  parallel  = TRUE,
  workers   = NULL,                  # NULL -> detectCores() - 0 (capped at 8)

  ## ---- cache integrity -------------------------------------------------------
  cache_check  = TRUE,
  source_files = file.path(PROJECT_ROOT, c("31_config.R", "32_core.R")),

  ## ---- frozen preferred benchmark -------------------------------------------
  # Locked after the first validated run of 33_run_main.R (see README).
  # A later run that fails to reproduce these numbers STOPS the pipeline.
  reference_h1          = -138.7157244,
  reference_h2          = -0.07896063417,
  reference_h1_ci       = c(-505.5585075, 228.1270586),
  reference_h2_ci       = c(-0.160000922, 0.002079653624),
  reference_n_segments  = 698L,
  reference_n_treated   = 44L,
  reference_n_countries = 19L,
  allow_unlocked_reference = FALSE,
  reproduce_tol    = 1e-8,
  reproduce_tol_ci = 1e-6
)
