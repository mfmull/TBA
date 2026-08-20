# =============================================================================
# 41_config.R -- CONFIGURATION (irrigation / land-use analysis, 41-series)
#
# Mirror of the wells pipeline (31-series) for the irrigated-cropland analysis.
# Every substantive analytical choice lives here and nowhere else. The defaults
# ARE the preferred specification: run_spec() with no overrides must reproduce
# the frozen benchmark below, and 43_run_main.R asserts that it does.
#
# The pipeline is self-contained on THREE input files under 1_data/:
#   _dataMain.csv            polygon-level master table (treat = 929 TBA-country
#                            segments; controls = HydroSHEDS basins), rebuilt
#                            from the current GEE exports (Feb 2026 build).
#   CtrlNoOverlapHYBAS_B.rds 500 non-overlapping control-basin realizations
#                            (greedy random selection; HYBAS_ID vectors).
#   DataS2_meta.csv          IGRAC metadata (aquifer names, codes, countries,
#                            region) keyed on aq_id x CountryName; used ONLY to
#                            regenerate Data_S2 and to label Lorenz examples --
#                            never as an analysis input.
# Plus, for the Fig. 3C Lorenz illustration only: landpct.csv, irpct.csv
# (distance-to-border percentile curves for the treated segments).
#
# CHANGES vs. the legacy irrigation pipeline (see CHANGELOG.md for the full log):
#   (i)   best ensemble member selected on lowest MAX standardized mean
#         difference (as stated in Methods), not lowest mean SMD;
#   (ii)  the propensity-score "distance" row is EXCLUDED from all SMD
#         balance statistics;
#   (iii) `RiverBorder` (which was TRUE when a polygon was NOT within 50 km of
#         a river border) is renamed near_river_border = !RiverBorder at load;
#   (iv)  per-covariate balance (pre/post) is stored for every realization,
#         feeding the SI balance diagnostics (reviewer R3.2);
#   (v)   country-clustered CR2 inference is reported for the best-balanced
#         member of every specification, mirroring the wells pipeline;
#   (vi)  a river-border exclusion robustness block (reviewer R5.4) and the
#         GW-attribution sensitivity (gini_gwpct_gt0; reviewer R5.5) are
#         canonical outputs;
#   (vii) The segment-level supplementary dataset is regenerated from the
#         CURRENT master table (the shipped Data_S2 predated the final GEE
#         rebuild; GW>0 counts changed 207->159).
#  (viii) 2026-08-19: that file becomes SI Dataset S3 (output/Dataset_S3.csv)
#         and is RESTRICTED to the variables the irrigation analysis actually
#         uses -- every outcome, sample filter, matching covariate and
#         regression covariate named in `specs` below, plus the identifiers
#         and the country grouping factor of the random intercept. Six
#         columns of the legacy Data_S2 are dropped:
#           type                 constant ("treat") -- the file is treated-only
#           IrNeed0              no specification uses it
#           gini_crpGWSpct_gt0   no specification uses it, AND it is an EXACT
#                                duplicate of G_IrrigNeed (verified: max abs
#                                difference 0 over all 929 rows)
#           G_IrrigNeed          descriptive only; no specification uses it
#           area_km2             descriptive only; the control-size caliper is
#                                applied upstream, in the control construction
#           Countries            list of riparians; already in Dataset S2,
#                                keyed on the same aq_id
#         The two descriptive drops (area_km2, G_IrrigNeed) remove their rows
#         from Table S (summary statistics); 47_tables.R is updated to match.
#
# Sourcing this file defines CFG and nothing else. It fits nothing and writes
# nothing.
# =============================================================================

# ---- project root -----------------------------------------------------------
.this_file <- (function() {
  for (fr in rev(sys.frames())) {
    f <- tryCatch(get("ofile", envir = fr), error = function(e) NULL)
    if (is.character(f) && grepl("41_config[.]R$", f))
      return(normalizePath(f))
  }
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) && grepl("41_config[.]R$", f[1])) return(normalizePath(f[1]))
  if (file.exists("41_config.R")) return(normalizePath("41_config.R"))
  stop("Cannot locate 41_config.R; run from the repository root.")
})()
PROJECT_ROOT <- dirname(.this_file)

CFG <- list(

  ## ---- paths (all relative to PROJECT_ROOT) --------------------------------
  root          = PROJECT_ROOT,
  data_file     = file.path(PROJECT_ROOT, "1_data", "_dataMain.csv"),
  ctrl_sets_file= file.path(PROJECT_ROOT, "1_data", "CtrlNoOverlapHYBAS_B.rds"),
  meta_file     = file.path(PROJECT_ROOT, "1_data", "DataS2_meta.csv"),
  landpct_file  = file.path(PROJECT_ROOT, "1_data", "landpct.csv"),
  irpct_file    = file.path(PROJECT_ROOT, "1_data", "irpct.csv"),
  out_dir       = file.path(PROJECT_ROOT, "output"),   # canonical tables/csv
  fig_dir       = file.path(PROJECT_ROOT, "figure"),   # canonical figures
  cache_dir     = file.path(PROJECT_ROOT, "derived", "cache"),
  audit_dir     = file.path(PROJECT_ROOT, "derived", "audit"),
  log_file      = file.path(PROJECT_ROOT, "derived", "audit", "run_log.txt"),

  ## ---- ensemble ------------------------------------------------------------
  n_realizations = 500,        # non-overlapping control realizations used
  sig_level      = 0.10,       # ensemble significance threshold (as in v3)
  # furrr derives its per-realization L'Ecuyer streams from the RNG state at
  # the call, so the state is fixed here. Nothing in fit_one() currently draws
  # -- full optimal matching and lmer are deterministic and the control
  # realizations are pre-generated -- so this changes no number today; it means
  # that adding a stochastic step later cannot break reproducibility silently.
  ensemble_seed  = 20260806L,

  ## ---- design ---------------------------------------------------------------
  # Full optimal matching (MatchIt::matchit, method = "full"), default
  # estimand: the AVERAGE TREATMENT EFFECT ON THE TREATED (ATT) -- stated
  # explicitly in the revised Methods (legacy text was silent).
  match_method  = "full",
  cov_int       = c("lat_c", "lon_c", "RDS_mainDist", "CS_max"),
  cov_gini      = c("lat_c", "lon_c", "RDS_mainDist", "giniCSI"),
  cov_geo       = c("lat_c", "lon_c", "RDS_mainDist"),
  river_cov     = "near_river_border",   # TRUE = within 50 km of a GSRB river border

  ## ---- balance / best-member selection --------------------------------------
  # FIX (i)+(ii): Methods and Table S5 caption say "lowest maximum standardized
  # mean difference"; the legacy code sorted on the MEAN SMD and included the
  # propensity 'distance' row. The defaults below implement the stated rule.
  best_rule             = "max_smd",   # max_smd (stated rule) | mean_smd (legacy)
  smd_exclude_distance  = TRUE,        # drop PS 'distance' row from SMD stats
  balance_tol_substantive = 0.10,      # |SMD| threshold flagged in balance table

  ## ---- second stage ----------------------------------------------------------
  # Weighted linear mixed model with a country random intercept, full-matching
  # weights carried into the regression (as in v3). Ensemble p-values are
  # Satterthwaite (lmerTest). For the best-balanced member of each
  # specification we ADDITIONALLY report CR2 country-clustered SEs
  # (clubSandwich), mirroring the wells pipeline's robust-inference tier.
  ci_level      = 0.90,

  ## ---- specifications --------------------------------------------------------
  # filter: expression on the master table (applied to treated AND controls);
  # group: figure grouping (main_intensity / main_gini / si / robust).
  specs = list(
    IrIntens     = list(outcome = "Ir",       filter = "Ir > 0",
                        cov_match = "cov_int",  cov_reg = NULL,
                        group = "main_intensity", title = "Irrigation"),
    Overdraft    = list(outcome = "OverIR3",  filter = "Ir > 0",
                        cov_match = "cov_int",  cov_reg = NULL,
                        group = "main_intensity", title = "Overdraft"),
    IrRivsInt    = list(outcome = "Ir",       filter = "Ir > 0",
                        cov_match = "cov_int",  cov_reg = "river_cov",
                        group = "main_intensity", title = "Irrig | Rivers"),
    IrGini       = list(outcome = "giniIr",   filter = "Ir > 0",
                        cov_match = "cov_gini", cov_reg = NULL,
                        group = "main_gini", title = "Irrigation"),
    IrRivsGini   = list(outcome = "giniIr",   filter = "Ir > 0",
                        cov_match = "cov_gini", cov_reg = "river_cov",
                        group = "main_gini", title = "Irrig | Rivers"),
    GWGini       = list(outcome = "giniGw",   filter = "GW > 0",
                        cov_match = "cov_gini", cov_reg = NULL,
                        group = "main_gini", title = "GW Irrig."),
    GWRivsGini   = list(outcome = "giniGw",   filter = "GW > 0",
                        cov_match = "cov_gini", cov_reg = "river_cov",
                        group = "main_gini", title = "GW Irrig. | Rivers"),
    ## -- SI --
    CropSuitInt  = list(outcome = "CS_max",   filter = "Ir > 0",
                        cov_match = "cov_geo",  cov_reg = NULL,
                        group = "si", title = "Crop Suitability"),
    IrNeedInt    = list(outcome = "IrNeed3",  filter = "Ir > 0",
                        cov_match = "cov_geo",  cov_reg = NULL,
                        group = "si", title = "Irrigation Need"),
    Overdraft_Irrig = list(outcome = "OverIR3", filter = "Ir > 0",
                        cov_match = c("lat_c", "lon_c", "RDS_mainDist", "Ir"),
                        cov_reg = NULL,
                        group = "si", title = "Overdraft | Irrig"),
    Overdraft6   = list(outcome = "OverIR6",  filter = "Ir > 0",
                        cov_match = "cov_int",  cov_reg = NULL,
                        group = "si", title = "Overdraft (>6 mo)"),
    IrIntDoubleRobust = list(outcome = "Ir",  filter = "Ir > 0",
                        cov_match = "cov_int",
                        cov_reg = c("RDS_mainDist", "CS_max"),
                        group = "si", title = "Irrigation (double robust)"),
    IrGiniDoubleRobust = list(outcome = "giniIr", filter = "Ir > 0",
                        cov_match = "cov_gini",
                        cov_reg = c("RDS_mainDist", "giniCSI"),
                        group = "si", title = "Irrig. Gini (double robust)"),
    GWShareGt0   = list(outcome = "gini_gwpct_gt0", filter = "GW > 0",
                        cov_match = "cov_gini", cov_reg = NULL,
                        group = "si", title = "GW Irrig. (GWf > 0)"),
    ## -- river-border exclusion robustness (reviewer R5.4) --
    IrIntens_exRiv = list(outcome = "Ir",     filter = "Ir > 0 & !near_river_border",
                        cov_match = "cov_int",  cov_reg = NULL,
                        group = "robust", title = "Irrigation (excl. river borders)"),
    IrGini_exRiv   = list(outcome = "giniIr", filter = "Ir > 0 & !near_river_border",
                        cov_match = "cov_gini", cov_reg = NULL,
                        group = "robust", title = "Irrig. Gini (excl. river borders)"),
    GWGini_exRiv   = list(outcome = "giniGw", filter = "GW > 0 & !near_river_border",
                        cov_match = "cov_gini", cov_reg = NULL,
                        group = "robust", title = "GW Gini (excl. river borders)")
  ),

  ## ---- SI supplementary dataset (Dataset S3) --------------------------------
  # FIX (viii). The published segment-level dataset. `si_dataset_cols` is the
  # WHOLE contract: build_si_irrigation_dataset() emits exactly these columns
  # in exactly this order, and 48_run_all.R asserts it. Every entry is
  # justified by the `specs` block above -- if you add a specification that
  # needs a new variable, add it here too or the assertion will fire.
  si_dataset_name = "Dataset_S3.csv",
  si_dataset_cols = c(
    # -- identifiers (also: CountryName is the random-intercept grouping) -----
    "aq_id", "IGRAC_Code", "AquiferName", "CountryName", "CountryCode", "Region",
    # -- matching / regression covariates -------------------------------------
    "lat_c",           # cov_int, cov_gini, cov_geo
    "lon_c",           # cov_int, cov_gini, cov_geo
    "RDS_mainDist",    # cov_int, cov_gini, cov_geo; cov_reg (double robust)
    "CSI",             # = CS_max: cov_int; cov_reg (IrIntDoubleRobust);
                       #   outcome of CropSuitInt
    "G_CSI",           # = giniCSI: cov_gini; cov_reg (IrGiniDoubleRobust)
    "NearRiverBorder", # = near_river_border: river_cov; exRiv sample filter
    # -- outcomes and sample filters ------------------------------------------
    "Irrig",           # = Ir: outcome (IrIntens, IrRivsInt, IrIntDoubleRobust,
                       #   IrIntens_exRiv); "Ir > 0" filter; cov_match of
                       #   Overdraft_Irrig
    "GWIrrig",         # = GW: "GW > 0" filter (GWGini, GWRivsGini,
                       #   GWShareGt0, GWGini_exRiv)
    "IrrigNeed",       # = IrNeed3: outcome of IrNeedInt
    "Overdraft",       # = OverIR3: outcome of Overdraft, Overdraft_Irrig
    "OverIR6",         # outcome of Overdraft6
    "G_Irrig",         # = giniIr: outcome (IrGini, IrRivsGini,
                       #   IrGiniDoubleRobust, IrGini_exRiv)
    "G_IrrigGW",       # = giniGw: outcome (GWGini, GWRivsGini, GWGini_exRiv)
    "gini_gwpct_gt0"), # outcome of GWShareGt0

  ## ---- Lorenz illustration (Fig. 3C) ----------------------------------------
  # Exact aq_id x country segments; chosen (as in v3) to span border-proximate,
  # near-uniform, and interior concentration. Gini values are recomputed from
  # the current percentile curves, not hard-coded.
  lorenz_examples = list(
    list(aq_id = 286, country = "Saudi Arabia"),  # Saq-Ram / Disi   (border-proximate)
    list(aq_id =  63, country = "India"),         # Indo-Gangetic    (near-uniform)
    list(aq_id = 314, country = "Kazakhstan")),   # Syr Darya        (interior)

  ## ---- parallel execution ----------------------------------------------------
  parallel  = TRUE,
  workers   = NULL,                  # NULL -> detectCores() (capped at 8)

  ## ---- cache integrity -------------------------------------------------------
  cache_check  = TRUE,
  source_files = file.path(PROJECT_ROOT, c("41_config.R", "42_core.R")),

  ## ---- frozen preferred benchmark -------------------------------------------
  # Locked after the first validated run of 43_run_main.R. A later run that
  # fails to reproduce these numbers STOPS the pipeline. (Ensemble-mean treated
  # contrast and share of realizations significant at sig_level, for the two
  # headline specifications.)
  reference_IrIntens_meaneff  = 0.01273323725,   # locked 2026-08-06 run
  reference_IrIntens_sharesig = 0.996,
  reference_IrGini_meaneff    = -0.06395431676,
  reference_IrGini_sharesig   = 0.654,
  reference_n_treat_ir        = 562L,   # model sample, Ir>0 specifications
  reference_n_treat_gw        = 159L,   # DATA subset with GW>0 (Dataset S3 check)
  reference_n_si_dataset_rows = 929L,   # treated segments in Dataset_S3.csv
  reference_n_treat_gw_model  = 72L,    # GWGini model sample after complete
                                        # cases on giniGw + matching covariates
                                        # (v3 Table S5: Num.Obs 191 = 72 + 119)
  allow_unlocked_reference    = TRUE,
  reproduce_tol               = 1e-8
)
