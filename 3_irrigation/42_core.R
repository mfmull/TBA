# =============================================================================
# 42_core.R -- CORE LIBRARY (irrigation / land-use analysis)
#
# Data loading, ensemble matching + outcome models, balance diagnostics,
# best-member refit with CR2 country-clustered inference, caching, logging.
# Sourcing this file defines functions only; it fits nothing and writes
# nothing until run_spec() is called.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr);  library(tidyr);   library(purrr);  library(tibble)
  library(readr);  library(stringr); library(MatchIt)
  library(lme4);   library(lmerTest)
  library(broom.mixed)
  library(clubSandwich)
  library(future); library(furrr)
})

# =============================================================================
# SMALL UTILITIES + LOGGING
# =============================================================================
cfg_with <- function(cfg = CFG, ...) {
  ov <- list(...)
  for (nm in names(ov)) cfg[[nm]] <- ov[[nm]]
  cfg
}

.log <- new.env(); .log$path <- NULL
log_open <- function(path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  .log$path <- path
}
say <- function(...) {
  msg <- paste0(...)
  cat(msg)
  if (!is.null(.log$path)) cat(msg, file = .log$path, append = TRUE)
  invisible(NULL)
}
sayf <- function(...) say(sprintf(...))

pstar <- function(p) {
  ifelse(is.na(p), "NA",
  ifelse(p < 0.001, "<0.001", formatC(p, digits = 3, format = "f")))
}
ess <- function(w) (sum(w)^2) / sum(w^2)

resolve_cov <- function(x, cfg) {
  # A cov_match / cov_reg entry may name a CFG covariate set ("cov_int",
  # "cov_gini", "cov_geo", "river_cov") or be a literal character vector.
  if (is.null(x)) return(NULL)
  out <- unlist(lapply(x, function(v)
    if (length(v) == 1 && v %in% names(cfg)) cfg[[v]] else v))
  unname(out)
}

# Run `expr` under a fixed RNG seed and restore the caller's stream afterwards,
# so seeding the ensemble cannot perturb anything downstream of it.
with_seed <- function(seed, expr, kind = "Mersenne-Twister") {
  if (is.null(seed)) seed <- 20260806L
  had  <- exists(".Random.seed", envir = globalenv())
  old  <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  okind <- RNGkind()
  on.exit({
    suppressWarnings(RNGkind(okind[1], okind[2], okind[3]))
    if (had) assign(".Random.seed", old, envir = globalenv())
  }, add = TRUE)
  set.seed(seed, kind = kind)
  force(expr)
}

par_workers <- function(cfg = CFG) {
  if (!isTRUE(cfg$parallel)) return(1L)
  w <- cfg$workers
  if (is.null(w)) w <- max(1L, parallel::detectCores())
  min(as.integer(w), 8L)
}

# =============================================================================
# CACHING (content-stamped, template-style)
# =============================================================================
.file_stamp <- function(paths) {
  paths <- paths[file.exists(paths)]
  paste(unname(tools::md5sum(paths)), collapse = "|")
}
# Hash the OBJECT, never a printed rendering of it. The previous digest was
# paste(capture.output(str(key))), and str() is a display function: it
# truncates character vectors after four elements, rounds numerics to three
# significant digits, and cuts lists at 99 entries. Verified against that
# code, all three of these collided -- identical digests for cov_match vectors
# differing in the 5th element, for 0.10 vs 0.1004, and for 120-element lists
# differing at element 120. Serialising and md5-summing has no such blind spot.
.md5_of <- function(obj) {
  f <- tempfile(); on.exit(unlink(f), add = TRUE)
  con <- file(f, "wb"); serialize(obj, con, ascii = FALSE); close(con)
  unname(tools::md5sum(f))
}
cache_stamp <- function(cfg = CFG, spec = NULL, overrides = list()) {
  # The code stamp covers 42_core.R only, so that editing 41_config.R comments
  # or locking the frozen benchmark does not invalidate caches. Everything in
  # 41_config.R that can change a NUMBER must therefore be named explicitly
  # here -- and that includes the covariates.
  #
  # A spec stores an ALIAS ("cov_match = 'cov_int'"), which resolve_cov()
  # expands against cfg at run time. Stamping the raw spec therefore recorded
  # the alias and not its contents: editing cov_int in 41_config.R changed
  # what every intensity specification matched on while leaving the stamp
  # identical, so the cache served results computed under the old covariate
  # set. Resolve before stamping.
  if (!is.null(spec) && is.list(spec)) {
    if (!is.null(spec$cov_match)) spec$cov_match <- resolve_cov(spec$cov_match, cfg)
    if (!is.null(spec$cov_reg))   spec$cov_reg   <- resolve_cov(spec$cov_reg,   cfg)
  }
  key <- list(
    data   = .file_stamp(c(cfg$data_file, cfg$ctrl_sets_file)),
    code   = .file_stamp(file.path(cfg$root, "42_core.R")),
    n      = cfg$n_realizations,
    rule   = cfg$best_rule,
    xdist  = cfg$smd_exclude_distance,
    # sig_level enters counts$share_sig -- the ensemble percentages reported in
    # Fig. 3 -- and match_method defines the design itself. Both were absent.
    alpha  = cfg$sig_level,
    method = cfg$match_method,
    spec   = spec,
    ov     = overrides)
  key$digest <- .md5_of(key)
  key
}
stamp_matches <- function(a, b) identical(a$digest, b$digest)
cache_read <- function(path, stamp, cfg = CFG) {
  if (!file.exists(path)) return(NULL)
  # A truncated or corrupt .rds (interrupted write, cloud-sync conflict) should
  # be discarded and recomputed, not halt the pipeline.
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(obj)) {
    say("  [cache] unreadable: ", basename(path), " -- recomputing\n")
    return(NULL)
  }
  if (isTRUE(cfg$cache_check) && !stamp_matches(obj$stamp, stamp)) {
    say("  [cache] stale stamp for ", basename(path), " -- recomputing\n")
    return(NULL)
  }
  obj$value
}
cache_write <- function(value, path, stamp) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(stamp = stamp, value = value), path)
  invisible(value)
}

# =============================================================================
# DATA LOADING
# =============================================================================
load_main_data <- function(cfg = CFG) {
  d <- read.csv(cfg$data_file) %>% as_tibble()
  # FIX (iii): the legacy column `RiverBorder` was TRUE when a polygon was NOT
  # within 50 km of a river-defined international border (built as
  # `!aq_id %in% aq_idsRivs` in 2_buildDataset.R). Rename to an unambiguous
  # flag; TRUE = within 50 km of a GSRB river border.
  stopifnot("RiverBorder" %in% names(d))
  d %>% mutate(near_river_border = !RiverBorder) %>% select(-RiverBorder)
}
load_control_sets <- function(cfg = CFG) {
  sets <- readRDS(cfg$ctrl_sets_file)
  stopifnot(length(sets) >= cfg$n_realizations)
  sets[seq_len(cfg$n_realizations)]
}
load_meta <- function(cfg = CFG) read.csv(cfg$meta_file) %>% as_tibble()

# =============================================================================
# BALANCE EXTRACTION
# =============================================================================
# Per-covariate standardized mean differences before/after matching, from a
# matchit object. FIX (ii): the propensity 'distance' row is excluded from the
# scalar summaries when cfg$smd_exclude_distance (it is retained, labelled, in
# the per-covariate table so the SI can show it transparently).
balance_of <- function(m.out, cfg = CFG) {
  s <- summary(m.out, standardize = TRUE)
  pre  <- s$sum.all[,     "Std. Mean Diff.", drop = TRUE]
  post <- s$sum.matched[, "Std. Mean Diff.", drop = TRUE]
  tab <- tibble(covariate = rownames(s$sum.matched),
                smd_pre  = unname(pre[rownames(s$sum.matched)]),
                smd_post = unname(post))
  keep <- if (isTRUE(cfg$smd_exclude_distance))
    tab$covariate != "distance" else rep(TRUE, nrow(tab))
  list(table    = tab,
       mean_smd = mean(abs(tab$smd_post[keep]), na.rm = TRUE),
       max_smd  = max(abs(tab$smd_post[keep]),  na.rm = TRUE))
}

# =============================================================================
# ONE REALIZATION: match + fit
# =============================================================================
fit_one <- function(trdat, cndat, control_ids, spec, cfg = CFG,
                    keep_fit = FALSE) {
  tryCatch({
    cov_match <- resolve_cov(spec$cov_match, cfg)
    cov_reg   <- resolve_cov(spec$cov_reg,   cfg)

    dat <- bind_rows(trdat, cndat %>% filter(aq_id %in% control_ids)) %>%
      mutate(type = factor(type))

    form_match <- reformulate(cov_match, response = "type")
    m.out <- matchit(form_match, data = dat, method = cfg$match_method)
    bal   <- balance_of(m.out, cfg)

    mdat <- match.data(m.out) %>% mutate(weights_combined = weights)

    # same-country share among matched controls (design diagnostic)
    sub_ctry <- mdat %>% filter(type == "treat") %>% group_by(subclass) %>%
      summarise(treat_country = first(CntrName), .groups = "drop")
    ctrl <- mdat %>% filter(type != "treat") %>%
      left_join(sub_ctry, by = "subclass") %>%
      mutate(same_country = CntrName == treat_country)
    sc_unw <- mean(ctrl$same_country, na.rm = TRUE)
    sc_w   <- weighted.mean(ctrl$same_country, w = ctrl$weights_combined,
                            na.rm = TRUE)

    n_treat   <- sum(mdat$type == "treat")
    n_control <- sum(mdat$type != "treat")
    ess_ratio <- ess(mdat$weights_combined[mdat$type != "treat"]) /
                 ess(mdat$weights_combined[mdat$type == "treat"])

    form_lm <- reformulate(c("type", cov_reg), response = spec$outcome)
    fit <- lmerTest::lmer(update(form_lm, . ~ . + (1 | CntrName)),
                          data = mdat, weights = weights_combined)
    tf  <- broom.mixed::tidy(fit, effects = "fixed")
    trow <- tf %>% filter(grepl("^type", term))
    irow <- tf %>% filter(term == "(Intercept)")

    out <- list(
      success   = TRUE,
      mean_smd  = bal$mean_smd, max_smd = bal$max_smd,
      balance   = bal$table,
      ess_ratio = ess_ratio,
      treat_eff = trow$estimate, treat_se = trow$std.error,
      treat_p   = trow$p.value,
      int_eff   = irow$estimate, int_se = irow$std.error,
      int_p     = irow$p.value,
      n_treat   = n_treat, n_control = n_control,
      n_total   = n_treat + n_control,
      same_country_unweighted = sc_unw,
      same_country_weighted   = sc_w)
    if (keep_fit) { out$fit <- fit; out$match <- m.out; out$mdat <- mdat }
    out
  }, error = function(e) list(success = FALSE, error = conditionMessage(e)))
}

# =============================================================================
# ONE SPECIFICATION ACROSS THE ENSEMBLE
# =============================================================================
prepare_spec_data <- function(spec, d, cfg = CFG) {
  cov_match <- resolve_cov(spec$cov_match, cfg)
  cov_reg   <- resolve_cov(spec$cov_reg,   cfg)
  dd <- d %>% filter(!!rlang::parse_expr(spec$filter))
  vars_keep <- unique(c(spec$outcome, "type", "aq_id", "CntrName",
                        cov_match, cov_reg))
  vars_keep <- vars_keep[vars_keep %in% names(dd)]
  dd <- dd %>% select(all_of(vars_keep)) %>%
    filter(complete.cases(across(everything())))
  list(trdat = dd %>% filter(type == "treat"),
       cndat = dd %>% filter(type != "treat"))
}

summarise_ensemble <- function(results, cfg = CFG) {
  summary_df <- imap_dfr(results, function(x, i) {
    if (isTRUE(x$success)) tibble(
      idx = as.integer(i),
      mean_smd = x$mean_smd, max_smd = x$max_smd, ess_ratio = x$ess_ratio,
      treat_eff = x$treat_eff, treat_se = x$treat_se, treat_p = x$treat_p,
      int_eff = x$int_eff, int_se = x$int_se, int_p = x$int_p,
      n_treat = x$n_treat, n_control = x$n_control, n_total = x$n_total,
      same_country_unweighted = x$same_country_unweighted,
      same_country_weighted   = x$same_country_weighted)
    else tibble(idx = as.integer(i), mean_smd = NA_real_)
  }) %>% filter(!is.na(mean_smd))

  # Every share below is computed over the realizations that CONVERGED. Record
  # the denominator explicitly: the assertions tolerate up to 5% loss, so
  # "100% of realizations significant" could mean 100% of 476 of 500, and
  # without these fields nothing downstream could tell the difference.
  n_requested <- length(results)
  n_used      <- nrow(summary_df)
  if (n_used < n_requested)
    sayf("    %d of %d realizations failed and are excluded from the shares.\n",
         n_requested - n_used, n_requested)

  a <- cfg$sig_level
  counts <- summary_df %>%
    summarise(
      total = n(),
      n_requested = n_requested,
      n_failed    = n_requested - n_used,
      n_sig       = sum(treat_p < a),
      n_pos_sig   = sum(treat_p < a & treat_eff > 0),
      n_neg_sig   = sum(treat_p < a & treat_eff < 0),
      share_sig     = n_sig / total,
      share_pos_sig = n_pos_sig / total,
      share_neg_sig = n_neg_sig / total,
      n_int_sig     = sum(int_p < a),
      n_int_pos_sig = sum(int_p < a & int_eff > 0),
      n_int_neg_sig = sum(int_p < a & int_eff < 0),
      mean_eff = mean(treat_eff), med_eff = median(treat_eff),
      q05_eff = quantile(treat_eff, 0.05), q95_eff = quantile(treat_eff, 0.95))
  list(summary_df = summary_df, counts = counts)
}

best_index <- function(summary_df, cfg = CFG) {
  ord <- switch(cfg$best_rule,
    max_smd  = order(summary_df$max_smd,  summary_df$mean_smd),
    mean_smd = order(summary_df$mean_smd, summary_df$max_smd),
    stop("unknown best_rule"))
  summary_df$idx[ord[1]]
}

# Best-balanced member: refit with objects retained + CR2 country-clustered
# inference (FIX (v)). clubSandwich does not support lmer fits with prior
# weights, so the CR2 tier follows the standard post-matching route (MatchIt
# documentation): the SAME weighted fixed-effect regression estimated by OLS,
# with CR2 small-sample-corrected SEs clustered on country -- the coarsest
# dependence structure, mirroring the wells pipeline's country clustering and
# subsuming matched-set (subclass) dependence within countries.
refit_best <- function(spec, trdat, cndat, sets, best_idx, cfg = CFG) {
  b <- fit_one(trdat, cndat, sets[[best_idx]], spec, cfg, keep_fit = TRUE)
  stopifnot(isTRUE(b$success))
  cov_reg <- resolve_cov(spec$cov_reg, cfg)
  form_lm <- reformulate(c("type", cov_reg), response = spec$outcome)
  ols <- lm(form_lm, data = b$mdat, weights = weights_combined)
  cr2 <- tryCatch(as.data.frame(
    clubSandwich::coef_test(ols, vcov = "CR2", cluster = b$mdat$CntrName,
                            test = "Satterthwaite")),
    error = function(e) NULL)
  if (!is.null(cr2)) {
    rn <- if ("Coef" %in% names(cr2)) cr2$Coef else rownames(cr2)
    est <- if ("beta" %in% names(cr2)) cr2$beta else cr2$Estimate  # clubSandwich
    dfc <- if ("df" %in% names(cr2)) cr2$df else cr2$df_Satt       # version compat
    tr <- grepl("^type", rn); it <- rn == "(Intercept)"
    b$treat_eff_ols <- est[tr]
    b$treat_se_cr2  <- cr2$SE[tr]; b$treat_p_cr2 <- cr2$p_Satt[tr]
    b$int_eff_ols   <- est[it]
    b$int_se_cr2    <- cr2$SE[it]; b$int_p_cr2 <- cr2$p_Satt[it]
    b$cr2_df        <- dfc[tr]
    b$n_country     <- dplyr::n_distinct(b$mdat$CntrName)
  } else {
    b$treat_eff_ols <- NA_real_; b$treat_se_cr2 <- NA_real_
    b$treat_p_cr2 <- NA_real_;   b$int_eff_ols <- NA_real_
    b$int_se_cr2 <- NA_real_;    b$int_p_cr2 <- NA_real_
    b$cr2_df <- NA_real_;        b$n_country <- NA_integer_
  }
  b$idx <- best_idx
  b
}

run_spec <- function(spec_name, cfg = CFG, d = NULL, sets = NULL,
                     force = FALSE) {
  spec <- cfg$specs[[spec_name]]
  stopifnot(!is.null(spec))
  stamp <- cache_stamp(cfg, spec = c(spec_name, spec))
  path  <- file.path(cfg$cache_dir, paste0("spec_", spec_name, ".rds"))
  if (!force) {
    hit <- cache_read(path, stamp, cfg)
    if (!is.null(hit)) { say("  [cache] ", spec_name, "\n"); return(hit) }
  }
  if (is.null(d))    d    <- load_main_data(cfg)
  if (is.null(sets)) sets <- load_control_sets(cfg)

  pd <- prepare_spec_data(spec, d, cfg)
  sayf("  %s: outcome=%s  n_treat=%d  n_ctrl_pool=%d  realizations=%d\n",
       spec_name, spec$outcome, nrow(pd$trdat), nrow(pd$cndat), length(sets))
  t0 <- Sys.time()

  oplan <- future::plan(future::multisession, workers = par_workers(cfg))
  on.exit(future::plan(oplan), add = TRUE)
  # furrr_options(seed = TRUE) derives its per-element L'Ecuyer streams from
  # the CURRENT RNG state. Nothing in fit_one() draws today -- full optimal
  # matching and lmer are deterministic, and the control realizations are
  # pre-generated -- so the ensemble is reproducible either way. Fixing the
  # state anyway costs nothing and means that if a stochastic step is ever
  # added inside fit_one(), reproducibility does not break silently.
  results <- with_seed(cfg$ensemble_seed, future_map(sets, function(ids)
    fit_one(pd$trdat, pd$cndat, ids, spec, cfg),
    .options = furrr_options(seed = TRUE), .progress = FALSE))

  ens  <- summarise_ensemble(results, cfg)
  bidx <- best_index(ens$summary_df, cfg)
  # legacy comparison: which member the old (mean-SMD) rule would have chosen
  bidx_legacy <- ens$summary_df$idx[order(ens$summary_df$mean_smd)][1]
  best <- refit_best(spec, pd$trdat, pd$cndat, sets, bidx, cfg)

  # thin per-realization balance detail (per-covariate post-match SMD)
  bal_long <- imap_dfr(results, function(x, i) {
    if (!isTRUE(x$success)) return(tibble())
    x$balance %>% mutate(idx = as.integer(i))
  })

  res <- list(
    label       = spec_name,
    spec        = spec,
    summary_df  = ens$summary_df,
    counts      = ens$counts,
    balance_long= bal_long,
    best_idx    = bidx,
    best_idx_legacy_rule = bidx_legacy,
    best        = best[setdiff(names(best), "mdat")],
    n_treat     = nrow(pd$trdat),
    runtime_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  sayf("    done in %.1f min  (mean eff %.4g, %d%% of members p<%.2f, best idx %d)\n",
       res$runtime_min, res$counts$mean_eff,
       round(100 * res$counts$share_sig), cfg$sig_level, bidx)
  cache_write(res, path, stamp)
  res
}

run_specs <- function(names, cfg = CFG, force = FALSE) {
  d    <- load_main_data(cfg)
  sets <- load_control_sets(cfg)
  out  <- list()
  for (nm in names) out[[nm]] <- run_spec(nm, cfg, d, sets, force = force)
  out
}

# =============================================================================
# LORENZ CURVES FOR THE FIG. 3C ILLUSTRATION
# =============================================================================
# Reproduces the Gini construction of 2_buildDataset.R exactly (cummax
# monotonicity, bidirectional interpolation, trapezoid area), so the plotted
# curves and printed Gini values are the analysis quantities.
# Sign convention (and the corrected Methods formula):
#   G = 1 - 2 * integral(L(p) dp)  =  -2 * integral( (L(p) - p) dp )
# Negative = concentration toward the border; positive = interior.
lorenz_curve <- function(l, i, p = seq(0, 100, 10)) {
  keep <- !(is.na(l) | is.na(i))
  p <- p[keep]; l <- l[keep]; i <- i[keep]
  if (length(l) < 3) return(NULL)
  l <- cummax(l); i <- cummax(i)
  IrrPct_for_Land <- approx(x = i, y = p, xout = l, rule = 2, ties = "ordered")$y
  LandPct_for_Irr <- approx(x = l, y = p, xout = i, rule = 2, ties = "ordered")$y
  df <- rbind(data.frame(LandPct = p, IrrPct = IrrPct_for_Land),
              data.frame(LandPct = LandPct_for_Irr, IrrPct = p))
  df <- df[order(df$LandPct, df$IrrPct), ]
  df$LandPct <- pmax(pmin(df$LandPct, 100), 0)
  df$IrrPct  <- pmax(pmin(df$IrrPct, 100), 0)
  df <- df[!duplicated(round(df[, c("LandPct", "IrrPct")], 6)), ]
  x <- df$LandPct / 100; y <- df$IrrPct / 100
  area <- sum(diff(x) * (head(y, -1) + tail(y, -1))) / 2
  list(curve = df, gini = 1 - 2 * area)
}

lorenz_examples <- function(cfg = CFG) {
  meta <- load_meta(cfg)
  lp <- read.csv(cfg$landpct_file) %>% as_tibble()
  ip <- read.csv(cfg$irpct_file)   %>% as_tibble()
  pcols <- grep("^p\\d+$", names(lp), value = TRUE)
  pvals <- as.numeric(sub("p", "", pcols))
  out <- list()
  for (ex in cfg$lorenz_examples) {
    hit <- meta %>%
      filter(aq_id == ex$aq_id, CountryName == ex$country) %>% slice(1)
    if (!nrow(hit)) { say("  [lorenz] no metadata for aq_id ", ex$aq_id, "\n"); next }
    lrow <- lp %>% filter(aq_id == ex$aq_id, CntryNm == ex$country)
    irow <- ip %>% filter(aq_id == ex$aq_id, CntryNm == ex$country)
    if (!nrow(lrow) || !nrow(irow)) { say("  [lorenz] no curve for aq_id ", ex$aq_id, "\n"); next }
    lc <- lorenz_curve(as.numeric(lrow[1, pcols])[order(pvals)],
                       as.numeric(irow[1, pcols])[order(pvals)],
                       p = sort(pvals))
    if (is.null(lc)) next
    key <- paste(hit$AquiferName, ex$country, sep = " | ")
    out[[key]] <- c(lc, list(name = hit$AquiferName, country = ex$country))
  }
  out
}

# =============================================================================
# SI SUPPLEMENTARY DATASET S3 -- SEGMENT LEVEL (FIX (vii) + (viii))
# =============================================================================
# Rebuilds the treatment-only supplement table from the CURRENT master table,
# carrying over IGRAC metadata (names, codes, region) from the archived
# Data_S2 by aq_id x CountryName. Values (Irrig, GWIrrig, ginis, ...) come
# exclusively from the current build.
#
# FIX (viii): the emitted columns are exactly cfg$si_dataset_cols -- the
# variables the irrigation analysis actually uses. See the annotated manifest
# in 41_config.R for the per-column justification and for the six legacy
# columns that were dropped. Nothing here selects columns by hand: the
# transmute() builds the full renamed frame and the manifest subsets it, so
# the config is the single point of truth and drift is impossible.
build_si_irrigation_dataset <- function(cfg = CFG) {
  d    <- load_main_data(cfg)
  meta <- load_meta(cfg)
  full <- d %>% filter(type == "treat") %>%
    left_join(meta, by = c("aq_id", "CntrName" = "CountryName")) %>%
    transmute(
      aq_id, CountryName = CntrName, type, area_km2, lat_c, lon_c,
      CSI = CS_max, RDS_mainDist,
      GWIrrig = GW, Irrig = Ir, IrrigNeed = IrNeed3, IrNeed0,
      Overdraft = OverIR3, OverIR6,
      G_Irrig = giniIr, G_IrrigGW = giniGw, G_IrrigNeed = giniIr_Need,
      G_CSI = giniCSI, gini_crpGWSpct_gt0, gini_gwpct_gt0,
      NearRiverBorder = near_river_border,
      IGRAC_Code, CountryCode, Countries, AquiferName, Region)

  missing <- setdiff(cfg$si_dataset_cols, names(full))
  if (length(missing))
    stop("si_dataset_cols names columns the build does not produce: ",
         paste(missing, collapse = ", "))
  dropped <- setdiff(names(full), cfg$si_dataset_cols)
  if (length(dropped))
    say("  [Dataset S3] dropping ", length(dropped), " non-analysis column(s): ",
        paste(dropped, collapse = ", "), "\n")

  out <- full[, cfg$si_dataset_cols, drop = FALSE]

  # The dataset is the published record of the analysis sample: one row per
  # treated segment, no duplicate keys, and the two sample filters must still
  # reproduce the frozen counts.
  stopifnot(!any(duplicated(out[, c("aq_id", "CountryName")])))
  stopifnot(sum(out$Irrig   > 0, na.rm = TRUE) == cfg$reference_n_treat_ir,
            sum(out$GWIrrig > 0, na.rm = TRUE) == cfg$reference_n_treat_gw)
  out
}

# Back-compatible alias: the legacy name still works, so older drivers and
# any user script that calls regenerate_data_s2() keep running.
regenerate_data_s2 <- function(cfg = CFG) build_si_irrigation_dataset(cfg)
