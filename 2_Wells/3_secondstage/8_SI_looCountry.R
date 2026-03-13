# ==============================================================================
# Script purpose (preamble)
# ==============================================================================
# This script implements an overlap-country leave-one-out (LOO) robustness check
# for the preferred second-stage analysis comparing transboundary (TB) and
# non-transboundary aquifer–country segments.
#
# Context
# - Input data are unit-level (aquifer×country) first-stage estimates from
#   ../2_firststage/firststageMain.csv, including:
#     * beta_0, se_0 : mean depletion (intercept) estimate and SE
#     * beta_1, se_1 : depletion–distance-to-border slope estimate and SE
#     * TB           : treatment indicator (transboundary = TRUE)
#     * CC           : country identifier (cluster)
#     * n            : number of wells contributing to the unit estimate
#     * matching covariates (e.g., lat_c, lon_c, urbkHaKm2, CS_max, LB_river)
#
# What the script does
# 1) Filters the analytic sample (minimum wells per segment, valid betas/SEs).
# 2) Defines the “overlap country” set: countries that contain BOTH TB and
#    non-TB segments (i.e., countries that contribute within-country variation
#    for identifying the TB effect).
# 3) Fits the baseline specification on the full sample:
#    - Full matching (MatchIt; method = "full") of TB vs non-TB units using the
#      specified covariates.
#    - Two multilevel meta-regressions (metafor::rma.mv) using matched weights
#      and aquifer precision:
#        * Model 0: yi = beta_0 with V = se_0^2 / matched_weight
#        * Model 1: yi = beta_1 with V = se_1^2 / matched_weight
#      Both models include a random intercept by country (CC).
# 4) Performs LOO over overlap countries:
#    - For each overlap country, drops all units in that country,
#      reruns matching, refits both meta-regressions, and extracts the TBTRUE
#      coefficient, SE, and p-value for both models.
# 5) Adds metadata for plotting/interpretation:
#    - wells_in_CC: total wells in the omitted country (bubble size)
#    - delta_b_*  : change in TB coefficient relative to baseline
#    - sig_*      : significance indicator at alpha (default 0.10)
#    - (optional) standardized influence measures (infl_std_*)
# 6) Writes the LOO results table to loo_country_overlap_Out.csv.
#
# Outputs
# - loo_country_overlap_Out.csv: one row per omitted overlap country, containing
#   baseline and LOO TB coefficients for both models plus diagnostics.
#
# Notes
# - Variances are winsorized (1%–99%) after weight adjustment to reduce leverage
#   from extreme precision/weight combinations.
# - This LOO exercise is designed to assess whether results are sensitive to
#   a small number of countries contributing a disproportionate share of the
#   within-country identifying variation.
# ==============================================================================

library(tidyverse)
library(MatchIt)
library(metafor)
library(progressr)
library(ggplot2)
library(patchwork)
# ==============================
# Params
# ==============================
cov_match <- c("lat_c","lon_c","urbkHaKm2","CS_max","LB_river")
match <- "full"
r <- 1
nMin <- 20
alpha <- 0.10  # significance threshold shown in plots

# ==============================
# Helpers
# ==============================
winsor <- function(x, p = c(0.01, 0.99)) {
  qs <- quantile(x, p, na.rm = TRUE)
  pmin(pmax(x, qs[1]), qs[2])
}

extract_TBTRUE <- function(fit) {
  b <- coef(fit)
  se <- sqrt(diag(vcov(fit)))
  nm <- names(b)
  i <- which(nm == "TBTRUE")
  if (!length(i)) return(list(est = NA_real_, se = NA_real_))
  list(est = unname(b[i]), se = unname(se[i]))
}

# ==============================
# Load + filter
# ==============================
aqf <- read.csv('../2_firststage/firststageMain.csv') %>%
  filter(n > nMin) %>%
  filter(!is.na(beta_0), !is.na(se_0), se_0 > 0,
         !is.na(beta_1), !is.na(se_1), se_1 > 0) %>%
  filter(!is.na(TB), !is.na(CC))

# overlap countries only (both TB and non-TB)
tab_cc <- aqf %>% count(CC, TB) %>% pivot_wider(names_from = TB, values_from = n, values_fill = 0)
cc_frame <- tab_cc %>% filter(`TRUE` > 0, `FALSE` > 0) %>% pull(CC)

# wells per country (used for point size)
cc_wells <- aqf %>%
  group_by(CC) %>%
  summarise(wells_in_CC = sum(n, na.rm = TRUE), .groups = "drop")

# ==============================
# Core pipeline: matching + two meta models
# ==============================
fit_pipeline <- function(dat) {
  
  form_match <- reformulate(cov_match, response = "TB")
  m.out <- MatchIt::matchit(form_match, data = dat, method = match, ratio = r)
  
  dat_w <- MatchIt::match.data(m.out) %>%
    mutate(
      vi_0 = winsor((se_0)^2 / weights),
      vi_1 = winsor((se_1)^2 / weights)
    )
  
  fit_mv_0 <- metafor::rma.mv(
    yi = beta_0, V = vi_0, mods = ~ 1 + TB,
    random = ~ 1 | CC,
    data = dat_w,
    method = "REML",
    control = list(optimizer = "optim", rel.tol = 1e-8)
  )
  
  fit_mv_1 <- metafor::rma.mv(
    yi = beta_1, V = vi_1, mods = ~ 1 + TB,
    random = ~ 1 | CC,
    data = dat_w,
    method = "REML",
    control = list(optimizer = "optim", rel.tol = 1e-8)
  )
  
  tb0 <- extract_TBTRUE(fit_mv_0)
  tb1 <- extract_TBTRUE(fit_mv_1)
  
  tibble(
    tb_b_0  = tb0$est,
    tb_se_0 = tb0$se,
    tb_p_0  = 2 * pnorm(abs(tb_b_0 / tb_se_0), lower.tail = FALSE),
    
    tb_b_1  = tb1$est,
    tb_se_1 = tb1$se,
    tb_p_1  = 2 * pnorm(abs(tb_b_1 / tb_se_1), lower.tail = FALSE)
  )
}

# ==============================
# Baseline (full sample) WITH SEs
# ==============================
baseline <- fit_pipeline(aqf) %>%
  transmute(
    base_b_0  = tb_b_0,
    base_se_0 = tb_se_0,
    base_p_0  = tb_p_0,
    base_b_1  = tb_b_1,
    base_se_1 = tb_se_1,
    base_p_1  = tb_p_1
  )

base_b0  <- baseline$base_b_0[1]
base_se0 <- baseline$base_se_0[1]
base_b1  <- baseline$base_b_1[1]
base_se1 <- baseline$base_se_1[1]

# ==============================
# LOO runner (drop one overlap country)
# ==============================
run_loo_cc <- function(cc_drop) {
  dat <- aqf %>% filter(CC != cc_drop)
  out <- fit_pipeline(dat)
  out$drop_CC <- cc_drop
  out
}

# ==============================
# Run LOO with progress
# ==============================
handlers(global = TRUE)
handlers("progress")

resultsLOO <- with_progress({
  p <- progressor(steps = length(cc_frame))
  res <- vector("list", length(cc_frame))
  
  for (i in seq_along(cc_frame)) {
    cc <- cc_frame[i]
    p(message = paste("Dropping CC =", cc))
    res[[i]] <- tryCatch(run_loo_cc(cc), error = function(e) NULL)
  }
  bind_rows(res)
})

# ==============================
# Add metadata: wells, baseline, deltas, sig flags, standardized influence
# ==============================
resultsLOO <- resultsLOO %>%
  left_join(cc_wells, by = c("drop_CC" = "CC")) %>%
  mutate(
    sig_0 = tb_p_0 < alpha,
    sig_1 = tb_p_1 < alpha,
    
    base_b_0 = base_b0,
    base_b_1 = base_b1,
    
    delta_b_0 = tb_b_0 - base_b0,
    delta_b_1 = tb_b_1 - base_b1,
    
    infl_std_0 = (tb_b_0 - base_b0) / base_se0,
    infl_std_1 = (tb_b_1 - base_b1) / base_se1
  )

write.csv(resultsLOO, "loo_country_overlap_Out.csv", row.names = FALSE)

