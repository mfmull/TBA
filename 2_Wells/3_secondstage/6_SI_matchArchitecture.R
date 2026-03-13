

library(tidyverse)
library(MatchIt)
library(metafor)
# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script runs a sensitivity analysis for the second-stage pipeline by (i) implementing 
#nearest neighbor instead of full optimal matching and (ii) varying
# the matching ratio used in nearest-neighbor matching. For each ratio r (number
# of controls matched per treated unit), it re-estimates the pooled TB vs non-TB
# effects using the same matched-sample construction and the same multilevel
# meta-regression specification, then exports the results across ratios.
#
# Inputs
# - ../2_firststage/firststageMain.csv
#     Unit-level first-stage results for Aquifer × Country units, including:
#       * n, beta_0/se_0 (baseline/intercept-only estimate and SE)
#       * beta_1/se_1 (dist_LB_km slope estimate and SE)
#       * TB indicator, CC cluster id, and matching covariates
#
# Parameters
# - nMin: minimum wells per unit retained
# - cov_match: covariates used to balance TB vs non-TB units via matching
# - match: MatchIt method (here: nearest-neighbor matching)
# - r_seq: matching ratios evaluated (1 through 10)
#
# Helper logic
# - winsor(): winsorizes sampling variances to reduce the influence of extreme values
# - getstars(): converts p-values to significance stars for reporting
#
# Core routine (run_r)
# For a given matching ratio r, run_r():
# 1) Performs nearest-neighbor matching:
#    - Fits TB ~ cov_match using MatchIt with method = "nearest" and ratio = r
#    - Extracts the matched dataset and matching weights
# 2) Constructs meta-analytic sampling variances adjusted by matching weights:
#      vi_0 = (se_0^2) / weights,  vi_1 = (se_1^2) / weights
#    and winsorizes vi_0 and vi_1 (1%–99%).
# 3) Fits two multilevel meta-regressions (REML) with a country random intercept:
#    - Model 0: yi = beta_0, V = vi_0, mods = ~ 1 + TB, random = ~ 1 | CC
#    - Model 1: yi = beta_1, V = vi_1, mods = ~ 1 + TB, random = ~ 1 | CC
# 4) Returns a results tibble including the ratio, coefficients, SEs, p-values,
#    stars, tau^2 by CC, matched sample size, and number of CC clusters.
#
# Main workflow
# 1) Load and filter the unit-level dataset (n > nMin, valid coefficients/SEs, TB observed).
# 2) Loop over r_seq = 1:10:
#    - Run run_r(r) for each ratio and stack results into one table.
# 3) Export:
#    - Write the stacked results to matchArchOut.csv
#
# Output
# - matchArchOut.csv:
#     Pooled TB vs non-TB estimates (for beta_0 and beta_1 models) as a function
#     of the nearest-neighbor matching ratio.
# ------------------------------------------------------------------------------
##############################
#Params
##############################
nMin=20
cov_match =c("lat_c","lon_c","urbkHaKm2",'CS_max','LB_river')
match='nearest'

##############################
#Helpers
##############################
winsor <- function(x, p = c(0.01, 0.99)) {
  qs <- quantile(x, p, na.rm = TRUE)
  pmin(pmax(x, qs[1]), qs[2])
}
getstars=function(p){
  
  dplyr::case_when(
    p < 0.01 ~ "***",
    p < 0.05  ~ "**",
    p < 0.1  ~ "*",
    TRUE      ~ ""
  )
}
run_r <- function(r) {
  
  # Matching
  form_match <- reformulate(cov_match, response = "TB")
  
  m.out <- MatchIt::matchit(
    form_match,
    data   = aqf,
    method = match,
    ratio  = r
  )
  
  aqf_full <- MatchIt::match.data(m.out) %>%
    dplyr::mutate(
      vi_0 = winsor((se_0)^2 / weights),
      vi_1 = winsor((se_1)^2 / weights)
    )
  
  # Meta regressions
  fit_mv_0 <- metafor::rma.mv(
    yi = beta_0, V = vi_0, mods = ~ 1 + TB,
    random = ~ 1 | CC,
    data = aqf_full,
    method = "REML",
    control = list(optimizer = "optim", rel.tol = 1e-8)
  )
  
  fit_mv_1 <- metafor::rma.mv(
    yi = beta_1, V = vi_1, mods = ~ 1 + TB,
    random = ~ 1 | CC,
    data = aqf_full,
    method = "REML",
    control = list(optimizer = "optim", rel.tol = 1e-8)
  )
  
  # Output table
  out <- dplyr::tibble(
    r = r,
    
    b_0  = as.numeric(coef(fit_mv_0)),
    se_0 = sqrt(diag(vcov(fit_mv_0))),
    p_0  = 2 * pnorm(abs(b_0 / se_0), lower.tail = FALSE),
    nm_0 = names(coef(fit_mv_0)),
    st0  = getstars(p_0),
    tau2_CC_0 = fit_mv_0$sigma2[1],
    
    b_1  = as.numeric(coef(fit_mv_1)),
    se_1 = sqrt(diag(vcov(fit_mv_1))),
    p_1  = 2 * pnorm(abs(b_1 / se_1), lower.tail = FALSE),
    nm_1 = names(coef(fit_mv_1)),
    st1  = getstars(p_1),
    tau2_CC_1 = fit_mv_1$sigma2[1],
    
    n         = nrow(aqf_full),
    k_cluster = dplyr::n_distinct(aqf_full$CC)
  )
  
  return(out)
}

##############################
#Fetch and filter
##############################
aqf=read.csv('../2_firststage/firststageMain.csv')%>%
  filter(n>nMin)%>%
  filter(!is.na(beta_0), !is.na(se_0), se_0 > 0,
         !is.na(beta_1), !is.na(se_1), se_1 > 0)%>%
  filter(!is.na(TB)) 

##############################
#Loop through match ratio, using match method nearest
##############################
r_seq <- 1:10
results_ratio <- map_dfr(r_seq, function(r) {
  cat("Running ratio =", r, "\n")
  res <- run_r(
    r = r
  )
  res$ratio <- r
  res
})

future::plan(sequential)
write.csv(results_ratio,'matchArchOut.csv')
