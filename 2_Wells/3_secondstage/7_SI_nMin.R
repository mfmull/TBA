

library(tidyverse)
library(MatchIt)
library(metafor)
# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script runs a sensitivity analysis on the minimum unit sample-size
# threshold (nMin) used in the second-stage meta-regression. It repeatedly
# executes the same “matching + multilevel meta-regression” pipeline while
# varying the cutoff for the number of wells per Aquifer × Country unit, and
# exports how the pooled TB vs non-TB effects change as nMin increases.
#
# Inputs
# - ../2_firststage/firststageMain.csv
#     Unit-level first-stage results (Aquifer × CC), including:
#       * n (wells per unit)
#       * beta_0/se_0 (baseline/intercept-only estimate and SE)
#       * beta_1/se_1 (dist_LB_km slope estimate and SE)
#       * TB indicator, CC cluster id, and matching covariates
#
# Parameters
# - cov_match: covariates used to balance TB vs non-TB units
# - match: MatchIt method (here: full matching)
# - r: matching ratio (used by some methods)
# - nMin_seq: sequence of minimum-well thresholds to evaluate (10…100)
#
# Core routine (run_nMin)
# For a given threshold nMin, run_nMin():
# 1) Filters the full unit dataset to units with n > nMin.
# 2) Performs covariate matching on TB:
#    - Fits TB ~ cov_match using MatchIt
#    - Extracts matched sample and weights
#    - Constructs sampling variances adjusted by matching weights:
#        vi_0 = (se_0^2) / weights,  vi_1 = (se_1^2) / weights
#      then winsorizes variances (1%–99%) to limit extreme influence.
# 3) Fits two multilevel meta-regressions (REML) with a country random intercept:
#    - Model 0: yi = beta_0, V = vi_0, mods = ~ 1 + TB, random = ~ 1 | CC
#    - Model 1: yi = beta_1, V = vi_1, mods = ~ 1 + TB, random = ~ 1 | CC
# 4) Returns a compact results tibble with coefficients, SEs, p-values, stars,
#    and tau^2 (between-country heterogeneity), plus matched sample size and
#    number of CC clusters represented.
#
# Main workflow
# 1) Read first-stage results and apply base validity filters (non-missing,
#    positive SEs, TB observed).
# 2) Loop over nMin_seq (10, 20, …, 100), run run_nMin() for each threshold, and
#    append the current nMin to the returned results.
# 3) Export:
#    - Write the stacked sensitivity results to nMinOut.csv
#
# Output
# - nMinOut.csv:
#     Second-stage pooled TB vs non-TB estimates (for beta_0 and beta_1 models)
#     as a function of the minimum wells-per-unit threshold nMin.
# ------------------------------------------------------------------------------
##############################
#Params
##############################
cov_match =c("lat_c","lon_c","urbkHaKm2",'CS_max','LB_river')
match='full'
r=1
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
run_nMin <- function(nMin) {
  aqfi=aqf%>%filter(n>nMin)
  # Matching
  form_match <- reformulate(cov_match, response = "TB")
  
  m.out <- MatchIt::matchit(
    form_match,
    data   = aqfi,
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
  filter(!is.na(beta_0), !is.na(se_0), se_0 > 0,
         !is.na(beta_1), !is.na(se_1), se_1 > 0)%>%
  filter(!is.na(TB)) 

##############################
#Loop through match ratio, using match method nearest
##############################
nMin_seq <- c(10,20,30,40,50,60,70,80,90,100)
results_nMin <- map_dfr(nMin_seq, function(nm) {
  cat("Running nMin =", nm, "\n")
  res <- run_nMin(
    nMin = nm
  )
  res$nMin <- nm
  res
})

write.csv(results_nMin,'nMinOut.csv')
