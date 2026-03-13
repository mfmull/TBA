

library(tidyverse)
library(MatchIt)
library(metafor)
library(clubSandwich)
# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script runs the second-stage meta-regression that pools Aquifer × Country
# first-stage estimates and compares transboundary (TB = 1) versus non-transboundary
# (TB = 0) units, after covariate balancing via matching. Relative to the basic
# version, it reports country-cluster-robust inference (CR2) for the fixed effects.
#
# Inputs
# - ../2_firststage/firststageMain.csv
#     Unit-level first-stage results (Aquifer × CC), including:
#       * n, beta_0/se_0 (intercept-only estimate and SE)
#       * beta_1/se_1 (dist_LB_km slope estimate and SE)
#       * TB indicator, CC cluster id, and matching covariates
#
# Parameters
# - nMin: minimum wells per unit retained
# - cov_match: covariates used to balance TB and non-TB units (matching on TB)
# - match: MatchIt method (here: full matching)
# - r: matching ratio (used by some matching methods)
#
# Workflow
# 1) Load and filter the unit-level dataset:
#    - Keep units with n > nMin
#    - Drop units with missing coefficients/SEs or non-positive SEs
#    - Require TB to be observed
# 2) Covariate matching / balancing:
#    - Fit a matching model TB ~ cov_match using MatchIt
#    - Extract matched data with MatchIt weights
#    - Build meta-analytic sampling variances adjusted by matching weights:
#        vi_0 = (se_0^2) / weights,  vi_1 = (se_1^2) / weights
#      then winsorize variances (1%–99%) to limit extreme influence.
# 3) Second-stage multilevel meta-regressions (REML):
#    - Fit two models with metafor::rma.mv, each with:
#        * fixed effects: intercept + TB
#        * random intercept: by country (CC)
#      Model 0 uses yi = beta_0, Model 1 uses yi = beta_1.
# 4) Cluster-robust inference (clubSandwich):
#    - Use coef_test(..., vcov = "CR2", cluster = CC) to compute CR2 robust SEs,
#      small-sample degrees of freedom, and p-values for the fixed effects.
#    - A helper (pick_p) selects the p-value column name returned by coef_test,
#      and stopifnot checks ensure coefficient ordering matches.
# 5) Export results:
#    - Output table includes point estimates, CR2 SEs, p-values, stars, df,
#      between-country heterogeneity (tau^2 for CC), plus matched sample size
#      and number of country clusters.
#
# Output
# - clustRobustSE_Out.csv:
#     Summary of fixed-effect estimates for both meta-regressions with CR2
#     cluster-robust inference and tau^2 by country.
# ------------------------------------------------------------------------------
##############################
#Params
##############################
nMin=20
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
##############################
#Fetch and filter
##############################

aqf=read.csv('../2_firststage/firststageMain.csv')%>%
  filter(n>nMin)%>%
  filter(!is.na(beta_0), !is.na(se_0), se_0 > 0,
         !is.na(beta_1), !is.na(se_1), se_1 > 0)%>%
  filter(!is.na(TB)) 

##############################
#Matching on covariates
##############################
form_match <- reformulate(cov_match, response = "TB")
m.out <- matchit(form_match, data = aqf, method = match,ratio=r)
aqf_full <- match.data(m.out) %>%
  dplyr::mutate(
    # Combine weights
    vi_0 = winsor((se_0)^2 / weights),
    vi_1 = winsor((se_1)^2 / weights)
  )


##############################
#meta-regress (2nd stage)
############################## 
fit_mv_0 <- 
  metafor::rma.mv(
    yi = beta_0, V = vi_0, mods = ~ 1 + TB,
    random = ~ 1 | CC,
    data = aqf_full, method = "REML",
    control = list(optimizer = "optim", rel.tol = 1e-8)
  )

fit_mv_1 <- 
  metafor::rma.mv(
    yi = beta_1, V = vi_1, mods = ~ 1 + TB,
    random = ~ 1 | CC,
    data = aqf_full, method = "REML",
    control = list(optimizer = "optim", rel.tol = 1e-8)
  )

rob_0 <- coef_test(fit_mv_0, vcov = "CR2", cluster = aqf_full$CC)
rob_1 <- coef_test(fit_mv_1, vcov = "CR2", cluster = aqf_full$CC)

pick_p <- function(x){
  if ("p_Satt" %in% names(x)) return(x$p_Satt)
  if ("p_val"  %in% names(x)) return(x$p_val)
  if ("p"      %in% names(x)) return(x$p)
  stop("No p-value column found in coef_test output.")
}

# sanity check: same coefficient ordering
stopifnot(identical(names(coef(fit_mv_0)), rownames(rob_0)))
stopifnot(identical(names(coef(fit_mv_1)), rownames(rob_1)))

out <- dplyr::tibble(
  # point estimates (same either way)
  b_0  = as.numeric(coef(fit_mv_0)),
  se_0 = rob_0$SE,
  p_0  = pick_p(rob_0),
  nm_0 = names(coef(fit_mv_0)),
  st0  = getstars(p_0),
  df_0 = rob_0$df,
  tau2_CC_0 = fit_mv_0$sigma2[1],
  
  b_1  = as.numeric(coef(fit_mv_1)),
  se_1 = rob_1$SE,
  p_1  = pick_p(rob_1),
  nm_1 = names(coef(fit_mv_1)),
  st1  = getstars(p_1),
  df_1 = rob_1$df,
  tau2_CC_1 = fit_mv_1$sigma2[1],
  
  n         = nrow(aqf_full),
  k_cluster = dplyr::n_distinct(aqf_full$CC)
)

write.csv(out, "clustRobustSE_Out.csv", row.names = FALSE)

