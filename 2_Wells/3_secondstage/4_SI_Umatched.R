
library(tidyverse)
library(metafor)
# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script estimates the pooled second-stage TB vs non-TB differences in the
# first-stage outcomes *without* any covariate matching. It directly meta-analyzes
# the Aquifer × Country unit estimates, weighting only by first-stage uncertainty,
# and includes a random intercept by country (CC) to account for within-country
# dependence/heterogeneity.
#
# Inputs
# - ../2_firststage/firststageMain.csv
#     Unit-level first-stage results for Aquifer × CC, including:
#       * n, beta_0/se_0 (baseline/intercept-only estimate and SE)
#       * beta_1/se_1 (dist_LB_km slope estimate and SE)
#       * TB indicator and CC identifier
#
# Parameters
# - nMin: minimum wells per unit retained (filters out small-sample units)
#
# Workflow
# 1) Load and filter unit-level data:
#    - Keep units with n > nMin
#    - Drop units with missing coefficients/SEs or non-positive SEs
#    - Require TB to be observed
# 2) Construct meta-analytic sampling variances (no matching):
#    - Use the first-stage variances only:
#        vi_0 = se_0^2,  vi_1 = se_1^2
#      with winsorization (1%–99%) to limit the influence of extreme values.
# 3) Second-stage multilevel meta-regressions (REML):
#    - Fit two metafor::rma.mv models with fixed effects (intercept + TB) and a
#      random intercept by CC:
#        * Model 0: yi = beta_0, V = vi_0, mods = ~ 1 + TB, random = ~ 1 | CC
#        * Model 1: yi = beta_1, V = vi_1, mods = ~ 1 + TB, random = ~ 1 | CC
# 4) Export results:
#    - Extract coefficients, SEs, p-values, stars, and tau^2 (between-country
#      heterogeneity), plus sample size and number of CC clusters.
#
# Output
# - unMachedOut.csv:
#     Summary table of pooled TB vs non-TB estimates (for beta_0 and beta_1 models)
#     using only first-stage uncertainty weights (no matching/balancing).
# ------------------------------------------------------------------------------
##############################
#Params
##############################
nMin=20


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
#No Matching on covariates
##############################
aqf_full <- aqf%>%
  dplyr::mutate(
    # first stage uncertainty weight
    vi_0 = winsor((se_0)^2 ),
    vi_1 = winsor((se_1)^2 )
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
  
  
out=dplyr::tibble(
    b_0  = as.numeric(coef(fit_mv_0)),
    se_0 = sqrt(diag(vcov(fit_mv_0))),
    p_0  = 2 * pnorm(abs(b_0 / se_0), lower.tail = FALSE),
    nm_0 = names(coef(fit_mv_0)),
    st0= getstars(p_0),
    tau2_CC_0  = fit_mv_0$sigma2[1],
    
    b_1  = as.numeric(coef(fit_mv_1)),
    se_1 = sqrt(diag(vcov(fit_mv_1))),
    p_1  = 2 * pnorm(abs(b_1 / se_1), lower.tail = FALSE),
    nm_1 = names(coef(fit_mv_1)),
    st1= getstars(p_1),
    tau2_CC_1  = fit_mv_1$sigma2[1],
    n         = nrow(aqf_full),
    k_cluster = dplyr::n_distinct(aqf_full$CC))
  
write.csv(out,'unMachedOut.csv')
  
  
