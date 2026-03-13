
  
library(tidyverse)
library(MatchIt)
library(metafor)
# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script runs the preferred second-stage “matching + multilevel meta-regression”
# analysis using first-stage estimates produced by robust regression (MASS::rlm),
# rather than OLS with HC standard errors. It evaluates whether the pooled TB vs
# non-TB differences are sensitive to outliers and high-leverage wells at the
# first-stage estimation step.
#
# Inputs
# - ../2_firststage/firstStageMainRobust.csv
#     Unit-level first-stage results for Aquifer × Country units estimated via
#     robust regression, including:
#       * n, beta_0/se_0 (baseline/intercept-only estimate and SE)
#       * beta_1/se_1 (dist_LB_km slope estimate and SE)
#       * TB indicator, CC cluster id, and matching covariates
#
# Parameters
# - nMin: minimum wells per unit retained
# - cov_match: covariates used to balance TB vs non-TB units via matching
# - match: MatchIt method (here: full matching)
# - r: matching ratio (used by some methods)
#
# Workflow
# 1) Load and filter the robust first-stage dataset:
#    - Keep units with n > nMin
#    - Drop units with missing coefficients/SEs or non-positive SEs
#    - Require TB to be observed
# 2) Covariate matching / balancing:
#    - Fit TB ~ cov_match using MatchIt (full matching)
#    - Extract matched sample and weights
#    - Construct meta-analytic sampling variances adjusted by matching weights:
#        vi_0 = (se_0^2) / weights,  vi_1 = (se_1^2) / weights
#      then winsorize variances (1%–99%) to limit extreme influence.
# 3) Second-stage multilevel meta-regressions (REML):
#    - Model 0: yi = beta_0, V = vi_0, mods = ~ 1 + TB, random intercept by CC
#    - Model 1: yi = beta_1, V = vi_1, mods = ~ 1 + TB, random intercept by CC
# 4) Export results:
#    - Extract coefficients, SEs, p-values, stars, and between-country heterogeneity
#      (tau^2 for CC), plus matched sample size and number of CC clusters.
#
# Output
# - robustOut.csv:
#     Summary table of pooled TB vs non-TB estimates (for beta_0 and beta_1 models)
#     using robust first-stage estimates with the standard matching + rma.mv
#     second-stage specification.
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
 
aqf=read.csv('../2_firststage/firstStageMainRobust.csv')%>%
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
  
write.csv(out,'robustOut.csv')
  
  
