
  
library(tidyverse)
library(MatchIt)
library(metafor)
# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script reruns the preferred second-stage “matching + multilevel meta-regression”
# analysis after excluding a predefined set of high-leverage / top-heavy aquifers.
# The goal is to check robustness of the pooled TB vs non-TB effects to removing
# a small number of aquifers that may dominate results because they contribute
# disproportionately large samples or influence.
#
# Inputs
# - ../2_firststage/firststageMain.csv
#     Unit-level first-stage results for Aquifer × Country units, including:
#       * n, beta_0/se_0 (baseline/intercept-only estimate and SE)
#       * beta_1/se_1 (dist_LB_km slope estimate and SE)
#       * TB indicator, CC cluster id, and matching covariates
#
# Parameters
# - nMin: minimum wells per unit retained (filters out small units)
# - cov_match: covariates used to balance TB vs non-TB units via matching
# - match: MatchIt method (here: full matching)
# - r: matching ratio (used by some methods)
# - drop list: a hard-coded set of aquifer names removed before matching and pooling
#
# Workflow
# 1) Load and filter unit-level dataset:
#    - Keep units with n > nMin
#    - Drop units with missing coefficients/SEs or non-positive SEs
#    - Require TB to be observed
# 2) Remove “top-heavy” aquifers:
#    - Exclude units whose Aquifer name is in the provided list (sensitivity check
#      targeting potentially influential aquifers).
# 3) Covariate matching / balancing:
#    - Fit TB ~ cov_match using MatchIt (full matching)
#    - Extract matched sample and weights
#    - Construct sampling variances adjusted by matching weights:
#        vi_0 = (se_0^2) / weights,  vi_1 = (se_1^2) / weights
#      then winsorize variances (1%–99%) to limit the influence of extremes.
# 4) Second-stage multilevel meta-regressions (REML):
#    - Model 0: yi = beta_0, V = vi_0, mods = ~ 1 + TB, random intercept by CC
#    - Model 1: yi = beta_1, V = vi_1, mods = ~ 1 + TB, random intercept by CC
# 5) Export results:
#    - Extract coefficients, SEs, p-values, stars, and between-country heterogeneity
#      (tau^2 for CC), plus matched sample size and number of CC clusters.
#
# Output
# - dropTopSEOur.csv:
#     Summary table of pooled TB vs non-TB estimates (for beta_0 and beta_1 models)
#     after excluding the specified “top-heavy” aquifers.
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
  
#drop top-heavy
aqf=aqf%>%filter(!Aquifer%in%c("Biscayne Aquifer","Central Minnesota Surficial and Buried Sand and Gravel Aquifers","East Bay Plain","Santa Clara Valley", "Orange County Coastal Plain" ,"Northern High Plains","Dakota Aquifer System","Delmarva Peninsula","Biscayne Aquifer", "Central Mississippi Embayment"))

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
  
# saveRDS(m.out, file = "match.rds")
write.csv(out,'dropTopSEOur.csv')
  
  
