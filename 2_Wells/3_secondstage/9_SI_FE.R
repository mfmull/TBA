library(tidyverse)
library(MatchIt)
library(metafor)
# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script runs a second-stage meta-regression that pools Aquifer × Country
# first-stage estimates and estimates how outcomes differ for transboundary
# (TB = 1) versus non-transboundary (TB = 0) units, after covariate balancing via
# matching. Unlike the multilevel random-effects version, this specification
# treats country (CC) as a fixed effect by including factor(CC) directly in the
# meta-regression.
#
# Inputs
# - ../2_firststage/firststageMain.csv
#     Unit-level first-stage results (Aquifer × CC), including:
#       * n, beta_0/se_0 (intercept-only estimate and SE)
#       * beta_1/se_1 (dist_LB_km slope estimate and SE)
#       * TB indicator, CC identifier, and matching covariates
#
# Parameters
# - nMin: minimum wells per unit retained
# - cov_match: covariates used for balancing TB vs non-TB units
# - match: MatchIt method (here: full matching)
# - r: matching ratio (used by some methods)
#
# Workflow
# 1) Load and filter unit-level data:
#    - Keep units with n > nMin
#    - Drop units with missing coefficients/SEs or non-positive SEs
#    - Require TB to be observed
# 2) Matching / balancing:
#    - Estimate TB ~ cov_match with MatchIt and generate matched-sample weights
#    - Extract matched dataset and combine weights into sampling variances:
#        vi_0 = (se_0^2) / weights,  vi_1 = (se_1^2) / weights
#      then winsorize variances (1%–99%) to limit extreme influence.
# 3) Fixed-effects meta-regression (country FE):
#    - Fit two metafor::rma.uni models (REML), each including:
#        * intercept, TB, and country fixed effects factor(CC)
#      Model 0 uses yi = beta_0 and vi = vi_0;
#      Model 1 uses yi = beta_1 and vi = vi_1.
# 4) Export results:
#    - Extract coefficients, SEs, normal-approx p-values, and significance stars
#    - Report matched sample size (n) and number of CC clusters represented
#
# Output
# - FE_Out.csv:
#     Coefficient table for both fixed-effects meta-regressions (baseline and
#     border-distance slope), including TB and CC fixed effects.
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
  metafor::rma.uni(
    yi = beta_0, vi = vi_0, mods = ~ 1 + TB + factor(CC),
    data = aqf_full, method = "REML",
    control = list(optimizer = "optim", rel.tol = 1e-8)
  )

fit_mv_1 <-
  metafor::rma.uni(
    yi = beta_1, vi = vi_1, mods = ~ 1 + TB + factor(CC),
    data = aqf_full, method = "REML",
    control = list(optimizer = "optim", rel.tol = 1e-8)
  )


out=dplyr::tibble(
  b_0  = as.numeric(coef(fit_mv_0)),
  se_0 = sqrt(diag(vcov(fit_mv_0))),
  p_0  = 2 * pnorm(abs(b_0 / se_0), lower.tail = FALSE),
  nm_0 = names(coef(fit_mv_0)),
  st0= getstars(p_0),
  
  b_1  = as.numeric(coef(fit_mv_1)),
  se_1 = sqrt(diag(vcov(fit_mv_1))),
  p_1  = 2 * pnorm(abs(b_1 / se_1), lower.tail = FALSE),
  nm_1 = names(coef(fit_mv_1)),
  st1= getstars(p_1),
  
  n         = nrow(aqf_full),
  k_cluster = dplyr::n_distinct(aqf_full$CC))

# saveRDS(m.out, file = "match.rds")
write.csv(out,'FE_Out.csv')