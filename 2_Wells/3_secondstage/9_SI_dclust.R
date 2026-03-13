

library(tidyverse)
library(MatchIt)
library(metafor)
# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script runs a sensitivity analysis for the second-stage meta-regression by
# repeating the full “matching + multilevel meta-regression” workflow separately
# for each spatial declustering setting (dclust). The input file contains the
# first-stage results already recomputed under different declustering radii (and
# optionally KDE-based declustering), and this script summarizes how the pooled
# TB vs non-TB effects change across those dclust values.
#
# Inputs
# - ../2_firststage/firststageDCKDE.csv
#     Unit-level first-stage outputs for Aquifer × Country units, repeated across
#     declustering settings `dclust`. Includes beta_0/se_0, beta_1/se_1, TB, CC,
#     matching covariates, and the dclust identifier.
#
# Parameters
# - nMin: minimum wells per unit retained
# - cov_match: covariates used for balancing TB vs non-TB units
# - match: MatchIt method (here: full matching)
# - r: matching ratio (used by some matching methods)
#
# Core routine (run_clst)
# For the currently active dataset `aqf` (set inside the dclust loop), run_clst():
# 1) Performs covariate matching:
#    - Fits TB ~ cov_match using MatchIt
#    - Extracts matched sample with MatchIt weights
#    - Builds sampling variances adjusted by matching weights:
#        vi_0 = (se_0^2) / weights,  vi_1 = (se_1^2) / weights
#      then winsorizes variances (1%–99%) to reduce the influence of extremes.
# 2) Fits two multilevel meta-regressions (REML):
#    - Model 0: yi = beta_0, V = vi_0, mods = ~ 1 + TB, random intercept by CC
#    - Model 1: yi = beta_1, V = vi_1, mods = ~ 1 + TB, random intercept by CC
# 3) Returns a compact results tibble with coefficients, SEs, p-values, stars,
#    and tau^2 (between-country heterogeneity), plus n and number of CC clusters.
#
# Main workflow
# 1) Load and filter the declustering first-stage dataset:
#    - Keep units with n > nMin
#    - Drop units with missing coefficients/SEs or non-positive SEs
#    - Require TB to be observed
# 2) Split by dclust and iterate:
#    - For each dclust value, set aqf to that subset and call run_clst()
#    - Append the current dclust value to the returned results
# 3) Export:
#    - Write the stacked results across all dclust settings to dclustOut.csv
#
# Output
# - dclustOut.csv:
#     Second-stage pooled TB vs non-TB estimates (for beta_0 and beta_1 models)
#     as a function of declustering radius setting `dclust`.
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
run_clst <- function() {
  
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

aqf=read.csv('../2_firststage/firststageDCKDE.csv')%>%
  filter(n>nMin)%>%
  filter(!is.na(beta_0), !is.na(se_0), se_0 > 0,
         !is.na(beta_1), !is.na(se_1), se_1 > 0)%>%
  filter(!is.na(TB)) 

res_all <- aqf %>% arrange(dclust)%>%
  dplyr::group_split(dclust) %>%
  purrr::map_dfr(function(dat) {
    message("Running dclust = ", dat$dclust[1], " (n = ", nrow(dat), ")")
    aqf <<- dat
    dplyr::mutate(run_clst(), dclust = dat$dclust[1])
  })


write.csv(res_all,'dclustOut.csv')


