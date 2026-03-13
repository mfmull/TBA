
  
library(tidyverse)
library(MatchIt)
library(metafor)
# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script performs the second-stage analysis that aggregates the unit-level
# first-stage estimates (Aquifer × Country units) into overall effects, comparing
# transboundary (TB = 1) versus non-transboundary (TB = 0) aquifers after
# covariate balancing via matching.
#
# Inputs
# - ../2_firststage/firststageMain.csv
#     Unit-level first-stage regression outputs, including:
#       * n (well count per unit)
#       * beta_0, se_0: intercept-only (baseline) estimate and SE
#       * beta_1, se_1: dist_LB_km slope estimate and SE
#       * TB indicator, CC (country code), and matching covariates (e.g., lat_c, lon_c, etc.)
#
# Parameterization
# - nMin: minimum wells per unit (drops small-sample units)
# - cov_match: covariates used to match TB and non-TB units
# - match: MatchIt matching method (here: full matching)
# - r: matching ratio (used by some methods; included for completeness)
#
# Main steps
# 1) Load and filter unit-level data:
#    - Keep units with n > nMin
#    - Drop units with missing or non-positive SEs for either first-stage effect
#    - Require TB to be observed
# 2) Covariate balancing via matching:
#    - Estimate a propensity model TB ~ covariates (form_match)
#    - Run MatchIt (method = "full") to create a matched sample with weights
#    - Extract matched dataset with weights (match.data)
# 3) Construct meta-analytic sampling variances:
#    - For each first-stage effect (0 and 1), compute vi = (se^2) / weights
#    - Winsorize vi to limit the influence of extreme variances (default 1%–99%)
# 4) Second-stage meta-regression (random effects by country):
#    - Fit two multilevel meta-regressions with metafor::rma.mv:
#        * Model 0: yi = beta_0, V = vi_0, mods = ~ 1 + TB, random intercept by CC
#        * Model 1: yi = beta_1, V = vi_1, mods = ~ 1 + TB, random intercept by CC
#      Both estimated by REML with tighter optimizer tolerances.
# 5) Compile results for export:
#    - Extract coefficients, SEs, p-values, and significance stars
#    - Report between-country heterogeneity (tau^2 for CC random effect)
#    - Report matched sample size (n) and number of country clusters (k_cluster)
#
# Outputs
# - match.rds:
#     Saved MatchIt object (matching specification + diagnostics)
# - preferredOut.csv:
#     One-row summary table with meta-regression coefficients for both models,
#     their uncertainty (SE, p, stars), heterogeneity (tau^2), and sample counts.
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
  
saveRDS(m.out, file = "match.rds")
write.csv(out,'preferredOut.csv')
  
# ------------------------------------------------------------------------------
# Share of matched controls in same country as treated unit
# ------------------------------------------------------------------------------

same_country_stats <- aqf_full %>%
  mutate(subclass = as.character(subclass)) %>%
  group_by(subclass) %>%
  mutate(
    treated_cc = ifelse(any(TB == 1), first(unique(CC[TB == 1])), NA_character_)
  ) %>%
  ungroup() %>%
  filter(TB == 0, !is.na(treated_cc)) %>%
  mutate(
    same_country = CC == treated_cc
  )

# Unweighted counts
n_ctrl_total <- nrow(same_country_stats)
n_ctrl_same  <- sum(same_country_stats$same_country, na.rm = TRUE)
share_ctrl_same <- n_ctrl_same / n_ctrl_total

# Weighted shares using MatchIt weights
w_ctrl_total <- sum(same_country_stats$weights, na.rm = TRUE)
w_ctrl_same  <- sum(same_country_stats$weights[same_country_stats$same_country], na.rm = TRUE)
share_ctrl_same_w <- w_ctrl_same / w_ctrl_total

# Print to console
cat("Matched controls in same country as treated unit:\n")
cat("Unweighted:", n_ctrl_same, "/", n_ctrl_total,
    sprintf("(%.1f%%)\n", 100 * share_ctrl_same))
cat("Weighted:", round(w_ctrl_same, 2), "/", round(w_ctrl_total, 2),
    sprintf("(%.1f%%)\n", 100 * share_ctrl_same_w))
