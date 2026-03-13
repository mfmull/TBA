

library(tidyverse)
library(MatchIt)
library(metafor)
library(progressr)
# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script runs a robustness / influence sensitivity analysis for the preferred
# second-stage pipeline by repeatedly re-estimating the pooled TB vs non-TB effects
# after randomly dropping a subset of transboundary aquifers (“leave-k-out”).
# For each draw, it reruns the full matching + multilevel meta-regression workflow
# and records how the estimated TB effect and heterogeneity change as more TB
# aquifers are removed.
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
# - cov_match: covariates used to balance TB vs non-TB units in matching
# - match: MatchIt method (here: full matching)
# - r: matching ratio (used by some methods)
# - M: number of repetitions (Monte Carlo draws)
# - k_values: numbers of TB aquifers to drop each draw (k = 1, 5, 10, …)
#   (the script currently includes a smaller test setting that overwrites M and k_values)
#
# Core routine (run_loo)
# Given a vector `loo` of aquifer names to drop:
# 1) Filter the dataset to exclude those aquifers.
# 2) Re-run covariate matching:
#    - Fit TB ~ cov_match using MatchIt
#    - Extract matched sample and weights
#    - Construct sampling variances adjusted by matching weights:
#        vi_0 = (se_0^2) / weights,  vi_1 = (se_1^2) / weights
#      then winsorize variances (1%–99%) to limit extreme influence.
# 3) Re-fit the two multilevel meta-regressions (REML), with random intercept by CC:
#    - Model 0: yi = beta_0, V = vi_0, mods = ~ 1 + TB
#    - Model 1: yi = beta_1, V = vi_1, mods = ~ 1 + TB
# 4) Return a compact results tibble with coefficients, SEs, p-values, stars,
#    tau^2 by CC, and matched sample size / number of CC clusters.
#
# Main workflow
# 1) Load and filter the unit-level dataset (n > nMin, valid coefficients/SEs, TB observed).
# 2) Define the sampling frame of aquifers to drop:
#    - TBA = list of aquifer names with TB == TRUE.
# 3) Monte Carlo leave-k-out loop with progress reporting:
#    - For each repetition m = 1…M and each k in k_values:
#        * sample k TB aquifers (seeded for reproducibility)
#        * run run_loo() on the reduced dataset (tryCatch to skip failures)
#        * store which aquifers were dropped (looaqf), k, and iteration id
#    - Combine results across repetitions into one long table.
# 4) Export:
#    - Write the stacked sensitivity results to loo_Out.csv
#    - Reset future plan to sequential.
#
# Output
# - loo_Out.csv:
#     A draw-by-draw record of pooled TB vs non-TB estimates (for beta_0 and beta_1),
#     including which TB aquifers were dropped, how many (k), and repetition index.
# ------------------------------------------------------------------------------
##############################
#Params
##############################
cov_match =c("lat_c","lon_c","urbkHaKm2",'CS_max','LB_river')
match='full'
r=1
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
run_loo <- function(loo) {
  
  aqfi=aqf%>%filter(!Aquifer%in%loo)
  
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
  filter(n>nMin)%>%
  filter(!is.na(beta_0), !is.na(se_0), se_0 > 0,
         !is.na(beta_1), !is.na(se_1), se_1 > 0)%>%
  filter(!is.na(TB)) 

##############################
#Loop through match ratio, using match method nearest
##############################
# # Vector of aquifers to sample from
TBA <- aqf %>% filter(TB) %>% pull(Aquifer)
# 

M <- 250                               # Number of repetitions
k_values <- c(1, 5, 10, 15, 20, 25, 30, 35, 40)



M <- 2                               # Number of repetitions
k_values <- c(1, 5)

# 
handlers(global = TRUE)
handlers("progress")

resultsLOO <- with_progress({
  total_steps <- M * length(k_values)
  p <- progressor(steps = total_steps)

  results_all <- vector("list", M)

  for (m in seq_len(M)) {
    results_m <- vector("list", length(k_values))

    for (ki in seq_along(k_values)) {
      k <- k_values[ki]
      p(message = sprintf("Iteration m=%d, dropping k=%d aquifers", m, k))

      set.seed(1000 + k * 100 + m)
      TBA_sel <- sample(TBA, k)

      res <- tryCatch(run_loo(loo = TBA_sel), error = function(e) NULL)
      if (!is.null(res)) {
        res$looaqf <- paste(TBA_sel, collapse = "; ")
        res$k <- k
        res$iter <- m
      }

      results_m[[ki]] <- res
    }

    # Combine all results from this repetition
    results_all[[m]] <- bind_rows(results_m)

    # Temporary save after each repetition
    # save(results_all, file = 'resultsLOO_final.rdata')
  }

  bind_rows(results_all)
})



future::plan(sequential)
write.csv(resultsLOO,'loo_Out.csv')
