# =============================================================================
# Script: pilot_match_batch.R
#
# Purpose
#   Run a batch of matching + mixed-effects outcome models across many candidate
#   control sets (nonOverlapsB). For each run:
#     1) Keep complete cases for the variables used (outcome, matching covariates,
#        and optional regression covariates).
#     2) Construct a treated vs control dataset for each control-id set.
#     3) Perform full matching (MatchIt::matchit(method = "full")) on cov_match.
#     4) Compute balance diagnostics (mean/max standardized mean differences),
#        effective sample sizes (ESS) and ESS ratio.
#     5) Fit a weighted mixed-effects model with a random intercept by country
#        (CntrName) using lme4::lmer.
#     6) Save a lightweight result object (.qs) containing summary tables and the
#        top model fit for the best-balanced control-id set.
#
# Inputs
#   Files
#     - ../1_buildData/_dataMain.csv
#         Main analysis dataset. Must contain (at minimum) columns:
#         type, aq_id, CntrName, plus any covariates referenced in cov_match/cov_reg
#         and any outcome_var used below.
#     - ../1_buildData/CtrlNoOverlapHYBAS_B.rds
#         List of candidate control aq_id sets (nonOverlapsB). Each element should
#         be a vector of aq_id values used to define the control pool.
#
#   Objects (in-script)
#     - d: data.frame read from _dataMain.csv
#     - nonOverlapsB: list read from CtrlNoOverlapHYBAS_B.rds
#     - pilot_match(): function that runs matching + model fitting for one outcome
#       specification across all control sets in nonOverlapsB.
#
# Outputs
#   Files
#     - resultOut/<lab>.qs
#         One file per MAIN call. Each .qs contains a list with:
#           * label: lab
#           * summary_df: per-control-set diagnostics + estimated treatment effect
#           * summary_counts: counts of significant/positive/negative effects
#           * top_models: list of the top (best-balance) lmer model(s)
#
#   Objects (in-script)
#     - In MAIN: each named object (e.g., IrIntens, Overdraft, ...) stores either:
#         * the res_light list returned by pilot_match(), or
#         * an error object if that call failed (captured by tryCatch).
#
# Error handling
#   - Inside pilot_match(): each control-id set is wrapped in tryCatch so failures
#     do not abort the full batch; the error message is stored per control-id set.
#   - In the MAIN section: each pilot_match() call is wrapped in tryCatch so a
#     single outcome specification failing does not abort the whole script; the
#     error object is retained.
# =============================================================================

library(qs)
library(sf)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(purrr)
library(furrr)
library(future)
library(MatchIt)
library(lme4)
library(lmerTest)     # recommended for p-values
library(broom.mixed)
library(parallel)

future::plan(future::sequential)

dir <- "resultOut"

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
d <- read.csv("../1_buildData/_dataMain.csv")
nonOverlapsB <- readRDS(paste0("../1_buildData/CtrlNoOverlapHYBAS_B.rds"))

# -----------------------------------------------------------------------------
# Core routine: for a given outcome specification, iterate over candidate control
# sets, match + fit, summarize, and save results.
# -----------------------------------------------------------------------------
pilot_match <- function(
    d,                    # full dataset
    cov_match,            # matching covariates
    cov_reg = NULL,       # optional regression covariates
    outcome_var,          # outcome variable name (string)
    nonOverlaps,          # list of control aq_id sets
    n_cores = parallel::detectCores() - 1,
    lab = "label"
) {
  
  # --- 1. Keep only complete cases for all variables used ---
  vars_keep <- unique(c(outcome_var, "type", "aq_id", "CntrName", cov_match, cov_reg))
  vars_keep <- vars_keep[vars_keep %in% names(d)]  # keep only existing columns
  
  d <- d %>%
    dplyr::select(all_of(vars_keep)) %>%
    dplyr::filter(stats::complete.cases(across(everything())))
  
  # Split data
  cndat <- d %>% filter(type != "treat")
  trdat <- d %>% filter(type == "treat")
  
  plan(multisession, workers = n_cores)
  
  # -------------------------------------------------------------------------
  # Match + fit for one candidate control-id set (safe: returns error message)
  # -------------------------------------------------------------------------
  match_and_fit <- function(control_ids) {
    tryCatch({
      dat <- trdat %>%
        bind_rows(cndat %>% filter(aq_id %in% control_ids)) %>%
        mutate(type = as.factor(type))
      
      # Matching
      form_match <- reformulate(cov_match, response = "type")
      m.out <- matchit(form_match, data = dat, method = "full")
      
      # Balance metrics
      s <- summary(m.out, standardize = TRUE)
      smd <- abs(s$sum.matched[, "Std. Mean Diff."])
      mean_smd <- mean(smd, na.rm = TRUE)
      max_smd  <- max(smd, na.rm = TRUE)
      
      # Matched data + combined weights
      mdat <- match.data(m.out)
      mdat <- mdat %>%
        mutate(weights_combined = weights)
      # same-country share among matched controls
      subclass_treat_country <- mdat %>%
        filter(type == "treat") %>%
        group_by(subclass) %>%
        summarise(treat_country = first(CntrName), .groups = "drop")
      
      ctrl <- mdat %>%
        filter(type != "treat") %>%
        left_join(subclass_treat_country, by = "subclass") %>%
        mutate(same_country = (CntrName == treat_country))
      
      same_country_unweighted <- mean(ctrl$same_country, na.rm = TRUE)
      same_country_weighted   <- weighted.mean(
        ctrl$same_country,
        w = ctrl$weights_combined,
        na.rm = TRUE
      )
      # ESS calculations using combined weights
      n_treat   <- sum(mdat$type == "treat")
      n_control <- sum(mdat$type != "treat")
      
      ess_treat <- with(
        subset(mdat, type == "treat"),
        (sum(weights_combined)^2) / sum(weights_combined^2)
      )
      
      ess_control <- with(
        subset(mdat, type != "treat"),
        (sum(weights_combined)^2) / sum(weights_combined^2)
      )
      
      ess_ratio <- ess_control / ess_treat
      
      # Mixed model formula (type always included)
      rhs <- c("type", cov_reg)
      form_lm <- reformulate(rhs, response = outcome_var)
      
      fit <- lmer(
        update(form_lm, . ~ . + (1 | CntrName)),
        data = mdat,
        weights = weights_combined
      )
      
      tidy_fit <- broom.mixed::tidy(fit, effects = "fixed")
      
      treat_row <- tidy_fit %>%
        filter(grepl("^type", term)) %>%
        select(term, estimate, std.error, p.value)
      
      int_row <- tidy_fit %>%
        filter(term == "(Intercept)") %>%
        select(term, estimate, std.error, p.value)
      
      list(
        success      = TRUE,
        mean_smd     = mean_smd,
        max_smd      = max_smd,
        ess_ratio    = ess_ratio,
        treat_eff    = treat_row$estimate,
        treat_se     = treat_row$std.error,
        treat_p      = treat_row$p.value,
        int_eff      = int_row$estimate,
        int_se       = int_row$std.error,
        int_p        = int_row$p.value,
        n_treat      = n_treat,
        n_control    = n_control,
        n_total      = n_treat + n_control,
        match        = m.out,
        same_country_unweighted = same_country_unweighted,
        same_country_weighted   = same_country_weighted,
        fit          = fit
      )
    }, error = function(e) {
      list(success = FALSE, error = conditionMessage(e))
    })
  }
  
  # Run all matches in parallel
  results <- future_map(nonOverlaps, match_and_fit, .progress = TRUE)
  
  # Summarize successful runs
  summary_df <- results %>%
    imap_dfr(~{
      if (.x$success) {
        tibble(
          idx        = .y,
          mean_smd   = .x$mean_smd,
          max_smd    = .x$max_smd,
          ess_ratio  = .x$ess_ratio,
          treat_eff  = .x$treat_eff,
          treat_se   = .x$treat_se,
          treat_p    = .x$treat_p,
          int_eff    = .x$int_eff,
          int_se     = .x$int_se,
          int_p      = .x$int_p,
          n_treat    = .x$n_treat,
          n_control  = .x$n_control,
          same_country_unweighted = .x$same_country_unweighted,
          same_country_weighted   = .x$same_country_weighted,
          n_total    = .x$n_total
        )
      } else {
        tibble(
          idx = .y,
          mean_smd = NA, max_smd = NA,
          ess_ratio = NA,
          treat_eff = NA, treat_se = NA, treat_p = NA,
          int_eff = NA, int_se = NA, int_p = NA,
          n_treat = NA, n_control = NA, n_total = NA,          same_country_unweighted = NA,
          same_country_weighted = NA
        )
      }
    }) %>%
    filter(!is.na(mean_smd), !is.na(max_smd), !is.na(ess_ratio)) %>%
    arrange(mean_smd)
  
  # Significance categories (treatment)
  sig_breaks <- c("<0.01", "<0.05", "<0.1", "≥0.1")
  summary_df <- summary_df %>%
    mutate(
      sig_cat = case_when(
        treat_p < 0.01 ~ "<0.01",
        treat_p < 0.05 ~ "<0.05",
        treat_p < 0.1  ~ "<0.1",
        TRUE           ~ "≥0.1"
      ) %>% factor(levels = sig_breaks)
    )
  
  # Significance categories (intercept)
  sig_breaks_int <- c("<0.01", "<0.05", "<0.1", "≥0.1")
  summary_df <- summary_df %>%
    mutate(
      sig_cat_int = case_when(
        int_p < 0.01 ~ "<0.01",
        int_p < 0.05 ~ "<0.05",
        int_p < 0.1  ~ "<0.1",
        TRUE         ~ "≥0.1"
      ) %>% factor(levels = sig_breaks_int)
    )
  
  # Summary counts
  summary_counts <- summary_df %>%
    mutate(
      signif = treat_p < 0.1,
      sign   = case_when(
        treat_eff > 0  ~ "positive",
        treat_eff < 0  ~ "negative",
        TRUE           ~ "zero"
      )
    ) %>%
    summarise(
      total        = n(),
      n_signif     = sum(signif),
      n_pos        = sum(sign == "positive"),
      n_neg        = sum(sign == "negative"),
      n_pos_signif = sum(signif & sign == "positive"),
      n_neg_signif = sum(signif & sign == "negative")
    )
  
  res <- list(
    summary_df     = summary_df,
    summary_counts = summary_counts,
    results        = results
  )
  
  # Identify top matches by balance quality
  top_idx <- head(res$summary_df$idx, 1)
  top_results <- res$results[top_idx]
  
  # Extract top model fits
  top_models <- lapply(top_results, function(x) x$fit)
  
  # Keep only key elements
  res_light <- list(
    label          = lab,
    summary_df     = res$summary_df,
    summary_counts = res$summary_counts,
    top_models     = top_models
  )
  
  qs::qsave(res_light, file.path(dir, paste0(lab, ".qs")))
  
  return(res_light)
}

# -----------------------------------------------------------------------------
# MAIN
#   Each call is wrapped in tryCatch so failures do not abort the whole script.
#   The error object is retained in the corresponding variable.
# -----------------------------------------------------------------------------
# nonOverlapsB = nonOverlapsB[1:10]

IrIntens <- tryCatch(
  pilot_match(d = d%>%filter(Ir>0),cov_match = c('lat_c','lon_c','RDS_mainDist','CS_max'),cov_reg = NULL,outcome_var = "Ir",nonOverlaps = nonOverlapsB,lab='IrIntens'),
  error = function(e) e
)
IrIntens$summary_df %>%
  slice(1) %>%
  transmute(
    same_country_unweighted = 100 * same_country_unweighted,
    same_country_weighted   = 100 * same_country_weighted
  )
Overdraft <- tryCatch(
  pilot_match(d = d%>%filter(Ir>0),cov_match = c('lat_c','lon_c','RDS_mainDist','CS_max'),cov_reg =NULL,outcome_var = "OverIR3",nonOverlaps = nonOverlapsB,lab='Overdraft'),
  error = function(e) e
)

IrGini <- tryCatch(
  pilot_match(d = d%>%filter(Ir>0),cov_match = c('lat_c','lon_c','RDS_mainDist','giniCSI'),cov_reg = NULL,outcome_var = "giniIr",nonOverlaps = nonOverlapsB,lab='IrGini'),
  error = function(e) e
)

IrRivsGini <- tryCatch(
  pilot_match(d = d%>%filter(Ir>0)%>%mutate(riv=RiverBorder),cov_match = c('lat_c','lon_c','RDS_mainDist','giniCSI'),cov_reg = 'riv',outcome_var = "giniIr",nonOverlaps = nonOverlapsB,lab='IrRivsGini'),
  error = function(e) e
)

GWGini <- tryCatch(
  pilot_match(d = d%>%filter(GW>0),cov_match = c('lat_c','lon_c','RDS_mainDist','giniCSI'),cov_reg = NULL,outcome_var = "giniGw",nonOverlaps = nonOverlapsB,lab='GWGini'),
  error = function(e) e
)

GWRivsGini <- tryCatch(
  pilot_match(d = d%>%filter(GW>0)%>%mutate(riv=RiverBorder),cov_match = c('lat_c','lon_c','RDS_mainDist','giniCSI'),cov_reg = 'riv',outcome_var = "giniGw",nonOverlaps = nonOverlapsB,lab='GWRivsGini'),
  error = function(e) e
)

CropSuitInt <- tryCatch(
  pilot_match(d =  d%>%filter(Ir>0),cov_match = c('lat_c','lon_c','RDS_mainDist'),cov_reg =NULL,outcome_var = "CS_max",nonOverlaps = nonOverlapsB,lab='CropSuitInt'),
  error = function(e) e
)

IrNeedInt <- tryCatch(
  pilot_match(d =  d%>%filter(Ir>0),cov_match = c('lat_c','lon_c','RDS_mainDist'),cov_reg =NULL,outcome_var = "IrNeed3",nonOverlaps = nonOverlapsB,lab='IrNeedInt'),
  error = function(e) e
)

# IrNeed0Int <- tryCatch(
#   pilot_match(d =  d%>%filter(Ir>0),cov_match = c('lat_c','lon_c','RDS_mainDist'),cov_reg =NULL,outcome_var = "IrNeed0",nonOverlaps = nonOverlapsB,lab='IrNeed0Int'),
#   error = function(e) e
# )

Overdraft_Irrig <- tryCatch(
  pilot_match(d = d%>%filter(Ir>0),cov_match = c('lat_c','lon_c','RDS_mainDist','Ir'),cov_reg =NULL,outcome_var = "OverIR3",nonOverlaps = nonOverlapsB,lab='Overdraft_Irrig'),
  error = function(e) e
)

Overdraft6 <- tryCatch(
  pilot_match(d = d%>%filter(Ir>0),cov_match = c('lat_c','lon_c','RDS_mainDist','CS_max'),cov_reg =NULL,outcome_var = "OverIR6",nonOverlaps = nonOverlapsB,lab='Overdraft6'),
  error = function(e) e
)

# Overdraft6Irrig <- tryCatch(
#   pilot_match(d = d%>%filter(Ir>0),cov_match = c('lat_c','lon_c','RDS_mainDist','Ir'),cov_reg =NULL,outcome_var = "OverIR6",nonOverlaps = nonOverlapsB,lab='Overdraft6Irrig'),
#   error = function(e) e
# )

# Ir_RivInt <- tryCatch(
#   pilot_match(d = d%>%filter(Ir>0)%>%mutate(riv=RiverBorder),cov_match = c('lat_c','lon_c','RDS_mainDist','CS_max'),cov_reg = 'riv',outcome_var = "Ir",nonOverlaps = nonOverlapsB,lab='Ir_RivInt'),
#   error = function(e) e
# )

IrRivsInt <- tryCatch(
  pilot_match(d = d%>%filter(Ir>0)%>%mutate(riv=RiverBorder),cov_match = c('lat_c','lon_c','RDS_mainDist','CS_max'),cov_reg = 'riv',outcome_var = "Ir",nonOverlaps = nonOverlapsB,lab='IrRivsInt'),
  error = function(e) e
)
GWInt <- tryCatch(
  pilot_match(d = d%>%filter(GW>0),cov_match = c('lat_c','lon_c','RDS_mainDist','CS_max'),cov_reg = NULL,outcome_var = "GW",nonOverlaps = nonOverlapsB,lab='GWInt'),
  error = function(e) e
)

# GW_RivInt <- tryCatch(
#   pilot_match(d = d%>%filter(GW>0)%>%mutate(riv=RiverBorder),cov_match = c('lat_c','lon_c','RDS_mainDist','CS_max'),cov_reg = 'riv',outcome_var = "GW",nonOverlaps = nonOverlapsB,lab='GW_RivInt'),
#   error = function(e) e
# )

IrIntDoubleRobust <- tryCatch(
  pilot_match(d = d %>% filter(Ir > 0),cov_match = c('lat_c','lon_c','RDS_mainDist','CS_max'),cov_reg = c('RDS_mainDist','CS_max'),outcome_var = "Ir",nonOverlaps = nonOverlapsB,lab = 'ItIntdoublerobust'
  ),
  error = function(e) e
)

# CropSuitGini <- tryCatch(
#   pilot_match(d = d %>% filter(Ir > 0),cov_match = c('lat_c','lon_c','RDS_mainDist'),cov_reg = NULL,outcome_var = "giniCSI",nonOverlaps = nonOverlapsB,lab = 'CropSuitGini'
#   ),
#   error = function(e) e
# )

# IrNeedGini <- tryCatch(
#   pilot_match(d = d %>% filter(Ir > 0),cov_match = c('lat_c','lon_c','RDS_mainDist'),cov_reg = NULL,outcome_var = "giniIr_Need",nonOverlaps = nonOverlapsB,lab = 'IrNeedGini'
#   ),
#   error = function(e) e
# )

# IrNeed0Gini <- tryCatch(pilot_match(d = d %>% filter(Ir > 0),cov_match = c('lat_c','lon_c','RDS_mainDist'),cov_reg = NULL,outcome_var = "gini_crpGWSpct_gt0",nonOverlaps = nonOverlapsB,lab = 'IrNeed0Gini'
#   ),
#   error = function(e) e
# )

# IrGiniFullData <- tryCatch(
#   pilot_match(d = d,cov_match = c('lat_c','lon_c','RDS_mainDist','giniCSI'),cov_reg = NULL,outcome_var = "giniIr",nonOverlaps = nonOverlapsB,lab = 'IrGiniFullData'
#   ),
#   error = function(e) e
# )

IrGiniDoubleRobust <- tryCatch(
  pilot_match( d = d %>% filter(Ir > 0),cov_match = c('lat_c','lon_c','RDS_mainDist','giniCSI'),cov_reg = c('RDS_mainDist','giniCSI'),outcome_var = "giniIr",nonOverlaps = nonOverlapsB,lab = 'IrGiniDoubleRobust'
               ),
  error = function(e) e
)

GWShareGt0 <- tryCatch(
  pilot_match(d = d %>% filter(GW > 0),cov_match = c('lat_c','lon_c','RDS_mainDist','giniCSI'),cov_reg = NULL,outcome_var = "gini_gwpct_gt0",nonOverlaps = nonOverlapsB,lab = 'GWGWfgt0'
  ),
  error = function(e) e
)

# --- missing calls (added) ---

# IrIntFullData <- tryCatch(
#   pilot_match(d = d,cov_match = c('lat_c','lon_c','RDS_mainDist','CS_max'),cov_reg = c('RDS_mainDist','CS_max'),outcome_var = "Ir",nonOverlaps = nonOverlapsB,lab='IrIntFullData'),
#   error = function(e) e
# )


# GW_RivsInt <- tryCatch(
#   pilot_match(d = d%>%filter(GW>0)%>%mutate(riv=RiverBorder),cov_match = c('lat_c','lon_c','RDS_mainDist','CS_max'),cov_reg = 'riv',outcome_var = "GW",nonOverlaps = nonOverlapsB,lab='GWRivsInt'),
#   error = function(e) e
# )

# Overdraft6_Irrig <- tryCatch(
#   pilot_match(d = d%>%filter(Ir>0),cov_match = c('lat_c','lon_c','RDS_mainDist','Ir'),cov_reg = NULL,outcome_var = "OverIR6",nonOverlaps = nonOverlapsB,lab='Overdraft6_Irrig'),
#   error = function(e) e
# )
