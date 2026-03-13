
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

pilot_match_single <- function(
    d,                    # full dataset
    cov_match,            # matching covariates
    cov_reg = NULL,       # optional regression covariates
    outcome_var,          # outcome variable name (string)
    # control_ids,          # ONE control aq_id set (vector)
    n_cores = parallel::detectCores() - 1,  # kept for compatibility; not used
    lab = "label",
    dir = "."             # where to save .qs
) {
  
  # --- 1) Keep only complete cases for all variables used ---
  vars_keep <- unique(c(outcome_var, "type", "aq_id", "CntrName", cov_match, cov_reg))
  vars_keep <- vars_keep[vars_keep %in% names(d)]
  
  d2 <- d %>%
    dplyr::select(dplyr::all_of(vars_keep)) %>%
    dplyr::filter(stats::complete.cases(dplyr::across(dplyr::everything())))
  
  # Split data
  cndat <- d2 %>% dplyr::filter(type != "treat")
  trdat <- d2 %>% dplyr::filter(type == "treat")
  
  # Build analysis dataset for this single control set
  dat <- d2 %>%
    mutate(type = as.factor(type))
  # dat <- trdat %>%
  #   dplyr::bind_rows(cndat %>% dplyr::filter(aq_id %in% control_ids)) %>%
  #   dplyr::mutate(type = as.factor(type))
  
  # --- 2) Matching ---
  form_match <- reformulate(cov_match, response = "type")
  m.out <- MatchIt::matchit(form_match, data = dat, method = "full")
  
  # Balance metrics
  s <- summary(m.out, standardize = TRUE)
  smd <- abs(s$sum.matched[, "Std. Mean Diff."])
  mean_smd <- mean(smd, na.rm = TRUE)
  max_smd  <- max(smd, na.rm = TRUE)
  
  # Matched data + weights
  mdat <- MatchIt::match.data(m.out) %>%
    dplyr::mutate(weights_combined = weights)
  
  # ESS
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
  
  # --- 3) Mixed model ---
  rhs <- c("type", cov_reg)
  form_lm <- reformulate(rhs, response = outcome_var)
  
  fit <- lme4::lmer(
    update(form_lm, . ~ . + (1 | CntrName)),
    data = mdat,
    weights = weights_combined
  )
  
  tidy_fit <- broom.mixed::tidy(fit, effects = "fixed")
  
  # if p.value missing, compute from t-stat + residual df (or normal approx)
  if (!("p.value" %in% names(tidy_fit))) {
    if (all(c("statistic", "df") %in% names(tidy_fit))) {
      tidy_fit <- tidy_fit %>%
        mutate(p.value = 2 * pt(abs(statistic), df = df, lower.tail = FALSE))
    } else if ("statistic" %in% names(tidy_fit)) {
      tidy_fit <- tidy_fit %>%
        mutate(p.value = 2 * pnorm(abs(statistic), lower.tail = FALSE))
    } else {
      tidy_fit <- tidy_fit %>%
        mutate(p.value = NA_real_)
    }
  }
  
  treat_row <- tidy_fit %>%
    filter(grepl("^type", term)) %>%
    select(any_of(c("term", "estimate", "std.error", "p.value")))
  
  int_row <- tidy_fit %>%
    filter(term == "(Intercept)") %>%
    select(any_of(c("term", "estimate", "std.error", "p.value")))
  
  res <- list(
    label       = lab,
    mean_smd    = mean_smd,
    max_smd     = max_smd,
    ess_ratio   = ess_ratio,
    treat_eff   = treat_row$estimate,
    treat_se    = treat_row$std.error,
    treat_p     = treat_row$p.value,
    int_eff     = int_row$estimate,
    int_se      = int_row$std.error,
    int_p       = int_row$p.value,
    n_treat     = n_treat,
    n_control   = n_control,
    n_total     = n_treat + n_control,
    match       = m.out,
    fit         = fit
  )
  
  # Optional: quick one-row summary tibble
  res$summary_row <- tibble::tibble(
    label      = lab,
    mean_smd   = mean_smd,
    max_smd    = max_smd,
    ess_ratio  = ess_ratio,
    treat_eff  = res$treat_eff,
    treat_se   = res$treat_se,
    treat_p    = res$treat_p,
    int_eff    = res$int_eff,
    int_se     = res$int_se,
    int_p      = res$int_p,
    n_treat    = n_treat,
    n_control  = n_control,
    n_total    = n_treat + n_control
  )
  
  # qs::qsave(res, file.path(dir, paste0(lab, ".qs")))
  res
}


# dir <- "resultOut"

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
d <- read.csv("_dataMain_JJ.csv")%>%rename(CntrName=CntrNm)

IrIntens <-
  pilot_match_single(d = d%>%filter(Ir>0,CntrName!=unique(d$CntrName)[30]),cov_match = c('lat_c','lon_c','RDS_mainDist','CS_max'),cov_reg = NULL,outcome_var = "Ir",lab='IrIntens')
Overdraft <- 
  pilot_match_single(d = d%>%filter(Ir>0),cov_match = c('lat_c','lon_c','RDS_mainDist','CS_max'),cov_reg =NULL,outcome_var = "OverIR3",lab='Overdraft')


IrGini <-
  pilot_match_single(d = d%>%filter(Ir>0),cov_match = c('lat_c','lon_c','RDS_mainDist','giniCSI'),cov_reg = NULL,outcome_var = "giniIr",lab='IrGini')

IrRivsGini <- 
  pilot_match_single(d = d%>%filter(Ir>0)%>%mutate(riv=RiverBorder),cov_match = c('lat_c','lon_c','RDS_mainDist','giniCSI'),cov_reg = 'riv',outcome_var = "giniIr",lab='IrRivsGini')

GWGini <-
  pilot_match_single(d = d%>%filter(GW>0),cov_match = c('lat_c','lon_c','RDS_mainDist','giniCSI'),cov_reg = NULL,outcome_var = "giniGw",lab='GWGini')
GWRivsGini <- 
  pilot_match_single(d = d%>%filter(GW>0)%>%mutate(riv=RiverBorder),cov_match = c('lat_c','lon_c','RDS_mainDist','giniCSI'),cov_reg = 'riv',outcome_var = "giniGw",lab='GWRivsGini')
IrRivsInt <- 
  pilot_match_single(d = d%>%filter(Ir>0)%>%mutate(riv=RiverBorder),cov_match = c('lat_c','lon_c','RDS_mainDist','CS_max'),cov_reg = 'riv',outcome_var = "Ir",lab='IrRivsInt')



##############Table
generate_res_table_single <- function(
    res_list,
    output = "latex",
    stars = TRUE,
    labels = NULL,
    caption = NULL,
    gof_omit = "AIC|BIC|Log.Lik|Deviance"
) {
  library(modelsummary)
  library(purrr)
  
  # Extract fitted models (single-run objects store model in $fit)
  models <- map(res_list, "fit")
  
  # Column labels
  if (is.null(labels)) {
    names(models) <- map_chr(res_list, ~ .x$label %||% NA_character_)
    # fallback if label missing
    if (any(is.na(names(models)))) names(models) <- paste0("Model_", seq_along(models))
  } else {
    if (length(labels) != length(models)) {
      stop("Length of 'labels' must match number of models in 'res_list'.")
    }
    names(models) <- labels
  }
  
  # Ensure LaTeX output is tabular-based (avoid talltblr)
  options(modelsummary_factory_latex = "kableExtra")
  
  modelsummary(
    models,
    stars    = stars,
    gof_omit = gof_omit,
    title    = caption,
    output   = output
  )
}

# helper for NULL-coalescing without rlang
`%||%` <- function(x, y) if (is.null(x)) y else x


x <- list(IrIntens, Overdraft, IrRivsInt, IrGini, IrRivsGini, GWGini, GWRivsGini)

generate_res_table_single(
  x,
  output = "latex",
  stars = TRUE,
  labels = sapply(x, `[[`, "label"),
  caption = NULL
)
