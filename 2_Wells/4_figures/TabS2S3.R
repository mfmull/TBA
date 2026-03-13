# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script builds publication-ready regression tables (LaTeX by default) that
# compare the preferred second-stage meta-regression results against a set of
# robustness and alternative specifications. It reads pre-computed second-stage
# outputs from CSV files, assigns human-readable labels, and formats them into a
# single table per outcome (reg = 0 for mean depletion; reg = 1 for border-distance
# trend).
#
# Inputs (pre-computed model output CSVs)
# - ../3_secondstage/preferredOut.csv          : preferred matched rma.mv results
# - ../3_secondstage/unMachedOut.csv           : no matching (weights only from SEs)
# - ../3_secondstage/robustOut.csv             : second stage using robust first-stage estimates
# - ../3_secondstage/altCovOut.csv             : matching with alternative covariates
# - ../3_secondstage/dropTopSEOur.csv          : preferred spec after dropping “heavy” aquifers
# - ../3_secondstage/FE_Out.csv                : country fixed-effects meta-reg (rma.uni with factor(CC))
# - ../3_secondstage/clustRobustSE_Out.csv     : preferred spec with cluster-robust (CR2) inference
#
# Core function: make_regtable(model_list, reg, caption, fmt)
# - Expects each element of model_list to be a data frame with a `lab` column and
#   columns named like:
#     nm_{reg}, b_{reg}, se_{reg}, p_{reg}, st{reg}, tau2_CC_{reg}, n, k_cluster
# - Stacks models, then constructs a coefficient block for the Intercept and TB
#   term (Transboundary), displaying two rows per term:
#     * coefficient with standard error:   b (se)
#     * p-value with stars:                [p]***
# - Adds a stats panel with Tau^2, Observations, and Clusters per specification.
# - Uses knitr::kable + kableExtra styling, including a header row with numbered
#   specifications and a midrule separating coefficients from summary statistics.
# - escape = FALSE is used so LaTeX math formatting (e.g., brackets) is preserved.
#
# Table assembly
# 1) Read each robustness/specification CSV and set a short label in `lab`.
# 2) For the fixed-effects variant (FE_Out.csv), keep only the first two rows
#    (Intercept and TB effect) and blank out tau^2 fields (not defined there).
# 3) Combine all model outputs into a list `mods` in the desired column order.
# 4) Call make_regtable() twice:
#    - reg = 0: “Effect of TB status on mean depletion”
#    - reg = 1: “Effect of TB status on distance to border trend”
#    Each call supplies a long caption describing the mapping from columns (1)–(…)
#    to the robustness checks.
#
# Outputs
# - Two formatted tables returned to the R session (and printable/knittable):
#     * reg = 0 table: pooled TB effect on mean depletion (intercept-level outcome)
#     * reg = 1 table: pooled TB effect on the distance-to-border trend coefficient
#   When run inside an RMarkdown/Quarto workflow, these render as LaTeX tables
#   with consistent formatting across specifications.
# ------------------------------------------------------------------------------
make_regtable <- function(model_list, reg = 0, caption = "Regression results", fmt = "latex") {
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(knitr)
  library(kableExtra)
  
  b_col   <- paste0("b_",  reg)
  se_col  <- paste0("se_", reg)
  p_col   <- paste0("p_",  reg)
  nm_col  <- paste0("nm_", reg)
  st_col  <- paste0("st",  reg)
  tau_col <- paste0("tau2_CC_", reg)
  
  need_cols <- c(nm_col, b_col, se_col, p_col, st_col, tau_col, "n", "k_cluster", "lab")
  
  # Combine all list elements
  df <- purrr::map_dfr(model_list, ~ .x %>%
                         dplyr::select(dplyr::all_of(need_cols)))
  
  # Rename for simplicity
  names(df) <- c("term", "b", "se", "p", "stars", "tau2", "n", "k_cluster", "label")
  
  # For coefficients/SEs and p-values
  fmt_num <- function(x) ifelse(is.na(x), "",
                                ifelse(abs(x) < 0.001, "<0.001", sprintf("%.3f", x)))
  fmt_tau <- function(x) ifelse(is.na(x), "", sprintf("%.1f", x))
  
  ##############################
  # Main coefficients: 2 rows per term (coef+se; then p+stars)
  ##############################
  df_main <- df %>%
    dplyr::filter(term %in% c("intrcpt", "TBTRUE")) %>%
    dplyr::mutate(
      term_disp  = dplyr::recode(term, intrcpt = "Intercept", TBTRUE = "Transboundary"),
      entry_coef = sprintf("%s (%s)", fmt_num(b), fmt_num(se)),
      entry_p    = sprintf("$[%s]$%s", fmt_num(p), stars)  # Fix A: math mode, no $<$ token
    ) %>%
    dplyr::select(term, term_disp, label, entry_coef, entry_p) %>%
    tidyr::pivot_longer(
      cols = c(entry_coef, entry_p),
      names_to = "rowtype",
      values_to = "entry"
    ) %>%
    dplyr::mutate(
      row_id   = paste0(term, "_", rowtype),                 # ensures uniqueness in pivot_wider
      term_out = ifelse(rowtype == "entry_p", "", term_disp) # blank displayed term on p-value row
    ) %>%
    dplyr::select(row_id, term_out, label, entry) %>%
    tidyr::pivot_wider(
      id_cols = c(row_id, term_out),
      names_from = label,
      values_from = entry
    ) %>%
    dplyr::arrange(row_id) %>%
    dplyr::select(-row_id) %>%
    dplyr::rename(term = term_out)
  
  ##############################
  # Stats panel
  ##############################
  df_stats <- df %>%
    dplyr::distinct(label, tau2, n, k_cluster) %>%
    dplyr::mutate(
      tau2 = fmt_tau(tau2),
      n = as.character(n),
      k_cluster = as.character(k_cluster)
    ) %>%
    tidyr::pivot_longer(
      cols = c(tau2, n, k_cluster),
      names_to = "term",
      values_to = "entry"
    ) %>%
    dplyr::mutate(
      term = dplyr::recode(term, tau2 = "Tau^2", n = "Observations", k_cluster = "Clusters")
    ) %>%
    tidyr::pivot_wider(names_from = label, values_from = entry)
  
  df_main  <- df_main  %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  df_stats <- df_stats %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  
  df_out <- dplyr::bind_rows(df_main, df_stats)
  
  model_labels <- names(df_out)[-1]
  spec_nums <- paste0("(", seq_along(model_labels), ")")
  
  knitr::kable(
    df_out, format = fmt, booktabs = TRUE, escape = FALSE,
    caption = caption, align = "lcccc",
    col.names = c("", model_labels)
  ) %>%
    kableExtra::add_header_above(
      c(" " = 1, stats::setNames(rep(1, length(model_labels)), spec_nums))
    ) %>%
    kableExtra::kable_styling(latex_options = c("hold_position", "scale_down")) %>%
    # midrule between coefficient block and stats block
    kableExtra::row_spec(nrow(df_main), extra_latex_after = "\\midrule")
}



a2=read.csv('../3_secondstage/preferredOut.csv')
a2$lab='Preferred'


##Regression table with robustness tests
#Unmatched
a1=read.csv('../3_secondstage/unMachedOut.csv')
a1$lab='UnMatched'

#Robust lm #barely ns.
a3=read.csv('../3_secondstage/robustOut.csv')
a3$lab='Robust LM'
#alt covs
a4=read.csv('../3_secondstage/altCovOut.csv')
a4$lab='Alt Cov'


a6=read.csv('../3_secondstage/dropTopSEOur.csv')
a6$lab='Drop Heavy'

a8=read.csv('../3_secondstage/FE_Out.csv')[1:2,]
a8$lab='Fixed Effects'
a8$tau2_CC_0<-a8$tau2_CC_0<-NA
a8$tau2_CC_1<-a8$tau2_CC_1<-NA


a9=read.csv('../3_secondstage/clustRobustSE_Out.csv')
a9$lab='Cluster-Robust SE (CR2)'

mods=list(a1,a2,a3,a4,a6,a8,a9)
make_regtable(mods,0,'Effect of TB status on mean depletion. (1) No matching weight; (2) Preferred specification: full matching weight, with LatLon, urban cover, Crop Suitability and Propensity of river border as matching covariates; (3) Aquifer-level robust linear models; (4) Alternative matching covariates:  Soil and terrain suitability, precipitations, potential ET; (5) Double Robust specification where matching covariates are also used as aquifer-level covariates; (6) Exclusion of 10 control aquifers with highest matching weights')
make_regtable(mods,1,'Effect of TB status on distance to border trend. (1) No matching weight; (2) Preferred specification: full matching weight, with LatLon, urban cover, Crop Suitability and Propensity of river border as matching covariates; (3) Aquifer-level robust linear models; (4) Alternative matching covariates:  Soil and terrain suitability, precipitations, potential ET; (5) Double Robust specification where matching covariates are also used as aquifer-level covariates; (6) Exclusion of 10 control aquifers with highest matching weights')

