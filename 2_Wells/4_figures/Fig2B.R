# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script defines a compact visualization for second-stage meta-regression
# outputs and applies it to the preferred specification. It produces a two-panel
# coefficient plot (one panel per outcome) showing only the Intercept and TB
# effect estimates with 90% confidence intervals and significance stars.
#
# Inputs
# - A data frame `df` containing second-stage output columns for both outcomes:
#     nm_0, b_0, se_0, p_0   (Mean depletion model terms)
#     nm_1, b_1, se_1, p_1   (Distance Trend model terms)
#   In usage, `df` is read from:
#     ../3_secondstage/preferredOut.csv
#
# Core function: regressplot(df)
# 1) Reshape to a long plotting format:
#    - Start from wide second-stage output (both outcomes in one row-set).
#    - pivot_longer over est_0/est_1 so each term appears twice: once for each
#      outcome (“Mean depletion” and “Distance Trend”).
# 2) Attach the correct standard error and p-value by outcome:
#    - For each long row, pick se_0/p_0 if outcome_id == est_0, else se_1/p_1.
# 3) Compute uncertainty and annotations:
#    - 90% CI using ±1.64 × SE
#    - significance stars based on p-value thresholds (0.1, 0.05, 0.01)
#    - component label restricted to the two terms of interest:
#        intrcpt -> “Intercept”
#        otherwise -> “TB effect”
#    - enforce consistent facet ordering via an outcome factor
#
# Plot
# - Horizontal coefficient plot, faceted by outcome (free x-scale per panel):
#    * points: estimates (shape 21, black outline)
#    * error bars: 90% CI
#    * stars: placed slightly to the right of the estimate
#    * vertical reference line at x = 0 (dashed)
# - Styling:
#    * theme_classic, bold facet labels/title
#    * fill mapped to component (Intercept vs TB effect), legend removed
#
# Output
# - A two-panel figure titled “Meta-regression results (Intercept and TB effect only)”
#   for the preferred specification in preferredOut.csv.
# ------------------------------------------------------------------------------
regressplot <- function(df) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  
  # --- Reshape to long form for plotting ---
  plot_df <- df %>%
    transmute(
      term = nm_0,
      est_0 = b_0, se_0 = se_0, p_0 = p_0,
      est_1 = b_1, se_1 = se_1, p_1 = p_1
    ) %>%
    pivot_longer(
      cols = starts_with("est_"),
      names_to = "outcome_id",
      values_to = "est"
    ) %>%
    mutate(
      se = if_else(outcome_id == "est_0", se_0, se_1),
      p  = if_else(outcome_id == "est_0", p_0, p_1),
      outcome = if_else(outcome_id == "est_0", "Mean depletion", "Distance Trend"),
      ci_low  = est - 1.64 * se,
      ci_high = est + 1.64 * se,
      stars = case_when(
        p < 0.01 ~ "***",
        p < 0.05 ~ "**",
        p < 0.1  ~ "*",
        TRUE ~ ""
      ),
      component = ifelse(term == "intrcpt", "Intercept", "TB effect"),
      outcome = factor(outcome, levels = c("Mean depletion", "Distance Trend"))
    )
  
  # --- Plot ---
  ggplot(plot_df, aes(y = component, x = est, fill = component)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                  width = 0.15, linewidth = 0.6) +
    geom_point(size = 3, shape = 21, color = "black") +
    geom_text(aes(label = stars, x = est * 1.08),
              vjust = 0, size = 5) +
    scale_fill_manual(values = c("Intercept" = "white", "TB effect" = "black")) +
    facet_wrap(~ outcome, scales = "free_x") +
    labs(
      x = NULL,
      y = "Estimated effect (90% CI)",
      title = "Meta-regression results (Intercept and TB effect only)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold")
    )+geom_vline(xintercept=0,lty=2)
}


d=read.csv('../3_secondstage/preferredOut.csv')
regressplot(d)
