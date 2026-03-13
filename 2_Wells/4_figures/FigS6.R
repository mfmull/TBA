# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script visualizes how the pooled second-stage estimates change as a function
# of the spatial declustering radius (dclust). It reads the sensitivity outputs
# (dclustOut.csv), reshapes them into a long plotting dataset for two outcomes,
# computes 90% confidence intervals, and produces two side-by-side line plots for
# the Intercept and TB effect across declustering distances.
#
# Inputs
# - ../3_secondstage/dclustOut.csv
#     Second-stage results by declustering setting, containing (for reg=0 and reg=1):
#       * dclust: declustering radius (km)
#       * b_0, se_0, nm_0: coefficients/SEs/names for the mean depletion model
#       * b_1, se_1, nm_1: coefficients/SEs/names for the distance-trend model
#
# Processing steps
# 1) Reshape results for plotting:
#    - Build a long data frame with columns:
#        dclust, term (coefficient name), model (outcome label), est, se
#      by stacking the reg=0 and reg=1 outputs.
# 2) Construct confidence intervals:
#    - Compute 90% CI bands using ±1.64 × SE:
#        ci_low = est - 1.64 * se
#        ci_high = est + 1.64 * se
# 3) Restrict to the coefficients of interest:
#    - Keep only “intrcpt” (Intercept) and “TBTRUE” (TB effect),
#      and relabel intrcpt for display.
# 4) Define outcome-specific y-axis windows:
#    - Use fixed y-limits per panel (Mean depletion vs Distance Trend) so the
#      visual scale is stable across dclust values.
#
# Plotting
# - make_plot(model_name) generates a ggplot for a given outcome:
#    * ribbon: 90% CI band over dclust
#    * line + points: point estimates over dclust
#    * dashed horizontal zero line as a reference
#    * separate color/fill by term (Intercept vs TBTRUE)
#    * coord_cartesian enforces the predefined y-axis window
# - Two panels are produced:
#    * p1: Mean depletion (legend suppressed)
#    * p2: Distance Trend (legend shown)
# - patchwork combines them side-by-side: p1 + p2
#
# Output
# - A two-panel figure showing sensitivity of:
#    * pooled mean depletion level (Intercept) and TB effect
#    * pooled distance-to-border trend (slope) and TB effect
#   across declustering distances (km), with 90% confidence bands.
# ------------------------------------------------------------------------------
results_dclust=read.csv('../3_secondstage/dclustOut.csv')

library(dplyr)
library(ggplot2)
library(patchwork)

# --- Prepare data from numeric tibble ---
plot_df <- results_dclust %>%
  transmute(
    dclust,
    term = nm_0,
    model = "Mean depletion",
    est = b_0,
    se = se_0
  ) %>%
  bind_rows(
    results_dclust %>%
      transmute(
        dclust,
        term = nm_1,
        model = "Distance Trend",
        est = b_1,
        se = se_1
      )
  ) %>%
  mutate(
    ci_low = est - 1.64 * se,
    ci_high = est + 1.64 * se
  ) %>%
  filter(term %in% c("intrcpt", "TBTRUE")) %>%
  mutate(
    term = recode(term, intrcpt = "(Intercept)", TBTRUE = "TBTRUE")
  )

# --- y-axis limits (adjustable) ---
ylims <- data.frame(
  model = c("Distance Trend", "Mean depletion"),
  ymin = c(-1, -220),
  ymax = c(1, 200)
)

# --- Plotting function ---
make_plot <- function(model_name) {
  lims <- ylims %>% filter(model == model_name)
  
  ggplot(filter(plot_df, model == model_name),
         aes(x = dclust, y = est, color = term, fill = term)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.25, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("(Intercept)" = "darkorange", "TBTRUE" = "#1f78b4")) +
    scale_fill_manual(values = c("(Intercept)" = "darkorange", "TBTRUE" = "#1f78b4")) +
    coord_cartesian(ylim = c(lims$ymin, lims$ymax)) +
    labs(
      x = "declustering distance (km)",
      y = "Estimated effect",
      title = model_name
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.title = element_blank())
}

p1 <- make_plot("Mean depletion") + theme(legend.position = "none")
p2 <- make_plot("Distance Trend")

p1 + p2
