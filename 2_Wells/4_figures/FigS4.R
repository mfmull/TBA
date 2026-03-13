# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script visualizes sensitivity of the pooled second-stage estimates to the
# minimum wells-per-unit threshold (nMin) used when constructing Aquifer × Country
# units. It reads the nMin sensitivity outputs (nMinOut.csv), reshapes them into
# a long plotting dataset for two outcomes, computes 90% confidence intervals,
# and produces two side-by-side line plots showing how the Intercept and TB effect
# vary as nMin increases.
#
# Inputs
# - ../3_secondstage/nMinOut.csv
#     Second-stage results evaluated across different nMin thresholds, containing:
#       * nMin: minimum wells-per-unit cutoff used in estimation
#       * b_0, se_0, nm_0: coefficients/SEs/names for the mean depletion model
#       * b_1, se_1, nm_1: coefficients/SEs/names for the distance-trend model
#
# Processing steps
# 1) Reshape results for plotting:
#    - Build a long data frame with columns:
#        nMin, term (coefficient name), model (outcome label), est, se
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
#      visual scale is comparable across nMin values.
#
# Plotting
# - make_plot(model_name) generates a ggplot for a given outcome:
#    * ribbon: 90% CI band over nMin
#    * line + points: point estimates over nMin
#    * dashed horizontal zero line as a reference
#    * separate color/fill by term (Intercept vs TBTRUE)
#    * coord_cartesian enforces the predefined y-axis window
# - Two panels are produced:
#    * p1: Mean depletion (legend suppressed)
#    * p2: Distance Trend (legend shown)
# - patchwork combines them side-by-side: p1 + p2
#
# Output
# - A two-panel figure showing how:
#    * pooled mean depletion level (Intercept) and TB effect
#    * pooled distance-to-border trend (slope) and TB effect
#   change when imposing stricter minimum-observation cutoffs nMin.
# ------------------------------------------------------------------------------
results_nMin=read.csv('../3_secondstage/nMinOut.csv')
library(patchwork)
# --- Prepare data from numeric tibble ---
plot_df <- results_nMin %>%
  transmute(
    nMin,
    term = nm_0,
    model = "Mean depletion",
    est = b_0,
    se = se_0
  ) %>%
  bind_rows(
    results_nMin %>%
      transmute(
        nMin,
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
         aes(x = nMin, y = est, color = term, fill = term)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.25, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("(Intercept)" = "darkorange", "TBTRUE" = "#1f78b4")) +
    scale_fill_manual(values = c("(Intercept)" = "darkorange", "TBTRUE" = "#1f78b4")) +
    coord_cartesian(ylim = c(lims$ymin, lims$ymax)) +
    labs(
      x = "Minimum Obs. cutoff (nMin)",
      y = "Estimated effect",
      title = model_name
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.title = element_blank())
}

# --- Generate plots and combine ---
p1 <- make_plot("Mean depletion") + theme(legend.position = "none")
p2 <- make_plot("Distance Trend")

p1 + p2
