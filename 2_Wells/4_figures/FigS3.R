# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script visualizes sensitivity of the pooled second-stage estimates to the
# nearest-neighbor matching ratio (number of control units matched per treated
# unit). It reads the ratio-sensitivity outputs (matchArchOut.csv), reshapes the
# results into a long plotting dataset for two outcomes, computes 90% confidence
# intervals, and produces a faceted plot showing how the Intercept and TB effect
# change as the matching ratio increases.
#
# Inputs
# - ../3_secondstage/matchArchOut.csv
#     Second-stage results evaluated across matching ratios, containing:
#       * ratio: matching ratio (controls per treated)
#       * b_0, se_0, nm_0: coefficients/SEs/names for the mean depletion model
#       * b_1, se_1, nm_1: coefficients/SEs/names for the distance-trend model
#
# Processing steps
# 1) Reshape results for plotting:
#    - Stack reg=0 (“Mean depletion”) and reg=1 (“Distance Trend”) outputs into
#      one long data frame with columns:
#        ratio, term (coefficient name), model (outcome label), est, se
# 2) Restrict to the coefficients of interest:
#    - Keep only “intrcpt” (Intercept) and “TBTRUE” (TB effect),
#      and relabel intrcpt for display.
# 3) Construct confidence intervals:
#    - Compute 90% CI bands using ±1.64 × SE:
#        ci_low  = est - 1.64 * se
#        ci_high = est + 1.64 * se
# 4) Set facet ordering:
#    - Make `model` a factor so facets appear in a consistent order:
#      Mean depletion first, then Distance Trend.
#
# Plot
# - X-axis: matching ratio
# - Y-axis: estimated coefficient
# - Lines/points: point estimates by ratio
# - Ribbon: 90% confidence band
# - Dashed horizontal line at 0 as a reference
# - Facets: separate panels for Mean depletion and Distance Trend (free y-scales)
# - Colors/fills: distinguish Intercept vs TBTRUE consistently across panels
#
# Output
# - A two-panel (faceted) figure showing how pooled estimates for:
#    * Intercept and TB effect in the mean depletion model, and
#    * Intercept and TB effect in the distance-trend model
#   vary with the chosen nearest-neighbor matching ratio.
# ------------------------------------------------------------------------------
results_ratio=read.csv('../3_secondstage/matchArchOut.csv')
plot_df <- bind_rows(
  # Regression 0: "Mean depletion"
  results_ratio %>%
    transmute(
      ratio,
      term  = nm_0,
      model = "Mean depletion",
      est   = b_0,
      se    = se_0
    ),
  # Regression 1: "Distance Trend"
  results_ratio %>%
    transmute(
      ratio,
      term  = nm_1,
      model = "Distance Trend",
      est   = b_1,
      se    = se_1
    )
) %>%
  filter(term %in% c("intrcpt", "TBTRUE")) %>%
  mutate(
    term   = recode(term, intrcpt = "(Intercept)", TBTRUE = "TBTRUE"),
    ci_low = est - 1.64 * se,
    ci_high = est + 1.64 * se
  )%>%mutate(model=factor(model,levels=c('Mean depletion', 'Distance Trend')))

ggplot(plot_df, aes(x = ratio, y = est, color = term)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = term),
              alpha = 0.25, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ model, scales = "free_y") +
  scale_color_manual(values = c("(Intercept)" = "darkorange", "TBTRUE" = "#1f78b4")) +
  scale_fill_manual(values = c("(Intercept)" = "darkorange", "TBTRUE" = "#1f78b4")) +
  labs(
    x = "Matching ratio (controls per treated)",
    y = "Estimated effect"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank())

