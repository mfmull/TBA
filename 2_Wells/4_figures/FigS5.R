# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script visualizes the sensitivity of second-stage estimates to randomly
# dropping k transboundary aquifers (“leave-k-out” robustness). It reads the
# Monte Carlo leave-k-out output (loo_Out.csv), summarizes how the coefficient
# distribution changes with k, and plots (i) the mean effect ± 1 SD and (ii) the
# share of draws that remain statistically significant, using a dual y-axis.
#
# Inputs
# - ../3_secondstage/loo_Out.csv
#     Monte Carlo leave-k-out results containing, for each draw:
#       * k: number of aquifers dropped
#       * iter: repetition index
#       * looaqf: list of dropped aquifers (as a string)
#       * b_0, se_0, p_0, nm_0: reg=0 coefficients/SEs/p-values/names
#       * b_1, se_1, p_1, nm_1: reg=1 coefficients/SEs/p-values/names
#
# Core function: plot_loo_robustness(resultsLOO, reg, coef, alpha_sig)
# - reg selects which outcome to plot:
#     * reg = 0: “Ave Depletion” (mean depletion / intercept-level outcome)
#     * reg = 1: “Dist. Trend” (distance-to-border trend outcome)
# - coef selects which coefficient to plot:
#     * "treatment"  -> term code "TBTRUE"  (TB effect)
#     * "intercept"  -> term code "intrcpt" (Intercept)
# - alpha_sig defines the significance threshold used for the “share significant” line.
#
# Processing steps inside plot_loo_robustness()
# 1) Select the appropriate coefficient columns (b_{reg}, p_{reg}, nm_{reg}).
# 2) Filter the dataset to the requested term (intrcpt or TBTRUE).
# 3) Summarize across Monte Carlo draws for each k:
#    - mean_b: mean coefficient across draws
#    - sd_b:   standard deviation across draws
#    - share_sig: fraction of draws with p-value < alpha_sig
# 4) Create a dual-axis scaling:
#    - share_sig is rescaled to the left-axis range via scale_fac so it can be
#      overlaid; the right axis maps back to the original 0–1 share scale.
#
# Plot elements
# - Blue ribbon: mean_b ± 1 SD across draws at each k
# - Blue line:   mean_b
# - Dashed gray line: y = 0 reference
# - Red dotted line:  share_sig (scaled to overlay), with a right-side axis
# - Faint dotted horizontals: visual guides for secondary-axis breakpoints
#
# Main calls
# - Read loo_Out.csv
# - Generate:
#     p1: reg=1, TB effect
#     p2: reg=1, intercept
#     p3: reg=0, TB effect
# - Combine selected plots with patchwork: p3 + p1 + p2
#
# Output
# - A multi-panel figure showing how the estimated TB effect (and optionally the
#   intercept) shifts and how often it remains significant as more TB aquifers are
#   removed (increasing k).
# ------------------------------------------------------------------------------
plot_loo_robustness <- function(resultsLOO,
                                reg  = 0,
                                coef = c("treatment", "intercept"),
                                alpha_sig = 0.05) {
  library(dplyr)
  library(ggplot2)
  
  coef <- match.arg(coef)
  
  # Choose correct columns
  b_col  <- paste0("b_",  reg)
  p_col  <- paste0("p_",  reg)
  nm_col <- paste0("nm_", reg)
  
  # Filter term
  term_code  <- if (coef == "intercept") "intrcpt" else "TBTRUE"
  term_label <- if (coef == "intercept") "Intercept" else "TB effect"
  
  df_term <- resultsLOO %>%
    filter(.data[[nm_col]] == term_code)
  
  # dplyr::summarize
  summ <- df_term %>%
    group_by(k) %>%
    summarise(
      mean_b    = mean(.data[[b_col]], na.rm = TRUE),
      sd_b      = sd(.data[[b_col]],   na.rm = TRUE),
      share_sig = mean(.data[[p_col]] < alpha_sig, na.rm = TRUE),
      .groups   = "drop"
    )
  
  # Compute scaling for secondary axis (matching visible range)
  max_abs   <- max(abs(summ$mean_b + summ$sd_b), na.rm = TRUE)
  scale_fac <- max_abs / max(summ$share_sig, 1e-8)
  regTm=ifelse(reg==1,'Dist. Trend', 'Ave Depletion')
  # Secondary axis break positions mapped to primary scale
  sec_breaks <- seq(0, 1, 0.25)
  sec_lines  <- sec_breaks * scale_fac
  
  ggplot(summ, aes(x = k)) +
    # Coefficient (mean ± SD)
    geom_ribbon(aes(ymin = mean_b - sd_b, ymax = mean_b + sd_b),
                fill = "skyblue", alpha = 0.3) +
    geom_line(aes(y = mean_b), color = "blue", linewidth = 1.1) +
    
    # Horizontal reference line at y = 0
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    
    # Faint horizontal lines for secondary axis
    geom_hline(yintercept = sec_lines,
               color = "firebrick", linetype = "dotted", alpha = 0.15) +
    
    # Share significant (scaled to same range)
    geom_line(aes(y = share_sig * scale_fac),
              color = "firebrick", linewidth = 1, linetype = "dotted") +
    
    # Dual axis
    scale_y_continuous(
      name = sprintf("%s (mean ± 1 SD)", term_label),
      sec.axis = sec_axis(~ . / scale_fac,
                          name = sprintf("Share significant (p < %.2f)", alpha_sig),
                          breaks = sec_breaks)
    ) +
    
    labs(
      x = "Number of aquifers dropped (k)",
      title = sprintf("%s: %s", regTm,term_label)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor.y = element_blank(),
      axis.title.y.right = element_text(color = "firebrick"),
      axis.text.y.right  = element_text(color = "firebrick"),
      axis.title.y.left  = element_text(color = "blue"),
      axis.text.y.left   = element_text(color = "blue")
    )
}

resultsLOO=read.csv('../3_secondstage/loo_Out.csv')

p1=plot_loo_robustness(resultsLOO, reg = 1, coef = "treatment")
p2=plot_loo_robustness(resultsLOO, reg = 1, coef = "intercept")
p3=plot_loo_robustness(resultsLOO, reg = 0, coef = "treatment")
p3+p1+p2
