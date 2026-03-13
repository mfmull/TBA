# ==============================================================================
# Script purpose (preamble)
# ==============================================================================
# This script (i) visualizes the overlap-country leave-one-out (LOO) robustness
# results produced in the second-stage analysis, and (ii) summarizes how much of
# the analytic sample’s well coverage and identifying overlap is concentrated in
# the subset of countries that contain BOTH treated (TB) and control (non-TB)
# aquifer–country segments.
#
# Inputs
# - ../3_secondstage/loo_country_overlap_Out.csv
#     Output from the LOO pipeline, containing (per dropped overlap country):
#       * tb_b_0, tb_p_0 : TBTRUE coefficient and p-value for Model 0 (beta_0)
#       * tb_b_1, tb_p_1 : TBTRUE coefficient and p-value for Model 1 (beta_1)
#       * sig_0, sig_1   : significance flags (at alpha)
#       * wells_in_CC    : total wells in the dropped country (for point size)
#       * base_b_0/base_b_1 (or equivalents): full-sample baseline coefficients
#
# - ../2_firststage/firststageMain.csv
#     Unit-level aquifer×country dataset with well counts (n) and TB + CC labels,
#     used to compute overlap-country diagnostics and concentration metrics.
#
# Outputs
# - A two-panel patchwork figure of TBTRUE coefficients under LOO (saved optionally)
# - Console summaries of overlap-country well coverage and concentration
# ==============================================================================
# Libraries
# ==============================================================================
library(tidyverse)
library(ggplot2)
library(patchwork)

# ==============================================================================
# Settings
# ==============================================================================
alpha <- 0.10  # significance threshold shown in plots

# ==============================================================================
# Load data
# ==============================================================================
resultsLOO <- read.csv("../3_secondstage/loo_country_overlap_Out.csv") %>%
  as_tibble()

aqf <- read.csv("../2_firststage/firststageMain.csv") %>%
  as_tibble()

# Identify overlap countries (both TB and non-TB segments)
tab_cc <- aqf %>% count(CC, TB) %>% pivot_wider(names_from = TB, values_from = n, values_fill = 0)
cc_frame <- tab_cc %>% filter(`TRUE` > 0, `FALSE` > 0) %>% pull(CC)

# If baseline coefficients are not explicitly stored, derive from the stored columns
# (preferred is that resultsLOO already contains base_b_0 and base_b_1)
base_b0 <- if ("base_b_0" %in% names(resultsLOO)) resultsLOO$base_b_0[1] else NA_real_
base_b1 <- if ("base_b_1" %in% names(resultsLOO)) resultsLOO$base_b_1[1] else NA_real_

# ==============================================================================
# Plot maker
# ==============================================================================
make_coeff_plot <- function(df, b_col, sig_col, base_b, title, xl) {
  
  plot_df <- df %>%
    filter(!is.na(.data[[b_col]]), !is.na(wells_in_CC)) %>%
    arrange(abs(.data[[b_col]] - base_b)) %>%
    mutate(
      drop_CC = factor(drop_CC, levels = drop_CC),
      sig = .data[[sig_col]]
    )
  
  ggplot(plot_df, aes(x = .data[[b_col]], y = drop_CC)) +
    geom_point(
      aes(size = wells_in_CC, fill = sig),
      shape = 21, color = "black", alpha = 0.85, stroke = 0.6
    ) +
    geom_vline(xintercept = base_b, linetype = "dashed", color = "red") +
    coord_flip() +
    scale_size_continuous(name = "Wells in dropped country") +
    scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "black"), guide = "none") +
    labs(
      title = title,
      x = xl,
      y = "Dropped country (CC)",
      caption = paste0("Dashed red = full-sample estimate; filled indicates p < ", alpha)
    ) +
    theme_minimal()
}

# ==============================================================================
# Build plots (both models)
# ==============================================================================
p_coef0 <- make_coeff_plot(
  resultsLOO,
  b_col   = "tb_b_0",
  sig_col = "sig_0",
  base_b  = base_b0,
  title   = "LOO (overlap countries)",
  xl      = "TB effect coefficient, Average depletion"
)

p_coef1 <- make_coeff_plot(
  resultsLOO,
  b_col   = "tb_b_1",
  sig_col = "sig_1",
  base_b  = base_b1,
  title   = "",
  xl      = "TB effect coefficient, Distance-to-border trend"
)

p_combined <- (p_coef0 | p_coef1) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(p_combined)

# ggsave("loo_tbtrue_coeffs_shared_legend.pdf", p_combined, width = 12.5, height = 5.5)

# ==============================================================================
# Overlap-country diagnostics
# ==============================================================================
total_wells <- sum(aqf$n, na.rm = TRUE)

overlap_wells <- aqf %>%
  filter(CC %in% cc_frame) %>%
  summarise(w = sum(n, na.rm = TRUE)) %>%
  pull(w)

pct_overlap_wells <- 100 * overlap_wells / total_wells
n_overlap_countries <- length(cc_frame)

n_segments_total <- nrow(aqf)
n_segments_overlap <- nrow(aqf %>% filter(CC %in% cc_frame))
pct_segments_overlap <- 100 * n_segments_overlap / n_segments_total

cat("\n--- Overlap-country diagnostics ---\n")
cat("Overlap countries:", n_overlap_countries, "\n")
cat("Total wells:", total_wells, "\n")
cat("Wells in overlap countries:", overlap_wells, "\n")
cat(sprintf("Share of wells in overlap countries: %.1f%%\n", pct_overlap_wells))
cat(sprintf("Share of segments in overlap countries: %.1f%%\n", pct_segments_overlap))
cat("-----------------------------------\n")

# Concentration of wells within overlap countries
overlap_cc_wells <- aqf %>%
  filter(CC %in% cc_frame) %>%
  group_by(CC) %>%
  summarise(wells = sum(n, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(wells)) %>%
  mutate(
    share = wells / sum(wells),
    cum_share = cumsum(share)
  )

print(overlap_cc_wells)

top1 <- overlap_cc_wells$cum_share[1]
top3 <- overlap_cc_wells$cum_share[min(3, nrow(overlap_cc_wells))]
top5 <- overlap_cc_wells$cum_share[min(5, nrow(overlap_cc_wells))]

cat(sprintf("Top-1 overlap country share of overlap wells: %.1f%%\n", 100 * top1))
cat(sprintf("Top-3 overlap country share of overlap wells: %.1f%%\n", 100 * top3))
cat(sprintf("Top-5 overlap country share of overlap wells: %.1f%%\n", 100 * top5))