# =============================================================================
# 56_figures.R -- FIGURES
#
#   figure_main_typology.pdf        Fig. 4: alluvial + FGW/FBW scatter
#   figure_si_typology_bars.pdf     group means with Tukey letters (SI)
#   figure_si_landuse_scale.pdf     mean-normalised land-use extent by class
#   figure_si_mosaic_*.pdf          class x cooperation / distance / river
#   figure_si_robustness.pdf        class counts across the five variants
#
# Every panel is drawn from the cached objects; nothing is reclassified here.
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!exists("CFG")) source(file.path(.root, "51_config.R"))
if (!exists("classify_dyads")) source(file.path(CFG$root, "52_core.R"))
if (is.null(.log$path)) log_open(CFG$log_file)
suppressPackageStartupMessages(library(patchwork))
cfg <- CFG
t0 <- Sys.time()

say("\n================ FIGURES =====================\n")
pdf_device_note()
main <- readRDS(file.path(cfg$cache_dir, "main_objects.rds"))
rob  <- readRDS(file.path(cfg$cache_dir, "robustness_objects.rds"))
tp <- main$typology; tuk <- main$tukey
cols <- cfg$class_cols
labs_metric <- metric_labels(cfg)
dir.create(cfg$fig_dir, showWarnings = FALSE, recursive = TRUE)

theme_paper <- theme_minimal(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.5, colour = "grey30"),
        legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 7))

# =============================================================================
# A. Alluvial: current -> future class
# =============================================================================
# ggalluvial exports StatStratum for exactly this use. If a future version
# stops doing so, attach the package instead of failing the whole figure.
if (have_alluvial && !exists("StatStratum", where = asNamespace("ggalluvial"),
                             inherits = FALSE)) {
  suppressPackageStartupMessages(library(ggalluvial))
  say("  ggalluvial attached (StatStratum not exported by this version).\n")
}
p_all <- if (have_alluvial) {
  ggplot(tp, aes(axis1 = class, axis2 = classF)) +
    ggalluvial::geom_alluvium(aes(fill = class), width = 1/12, alpha = 0.85) +
    ggalluvial::geom_stratum(aes(fill = class), width = 1/6, colour = "black",
                             linewidth = 0.2) +
    # stat = ggproto object, NOT the string "stratum": the string form makes
    # ggplot2 look up StatStratum on the search path, which only works if
    # ggalluvial has been ATTACHED with library(). This file reaches the
    # package through ggalluvial:: (guarded by have_alluvial), so the object
    # is passed directly and no attach is needed.
    geom_text(stat = ggalluvial::StatStratum,
              aes(label = after_stat(stratum)), size = 2.4) +
    scale_x_discrete(limits = c("Current", "Future"), expand = c(.05, .05)) +
    scale_fill_manual(values = cols, name = NULL, guide = "none") +
    labs(title = "A   Configuration now and under water-stressed cropland",
         x = NULL, y = "Dyads") +
    theme_paper
} else {
  say("  ggalluvial not installed: panel A falls back to a transition bar.\n")
  ggplot(main$typology %>% count(class, classF),
         aes(class, n, fill = classF)) +
    geom_col(width = .7) +
    scale_fill_manual(values = cols, name = "Future") +
    labs(title = "A   Configuration now and under water-stressed cropland",
         x = NULL, y = "Dyads") + theme_paper
}

# =============================================================================
# B. Scatter -- the main Fig. 4B/C panel (reviewer R1.7)
# =============================================================================
# x = national relevance, y = overdraft exposure, one point per dyad, class
# centroids with t-based CI crosshairs. Tukey groupings are shown twice per
# centroid: latin letters below for the x metric, greek to the right for y --
# classes sharing a letter are not significantly different on THAT axis.
# Identity is never colour-alone: every centroid is directly labelled.
greek <- c("alpha", "beta", "gamma", "delta", "epsilon")
letters_x <- tuk %>% filter(Variable == "FGW") %>% select(class, lx = letters)
letters_y <- tuk %>% filter(Variable == "FBW") %>% select(class, ly = letters)
# Map the y-axis letter set onto greek symbols so the two axes cannot be confused.
ly_map <- sort(unique(unlist(strsplit(paste(letters_y$ly, collapse = ""), ""))))
to_greek <- function(s) {
  paste(vapply(strsplit(s, "")[[1]],
               function(ch) greek[match(ch, ly_map)], character(1)),
        collapse = "")
}
cent <- tuk %>%
  select(Variable, class, mean, ci_low, ci_high) %>%
  pivot_wider(names_from = Variable,
              values_from = c(mean, ci_low, ci_high)) %>%
  left_join(letters_x, by = "class") %>%
  left_join(letters_y, by = "class") %>%
  mutate(ly_greek = vapply(ly, to_greek, character(1)))

# The x axis is log10 and FGW is exactly 0 for dyads with no groundwater
# irrigation anywhere in either country -- log10(0) is -Inf, and ggplot drops
# those points with only a warning. Exclude them explicitly and say how many,
# so the panel's sample is stated rather than silently reduced. They are
# almost all interior dyads, where absence of irrigation is the defining
# feature, so the scatter is not the place they are characterised.
pts_all <- tp %>% filter(land_total > 0, is.finite(FGW), is.finite(FBW))
pts <- pts_all %>% filter(FGW > 0)
n_zero <- nrow(pts_all) - nrow(pts)
if (n_zero > 0)
  sayf("  scatter: %d of %d dyads have FGW = 0 and cannot be shown on a log axis (%s).\n",
       n_zero, nrow(pts_all),
       paste(sprintf("%s=%d", names(table(droplevels(pts_all$class[pts_all$FGW <= 0]))),
                     table(droplevels(pts_all$class[pts_all$FGW <= 0]))),
             collapse = ", "))
p_sc <- ggplot(pts, aes(FGW, FBW)) +
  geom_point(aes(colour = class), size = 1.1, alpha = 0.45, stroke = 0) +
  # orientation = "y" rather than geom_errorbarh(), which ggplot2 4.0.0
  # deprecated (and which was silently translating `height` to `width`).
  geom_errorbar(data = cent, inherit.aes = FALSE,
                aes(y = mean_FBW, xmin = ci_low_FGW, xmax = ci_high_FGW,
                    colour = class), orientation = "y", width = 0,
                linewidth = 0.5) +
  geom_errorbar(data = cent, inherit.aes = FALSE,
                aes(x = mean_FGW, ymin = ci_low_FBW, ymax = ci_high_FBW,
                    colour = class), width = 0, linewidth = 0.5) +
  geom_point(data = cent, inherit.aes = FALSE,
             aes(mean_FGW, mean_FBW, fill = class),
             shape = 21, size = 3, stroke = 0.6, colour = "white") +
  geom_text(data = cent, inherit.aes = FALSE,
            aes(mean_FGW, mean_FBW, label = paste0(class, " ", lx, "/",
                                                   ly_greek)),
            vjust = -1.4, size = 2.3, colour = "grey15") +
  scale_x_log10() +
  scale_colour_manual(values = cols, name = NULL) +
  scale_fill_manual(values = cols, guide = "none") +
  labs(title = "B   National relevance against overdraft exposure",
       subtitle = paste("Points are dyads; filled markers are class means with",
                        "95% confidence intervals on both axes.\nLetters are",
                        "Tukey groupings: latin for the x metric, greek for y;",
                        "classes sharing a letter\ndo not differ on that axis."),
       x = labs_metric[["FGW"]], y = labs_metric[["FBW"]]) +
  theme_paper

fig4 <- if (identical(cfg$fig4_main_panel, "scatter"))
  p_all / p_sc + plot_layout(heights = c(1, 1.6)) else p_all
# "Some strata appear at multiple axes" is expected and correct here: the two
# axes are the SAME five classes at two points in time, which is the whole
# point of the panel. Muffle only that message -- fired by to_lodes_form()
# during ggplot_build, so it has to be caught at the draw call, not at
# construction -- and let every other warning through.
withCallingHandlers(
  save_pdf(fig4, "figure_main_typology", 7.2, 8.0, cfg),
  warning = function(w) {
    if (grepl("strata appear at multiple axes", conditionMessage(w)))
      invokeRestart("muffleWarning")
  })

# =============================================================================
# SI: the bar panels, so the Tukey letters stay reported in a compact form
# =============================================================================
# NOTE: the legacy SI scripts named the first metric FIR while the factor
# levels were c("FGW","FBW"), which turned every row of that facet into NA and
# printed an "NA" strip. The metric names are fixed in 52_core.R and only the
# axis label varies, so that cannot recur.
bars <- tuk %>% mutate(Variable = factor(Variable, levels = c("FGW", "FBW"),
                                         labels = labs_metric))
p_bars <- ggplot(bars, aes(class, mean, fill = class)) +
  geom_col(width = .6) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .12,
                linewidth = 0.3) +
  geom_text(aes(label = letters, y = ci_high + 0.08 * (ci_high - ci_low)),
            size = 2.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~ Variable, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = cols, name = NULL) +
  labs(x = NULL, y = "Class mean",
       title = "Group means with Tukey groupings",
       subtitle = sprintf(paste("Error bars are standard errors; letters are",
                                "Tukey HSD groupings at %.0f%%.\nClasses",
                                "sharing a letter are not significantly",
                                "different."), 100 * cfg$tukey_conf)) +
  theme_paper +
  theme(strip.text = element_text(face = "bold", size = 7, hjust = 0))
save_pdf(p_bars, "figure_si_typology_bars", 7.0, 3.4, cfg)

# =============================================================================
# SI: how big are these systems? (mean-normalised extent, log scale)
# =============================================================================
lgw <- tp %>% filter(land_total > 0) %>%
  mutate(L = land_total / mean(land_total, na.rm = TRUE)) %>%
  group_by(class) %>%
  summarise(mean = mean(L, na.rm = TRUE), se = se_mean(L, cfg$se_drop_na),
            n = sum(!is.na(L)), .groups = "drop")
p_lgw <- ggplot(lgw, aes(class, mean, fill = class)) +
  geom_col(width = .6) +
  geom_errorbar(aes(ymin = pmax(mean - se, 1e-6), ymax = mean + se),
                width = .12, linewidth = 0.3) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_y_log10() +
  labs(x = NULL, y = "Mean-normalised extent",
       title = "Scale of the systems in each class") +
  theme_paper
save_pdf(p_lgw, "figure_si_landuse_scale", 3.2, 3.2, cfg)

# =============================================================================
# SI: mosaics (base graphics via vcd)
# =============================================================================
if (have_vcd) {
  for (nm in names(main$mosaics)) {
    m <- main$mosaics[[nm]]
    save_pdf_base(
      vcd::mosaic(m$table, shade = TRUE, legend = TRUE,
                  labeling = vcd::labeling_values,
                  main = sprintf("%s (%s p = %.3g)", m$title,
                                 if (isTRUE(cfg$chisq_simulate))
                                   "Monte Carlo" else "chi-square",
                                 m$test$p_value),
                  labeling_args = list(rot_labels = c(0, 0),
                                       gp_labels = grid::gpar(cex = 0.6))),
      paste0("figure_si_mosaic_", nm), 6.5, 5.2, cfg)
  }
} else say("  vcd not installed: mosaic panels skipped.\n")

# =============================================================================
# SI: robustness -- class composition across the registry
# =============================================================================
comp <- bind_rows(lapply(names(rob$variants), function(nm) {
  r <- rob$variants[[nm]]
  tibble(variant = sprintf("%s (%s)", r$label, r$si),
         class = factor(cfg$class_levels, levels = cfg$class_levels),
         n = r$counts$current)
}))
comp <- bind_rows(
  tibble(variant = "Preferred specification (Fig. 4)",
         class = factor(cfg$class_levels, levels = cfg$class_levels),
         n = as.integer(cfg$reference_class[cfg$class_levels])),
  comp)
comp$variant <- factor(comp$variant, levels = unique(comp$variant))
p_rob <- ggplot(comp, aes(n, variant, fill = class)) +
  geom_col(width = .68) +
  scale_fill_manual(values = cols, name = NULL) +
  labs(x = "Dyads", y = NULL,
       title = "Classification under alternative thresholds and zones",
       subtitle = paste("The interior / urban / bilateral partition is",
                        "insensitive; the BL-DA boundary inside the bilateral\n",
                        "group moves with the detection and concentration",
                        "thresholds, as a tail-defined category should.")) +
  theme_paper +
  theme(axis.text.y = element_text(size = 6))
save_pdf(p_rob, "figure_si_robustness", 8.0, 3.6, cfg)

sayf("\n56_figures.R done in %.2f min.\n",
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
