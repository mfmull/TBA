# =============================================================================
# 36_figures.R -- CANONICAL FIGURES (PDF ONLY)
#
#   output/figure_main_three_panel.pdf      Panels A, B, C
#   output/figure_si_robustness_h1_h2.pdf   forest, H1 and H2 registry
#   output/figure_si_localization_robustness.pdf
#   output/figure_si_influence.pdf
#
# All values plotted here are read from the same cached objects that 37 writes
# to tables; 38_run_all.R asserts figure-table agreement on the main results.
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!exists("CFG")) source(file.path(.root, "31_config.R"))
if (!exists("run_specification")) source(file.path(CFG$root, "32_core.R"))
if (is.null(.log$path)) log_open(CFG$log_file)
suppressPackageStartupMessages(library(patchwork))
cfg <- CFG
t0 <- Sys.time()

main <- readRDS(file.path(cfg$cache_dir, "main_objects.rds"))
rob  <- readRDS(file.path(cfg$cache_dir, "robustness_objects.rds"))
summ <- readRDS(file.path(cfg$cache_dir, "summaries.rds"))
pref <- main$pref; d <- pref$data; h3 <- main$h3; ident <- main$ident

dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
# PDF device. cairo_pdf embeds fonts (which journals require) but needs a
# working cairo build -- on macOS, XQuartz. Without it grDevices::cairo_pdf
# warns "failed to load cairo DLL" and writes NOTHING, so the figures on disk
# silently stay at whatever the last successful run produced while 38's
# manifest check -- which only tests existence -- still passes. Probe cairo
# once, fall back to base pdf(), and verify every file after writing it.
.cairo_ok <- local({
  f <- tempfile(fileext = ".pdf")
  ok <- tryCatch({
    suppressWarnings(grDevices::cairo_pdf(f, width = 1, height = 1))
    plot.new(); grDevices::dev.off()
    file.exists(f) && file.size(f) > 0
  }, error = function(e) FALSE)
  while (grDevices::dev.cur() > 1) grDevices::dev.off()
  unlink(f)
  isTRUE(ok)
})
if (.cairo_ok) {
  say("  PDF device: cairo_pdf (fonts embedded).\n")
} else {
  say("  PDF device: base pdf() -- cairo is unavailable in this R build.\n",
      "    Figures are still written, but fonts are NOT embedded; install\n",
      "    XQuartz (https://www.xquartz.org) and restart R before submission.\n")
  notes_add("cairo unavailable: figures written with base pdf(), fonts not embedded")
}
save_pdf <- function(p, name, w, h) {
  path <- file.path(cfg$out_dir, paste0(name, ".pdf"))
  if (file.exists(path)) unlink(path)   # so a failed write cannot leave a stale file
  if (.cairo_ok) {
    ggsave(path, p, width = w, height = h, device = grDevices::cairo_pdf)
  } else {
    ggsave(path, p, width = w, height = h, device = grDevices::pdf,
           encoding = "WinAnsi.enc")
  }
  if (!file.exists(path) || file.size(path) < 1000)
    stop("save_pdf: ", basename(path), " was not written (device failure). ",
         "Check the PDF device warnings above.")
  say("  wrote ", basename(path), "\n")
}

theme_paper <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        axis.line  = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        strip.text = element_text(face = "bold", size = 8, hjust = 0),
        plot.title = element_text(face = "bold", size = 8.5),
        plot.subtitle = element_text(size = 6.5, colour = "grey30"),
        plot.tag = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 7))

# ---- colour theme for the main three-panel figure only ("aquifer": teal for
# the transboundary contrast / primary quantities, amber for baseline and
# comparison quantities -- a colour-blind-safe blue/orange-family pair). SI
# figures keep the grey theme_paper above unchanged.
pal <- c(teal_dark = "#0B4F5C", teal = "#1B7A8C", teal_light = "#8FC4CC",
        amber_dark = "#8A5A12", amber = "#D98E04", amber_light = "#F0C46B",
        grey_ref = "grey55")

# =============================================================================
# PANEL A -- H1 and H2 main estimates
# =============================================================================
# Both intervals shown for the transboundary contrast are null-imposed wild
# cluster bootstrap (country, Rademacher weights, test inversion): a thick
# one-sided 95% bound (the pre-stated directional test at alpha_one = 0.05,
# inverting the SAME directional bootstrap p reported as p1) and a thin
# two-sided 90% interval. Unlike a parametric CR2/Wald interval, the wild
# bootstrap bound and p1 cannot disagree about significance by construction.
# The domestic baseline / placebo rows are descriptive, not the coefficient
# the bootstrap targets, and keep CR2-clustered inference.
rowsA <- local({
  s1 <- pref$h1; s2 <- pref$h2
  bind_rows(
    tibble(hyp = "H1: depletion level (mm/yr)", term = "transboundary contrast",
           est = s1$estimate, lo = s1$wcb_lo, hi = s1$wcb_hi,
           bound = s1$wcb_bound_one,
           alt = "greater", is_tb = TRUE,
           note = sprintf("one-sided p = %s", pstar(s1$p_one))),
    tibble(hyp = "H1: depletion level (mm/yr)", term = "matched domestic baseline",
           est = s1$baseline, lo = s1$baseline_lo, hi = s1$baseline_hi,
           bound = NA_real_,
           alt = NA_character_, is_tb = FALSE,
           note = "descriptive baseline"),
    tibble(hyp = "H2: borderward gradient (Fisher z)", term = "transboundary contrast",
           est = s2$estimate, lo = s2$wcb_lo, hi = s2$wcb_hi,
           bound = s2$wcb_bound_one,
           alt = "less", is_tb = TRUE,
           note = sprintf("one-sided p = %s", pstar(s2$p_one))),
    tibble(hyp = "H2: borderward gradient (Fisher z)", term = "intercept (placebo-type)",
           est = s2$baseline, lo = s2$baseline_lo, hi = s2$baseline_hi,
           bound = NA_real_,
           alt = NA_character_, is_tb = FALSE,
           note = sprintf("two-sided p = %s", pstar(s2$baseline_p)))) %>%
    mutate(term = factor(term, levels = c("transboundary contrast",
                                          "matched domestic baseline",
                                          "intercept (placebo-type)")),
           y = 2 - (term != "transboundary contrast") * 1,
           bdir  = case_when(alt %in% "less" ~ "left",
                             alt %in% "greater" ~ "right",
                             TRUE ~ NA_character_))
})
ar <- grid::arrow(length = unit(0.05, "in"), type = "closed")
pA <- ggplot(rowsA, aes(est, y)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  # thin, pale two-sided 90% wild cluster bootstrap interval
  geom_errorbar(aes(xmin = lo, xmax = hi, colour = is_tb), orientation = "y",
                width = 0, linewidth = 0.45, alpha = 0.55) +
  # thick one-sided 95% wild cluster bootstrap bound, arrow open toward the
  # predicted direction (transboundary-contrast rows only)
  geom_segment(data = rowsA %>% filter(bdir == "left"),
               aes(x = bound, xend = -Inf, y = y - 0.16, yend = y - 0.16),
               linewidth = 1.1, colour = pal[["teal_dark"]], arrow = ar) +
  geom_segment(data = rowsA %>% filter(bdir == "right"),
               aes(x = bound, xend = Inf, y = y - 0.16, yend = y - 0.16),
               linewidth = 1.1, colour = pal[["teal_dark"]], arrow = ar) +
  geom_point(data = rowsA %>% filter(is.finite(bound)),
             aes(x = bound, y = y - 0.16), shape = 124, size = 2.4,
             colour = pal[["teal_dark"]]) +
  geom_point(aes(shape = is_tb, size = is_tb, colour = is_tb)) +
  geom_text(aes(label = note), vjust = -1.6, size = 1.9, colour = "grey25") +
  geom_text(aes(x = -Inf, y = y + 0.34, label = term), hjust = -0.02,
            size = 2.0, colour = "grey40", fontface = "italic") +
  scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 2.8, `FALSE` = 1.7), guide = "none") +
  scale_colour_manual(values = c(`TRUE` = pal[["teal_dark"]],
                                 `FALSE` = pal[["amber_dark"]]),
                      guide = "none") +
  scale_y_continuous(limits = c(0.45, 2.75), breaks = NULL) +
  scale_x_continuous(expand = expansion(mult = 0.14)) +
  facet_wrap(~ hyp, scales = "free_x", nrow = 1) +
  labs(x = "estimate", y = NULL,
       title = "A   Transboundary contrasts under the preferred design",
       subtitle = paste0(
         "Teal diamond: transboundary contrast (thick bar: one-sided 95% wild ",
         "cluster bootstrap bound, country,\npre-stated direction, alpha = ",
         "0.05, null-imposed, open toward the predicted side; thin bar: two-",
         "sided\n90% wild cluster bootstrap interval). Amber circle: baseline ",
         "/ placebo, CR2 country-clustered, two-sided.")) +
  theme_paper

# =============================================================================
# PANEL B -- descriptive distance profile and dilution
# =============================================================================
prof <- main$profile
binlv <- levels(prof$bin)
# human-readable bin labels: "[0,10)" -> "0-10", "[200,Inf]" -> "200+"
binlab <- vapply(binlv, function(b) {
  v <- strsplit(gsub("\\[|\\)|\\]", "", b), ",")[[1]]
  if (grepl("Inf", v[2])) paste0(v[1], "+") else paste0(v[1], "–", v[2])
}, "")
n_in <- sum(cfg$profile_bins[-1] <= cfg$influence_zone_km)  # bins within zone
profB <- prof %>% mutate(xb = as.integer(bin))
pB1 <- ggplot(profB, aes(xb, mean_depletion)) +
  annotate("rect", xmin = 0.5, xmax = n_in + 0.5, ymin = -Inf, ymax = Inf,
           fill = "grey88", alpha = 0.55) +
  annotate("text", x = (0.5 + n_in + 0.5) / 2, y = Inf, vjust = 1.4,
           size = 2.0, colour = "grey35", fontface = "italic",
           label = "zone of potential cross-border\nhydraulic influence") +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.6,
                colour = pal[["teal_dark"]]) +
  geom_point(size = 2.1, colour = pal[["teal_dark"]]) +
  scale_x_continuous(breaks = seq_along(binlv), labels = binlab) +
  labs(x = NULL, y = "mean depletion (mm/yr)",
       title = "B   Depletion by distance to border (transboundary segments)",
       subtitle = paste0(
         "Overlap-weighted means with 90% country-cluster bootstrap intervals. ",
         "Descriptive within-TBA profile;\nno matched treated-control test by ",
         "bin (domestic controls lack near-border support).")) +
  theme_paper +
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
pB2 <- ggplot(profB, aes(xb, share_wells)) +
  annotate("rect", xmin = 0.5, xmax = n_in + 0.5, ymin = -Inf, ymax = Inf,
           fill = "grey88", alpha = 0.55) +
  geom_col(fill = pal[["amber"]], width = 0.7) +
  scale_x_continuous(breaks = seq_along(binlv), labels = binlab) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "distance to nearest international border (km)",
       y = "share of wells") +
  theme_paper +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 6))
pB <- pB1 / pB2 + plot_layout(heights = c(3, 1))

# =============================================================================
# PANEL C -- aquifer-level heterogeneity
# =============================================================================
sc <- h3$scatter %>%
  left_join(d %>% select(unit_id, share100 = share_within_100km),
            by = "unit_id") %>%
  left_join(ident$table %>% select(.gid, selected), by = ".gid") %>%
  mutate(selected = coalesce(selected, FALSE),
         share_cls = cut(share100, breaks = cfg$share100_breaks,
                         include.lowest = TRUE,
                         labels = c("< 25%", "25-75%", "> 75%")))
stopifnot(all(sc$share100 >= 0 & sc$share100 <= 1))
qlab <- tibble(
  x = c(Inf, Inf, -Inf, -Inf), y = c(-Inf, Inf, -Inf, Inf),
  hj = c(1.05, 1.05, -0.05, -0.05), vj = c(-0.6, 1.4, -0.6, 1.4),
  lab = c("fast + borderward", "fast, not borderward",
          "borderward, not fast", "neither"))
pC <- ggplot(sc, aes(exc_level, exc_gradient)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_point(aes(size = n_wells, fill = share_cls), shape = 21,
             colour = "grey25", stroke = 0.35, alpha = 0.9) +
  geom_text(data = qlab, aes(x, y, label = lab, hjust = hj, vjust = vj),
            size = 1.9, colour = "grey45", fontface = "italic",
            inherit.aes = FALSE) +
  scale_size_area(max_size = 4.5, name = "wells") +
  scale_fill_manual(values = c("< 25%" = pal[["teal_light"]],
                               "25-75%" = pal[["teal"]],
                               "> 75%" = pal[["teal_dark"]]),
                    name = "wells within 100 km") +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.10)) +
  labs(x = "excess depletion vs matched domestic baseline (mm/yr, shrunken)",
       y = "excess borderward gradient (Fisher z, shrunken)",
       title = "C   Heterogeneity across shared aquifer segments",
       subtitle = paste0(
         "Zero on each axis is the matched domestic baseline. Fill: share of ",
         "wells within 100 km of the border\n(a cue for potential dilution, ",
         "not a determinant of competition). Labelled: segments selected by ",
         "the coded\ncross-analysis rule.")) +
  theme_paper
lab_d <- sc %>% filter(selected)
if (nrow(lab_d) && have_repel)
  pC <- pC + ggrepel::geom_text_repel(
    data = lab_d, aes(label = label_txt), size = 2.0,
    min.segment.length = 0, segment.size = 0.2, box.padding = 0.4,
    max.overlaps = Inf, colour = "grey10")

figmain <- (pA / (pB | pC)) + plot_layout(heights = c(1, 1.6))
save_pdf(figmain, "figure_main_three_panel", 7.2, 8.2)

# =============================================================================
# SI FIGURE 1 -- H1 / H2 robustness forest
# =============================================================================
tab <- rob$tab
fam_levels <- c("Design", "Support diagnostics", "First stage",
                "Spatial uncertainty", "Second stage and parameters")
mk_forest <- function(df, est, lo, hi, pcol, alpha_p, xlab, pref_val, title,
                      wrap_at = 44) {
  df <- df %>%
    filter(ok, is.finite(.data[[est]])) %>%
    mutate(family = factor(family, levels = fam_levels),
           sig = is.finite(.data[[pcol]]) & .data[[pcol]] < alpha_p) %>%
    arrange(family, .data[[est]]) %>%
    mutate(row = factor(stringr::str_wrap(specification, wrap_at),
                        levels = rev(unique(stringr::str_wrap(specification, wrap_at)))))
  wid <- df[[hi]] - df[[lo]]
  typical <- stats::quantile(wid[is.finite(wid)], 0.75, na.rm = TRUE)
  keep <- is.finite(wid) & wid <= 10 * max(typical, .Machine$double.eps)
  rng <- range(c(df[[est]], df[[lo]][keep], df[[hi]][keep], 0, pref_val),
               na.rm = TRUE, finite = TRUE)
  span <- diff(rng); if (!is.finite(span) || span == 0) span <- 1
  lim <- rng + c(-1, 1) * 0.06 * span
  df <- df %>% mutate(off_lo = .data[[lo]] < lim[1], off_hi = .data[[hi]] > lim[2],
                      seg_lo = pmax(.data[[lo]], lim[1]),
                      seg_hi = pmin(.data[[hi]], lim[2]))
  p <- ggplot(df, aes(.data[[est]], row)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey30") +
    geom_vline(xintercept = pref_val, linewidth = 0.35, colour = "grey15",
               linetype = 2) +
    geom_errorbar(aes(xmin = seg_lo, xmax = seg_hi), orientation = "y",
                  width = 0, linewidth = 0.5, colour = "grey55") +
    geom_point(aes(shape = sig), size = 1.8, colour = "grey10", fill = "grey10") +
    scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 1),
                       labels = c(`TRUE` = sprintf("one-sided p < %.2f", alpha_p),
                                  `FALSE` = "not rejected"), name = NULL) +
    facet_grid(family ~ ., scales = "free_y", space = "free_y", switch = "y") +
    coord_cartesian(xlim = lim) +
    labs(x = xlab, y = NULL, title = title,
         subtitle = "Dashed line: preferred estimate. 90% two-sided intervals; intervals beyond the axis are clipped (arrows).") +
    theme_paper +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 0, size = 6.5),
          axis.text.y = element_text(size = 5.8, lineheight = 0.85),
          panel.spacing.y = unit(2, "pt"))
  arr <- grid::arrow(length = unit(0.045, "in"), type = "closed")
  if (any(df$off_lo, na.rm = TRUE))
    p <- p + geom_segment(data = df %>% filter(off_lo),
      aes(x = lim[1] + 0.03 * span, xend = lim[1], y = row, yend = row),
      arrow = arr, linewidth = 0.35, colour = "grey45", inherit.aes = FALSE)
  if (any(df$off_hi, na.rm = TRUE))
    p <- p + geom_segment(data = df %>% filter(off_hi),
      aes(x = lim[2] - 0.03 * span, xend = lim[2], y = row, yend = row),
      arrow = arr, linewidth = 0.35, colour = "grey45", inherit.aes = FALSE)
  p
}
fH1 <- mk_forest(tab %>% filter(!(h1_same %in% TRUE)),
                 "h1_estimate", "h1_ci90_low", "h1_ci90_high",
                 "h1_p_one", cfg$alpha_one,
                 "H1 contrast, mm/yr (90% interval)", rob$tab$h1_estimate[rob$tab$spec_id == "pref"],
                 "H1: full-aquifer depletion contrast")
fH2 <- mk_forest(tab %>% filter(!table_only, h2_scale %in% "z"),
                 "h2_estimate", "h2_ci90_low", "h2_ci90_high",
                 "h2_p_one", cfg$alpha_one,
                 "H2 contrast, Fisher z (90% interval)", rob$tab$h2_estimate[rob$tab$spec_id == "pref"],
                 "H2: borderward gradient contrast")
save_pdf(fH1 | fH2, "figure_si_robustness_h1_h2", 12, 8.5)

# =============================================================================
# SI FIGURE 2 -- localization robustness
# =============================================================================
profs <- bind_rows(rob$locC$profiles,
                   rob$locC$prof_unw, rob$locC$prof_ctr) %>%
  mutate(mid = map_dbl(as.character(bin), function(b) {
    v <- suppressWarnings(as.numeric(strsplit(gsub("\\[|\\)|\\]", "", b), ",")[[1]]))
    if (!is.finite(v[2])) v[1] * 1.5 else mean(v)
  }))
pL1 <- ggplot(profs, aes(mid, mean_depletion, group = bins, linetype = bins)) +
  annotate("rect", xmin = 0, xmax = cfg$influence_zone_km, ymin = -Inf,
           ymax = Inf, fill = "grey90", alpha = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_line(linewidth = 0.4, colour = "grey30") +
  geom_point(size = 1.4, colour = "grey15") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.35,
                colour = "grey55") +
  scale_x_continuous(trans = "sqrt",
                     breaks = c(10, 50, 100, 200, 400)) +
  labs(x = "distance to border, km (sqrt scale; bin midpoints)",
       y = "mean depletion (mm/yr)",
       title = "Localization profile under alternative constructions",
       subtitle = "Bin definitions, unweighted, and within-segment centred variants. Shaded: 0-100 km zone. 90% country-bootstrap intervals.") +
  theme_paper
loo <- rob$locC$loo %>%
  arrange(contrast_within_minus_beyond) %>%
  mutate(row = factor(specification, levels = specification))
pL2 <- ggplot(loo, aes(contrast_within_minus_beyond, row)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_vline(xintercept = rob$locC$table$contrast_within_minus_beyond[1],
             linetype = 2, linewidth = 0.35, colour = "grey15") +
  geom_errorbar(aes(xmin = contrast_lo, xmax = contrast_hi), orientation = "y",
                width = 0, linewidth = 0.45, colour = "grey55") +
  geom_point(size = 1.7, colour = "grey10") +
  labs(x = "within-100km minus beyond-100km mean depletion (mm/yr)",
       y = NULL, title = "Leave-one-country-out (descriptive contrast)",
       subtitle = "Dashed line: full-sample value. 90% country-bootstrap intervals.") +
  theme_paper
save_pdf(pL1 / pL2 + plot_layout(heights = c(1, 1.4)),
         "figure_si_localization_robustness", 7, 8.5)

# =============================================================================
# SI FIGURE 3 -- influence (leave-one-out)
# =============================================================================
mk_inf <- function(df, title, n_max = 40) {
  df <- df %>% filter(is.finite(h2_estimate)) %>%
    arrange(desc(abs(d_h2))) %>% slice_head(n = n_max) %>%
    arrange(h2_estimate) %>%
    mutate(row = factor(dropped, levels = dropped))
  ggplot(df, aes(h2_estimate, row)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey30") +
    geom_vline(xintercept = rob$tab$h2_estimate[rob$tab$spec_id == "pref"],
               linetype = 2, linewidth = 0.35, colour = "grey15") +
    geom_errorbar(aes(xmin = h2_ci90_low, xmax = h2_ci90_high),
                  orientation = "y", width = 0, linewidth = 0.4,
                  colour = "grey55") +
    geom_point(aes(shape = sign_reversal_h2), size = 1.7, colour = "grey10") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                       labels = c(`FALSE` = "same sign", `TRUE` = "sign reversal"),
                       name = NULL) +
    labs(x = "H2 contrast, Fisher z (90% interval)", y = NULL, title = title,
         subtitle = "Dashed line: preferred estimate.") +
    theme_paper +
    theme(axis.text.y = element_text(size = 5.6))
}
pI1 <- mk_inf(rob$inf_cc, "Leave one country out")
pI2 <- mk_inf(rob$inf_seg, "Leave one treated segment out")
pI3 <- mk_inf(rob$inf_aq, "Leave one physical aquifer out")
save_pdf(pI1 | pI2 | pI3, "figure_si_influence", 12, 7)

sayf("36_figures.R done in %.1f min.\n",
     as.numeric(difftime(Sys.time(), t0, units = "mins")))


# =============================================================================
# WELL-DENSITY MAP, OUTLIER INSETS AND THE CONNECTIVITY FIGURE
#
# Folded in from the reviewer-response scripts of the revision round.
#   figure_main_wells_map.pdf          Fig. 2D  global observation-well density
#   figure_main_wells_insets.pdf       Fig. 2E  the eight labelled segments
#   figure_si_wells_connectivity.pdf   Fig. S2  contrasts by hydraulic
#                                               connectivity (rob$conn_tab)
#
# The two map figures need the shapefiles under map/ (Jasechko aquifer
# polygons and geoBoundaries ADM0), which are third-party inputs and are not
# redistributed here -- see the data table in the top-level README. The block
# stops with a clear message if they are absent.
# =============================================================================
suppressPackageStartupMessages({library(sf); library(maps)})
sf::sf_use_s2(FALSE)
wells <- fs_wells(cfg)

# theme_conn preserves the exact appearance of the published Fig. S2: it is
# theme_paper without the axis lines/ticks and bold strips added for the main
# figure, and is used only by the connectivity panel below.
theme_conn <- theme_minimal(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 8.5),
        plot.subtitle = element_text(size = 6.5, colour = "grey30"),
        legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 7))

# =============================================================================
# 1. TASK I -- PANEL A WORLD MAP AND OUTLIER-AQUIFER INSETS
# =============================================================================
say("\n---- maps ----\n")
# Aquifer geometries: Jasechko et al. (2024) polygons split by country
# (map/jasechko_aquifs/jasechko_CountrySplit.shp; key Aquifer x shapeGroup
# matches every analysis segment). Political background: geoBoundaries CGAZ
# ADM0; land borders: map/geoBoundariesCGAZ_ADM0/LB.shp.
suppressPackageStartupMessages(library(sf))
sf::sf_use_s2(FALSE)
map_dir <- file.path(cfg$root, "map")
shp_aq  <- file.path(map_dir, "jasechko_aquifs", "jasechko_CountrySplit.shp")
shp_adm <- file.path(map_dir, "geoBoundariesCGAZ_ADM0", "geoBoundariesCGAZ_ADM0.shp")
shp_lb  <- file.path(map_dir, "geoBoundariesCGAZ_ADM0", "LB.shp")
for (f in c(shp_aq, shp_adm, shp_lb))
  if (!file.exists(f)) stop("Shapefile not found: ", f,
                            " (copy the Final/map folder next to the scripts).")
# ---- geometry: simplify once, cache, reuse -----------------------------------
# geoBoundariesCGAZ_ADM0.shp is 155 MB and LB.shp 51 MB. At the printed size of
# this figure (9.5 in for 350 deg of longitude, ~37 deg/in) one device pixel at
# 600 dpi spans ~0.06 deg, so vertices finer than that cannot be resolved: they
# only inflate the PDF and the run time. Two derived layers are built --
#   world: heavy simplification + small islands dropped, for the global map
#   inset: light simplification, cropped to the union of the inset boxes
# -- and both are cached, so the shapefiles are read and cleaned only when they
# or these settings change. Nothing here touches 31_config.R or 32_core.R, so
# the analysis caches are unaffected.
GEOM <- list(
  tol_world      = 0.08,   # degrees (~9 km) -- below one device pixel on the world map
  tol_inset      = 0.005,  # degrees (~550 m) -- insets span only a few degrees
  min_poly_km2   = 200,    # drop islands smaller than this from the world basemap
  inset_pad_frac = 0.22,   # must match inset_of() below
  inset_pad_min  = 0.30)

seg_info <- d %>% select(unit_id, Aquifer, CC, n, TBn)
sel <- ident$table %>% filter(selected) %>%
  arrange(rank_raw) %>%
  mutate(tag = paste0("(", letters[seq_len(n())], ")"))

# Stamp on the shapefiles themselves and on the settings above -- deliberately
# NOT on cfg$source_files, because the geometry does not depend on the analysis
# code and should survive edits to it.
geom_stamp <- function(files, params) list(
  cfg_md5   = .md5_of("geometry"),
  ov_md5    = .md5_of(params),
  code      = .file_stamp(files),
  data      = character(0),
  extra_md5 = .md5_of(NULL),
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  written_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

# Drop parts below an area threshold; keeps mainlands, discards islet noise.
drop_small_parts <- function(x, min_km2) {
  if (!is.finite(min_km2) || min_km2 <= 0) return(x)
  x <- x[!st_is_empty(st_geometry(x)), , drop = FALSE]
  parts <- suppressWarnings(st_cast(x, "POLYGON", warn = FALSE))
  keep <- as.numeric(st_area(parts)) / 1e6 >= min_km2
  keep[is.na(keep)] <- FALSE
  parts[keep, , drop = FALSE]
}
simplify_q <- function(x, tol) {
  y <- suppressWarnings(st_simplify(x, dTolerance = tol,
                                    preserveTopology = TRUE))
  y[!st_is_empty(st_geometry(y)), , drop = FALSE]
}
n_vert <- function(x) sum(vapply(st_geometry(x), function(g)
  nrow(sf::st_coordinates(g)), numeric(1)))

geom_path <- file.path(cfg$cache_dir, "basemap_geometry.rds")
geom_st <- geom_stamp(c(shp_aq, shp_adm, shp_lb),
                      list(GEOM = GEOM, units = sort(d$unit_id),
                           sel = sort(sel$.gid)))
GEO <- cache_read(geom_path, geom_st, cfg, force = FORCE)
if (is.null(GEO)) {
  say("  building simplified basemap (one-off; cached afterwards)\n")
  aq_full <- st_read(shp_aq, quiet = TRUE) %>%
    mutate(unit_id = paste0(Aquifer, "_", shapeGroup)) %>%
    filter(unit_id %in% d$unit_id) %>%
    st_make_valid()
  stopifnot(all(d$unit_id %in% aq_full$unit_id))
  # Area from the FULL polygons: well density must not depend on simplification.
  aq_full$area_km2_poly <- as.numeric(st_area(aq_full)) / 1e6
  aq_full <- aq_full %>%
    left_join(seg_info %>% select(unit_id, n, TBn), by = "unit_id") %>%
    mutate(dens = n / pmax(area_km2_poly, 1))

  adm0 <- st_read(shp_adm, quiet = TRUE) %>% st_make_valid()
  lb0  <- st_read(shp_lb,  quiet = TRUE)
  v0 <- c(adm = n_vert(adm0), lb = n_vert(lb0), aq = n_vert(aq_full))

  # World layers: simplify, then drop islets.
  adm_w <- drop_small_parts(simplify_q(adm0, GEOM$tol_world), GEOM$min_poly_km2)
  lb_w  <- simplify_q(lb0, GEOM$tol_world)
  aq_w  <- simplify_q(aq_full, GEOM$tol_world)
  # Simplification must not silently delete a segment from the world map.
  stopifnot(setequal(aq_w$unit_id, aq_full$unit_id))

  # Inset layers: crop ONCE to the union of the inset boxes, then simplify
  # lightly. Previously each inset clipped the full 155 MB layer -- eight
  # insets x two layers = sixteen intersections over global geometry.
  box_of <- function(gid) {
    bb <- st_bbox(aq_full %>% filter(unit_id == gid))
    pad <- pmax(GEOM$inset_pad_min,
                GEOM$inset_pad_frac * max(bb["xmax"] - bb["xmin"],
                                          bb["ymax"] - bb["ymin"]))
    st_bbox(c(xmin = unname(bb["xmin"]) - pad, xmax = unname(bb["xmax"]) + pad,
              ymin = unname(bb["ymin"]) - pad, ymax = unname(bb["ymax"]) + pad),
            crs = st_crs(adm0))
  }
  boxes <- lapply(sel$.gid, box_of)
  names(boxes) <- sel$.gid
  hull <- st_as_sfc(st_bbox(c(
    xmin = min(vapply(boxes, function(b) unname(b["xmin"]), 0)),
    xmax = max(vapply(boxes, function(b) unname(b["xmax"]), 0)),
    ymin = min(vapply(boxes, function(b) unname(b["ymin"]), 0)),
    ymax = max(vapply(boxes, function(b) unname(b["ymax"]), 0))),
    crs = st_crs(adm0)))
  # s2 is off (set above), so these lon/lat operations are planar by design --
  # correct for axis-aligned bbox clipping. Silence the standard notice.
  crop_q <- function(x, g) suppressMessages(suppressWarnings(st_crop(x, g)))
  adm_i <- simplify_q(crop_q(adm0, hull), GEOM$tol_inset)
  lb_i  <- simplify_q(crop_q(lb0,  hull), GEOM$tol_inset)

  v1 <- c(adm = n_vert(adm_w), lb = n_vert(lb_w), aq = n_vert(aq_w))
  sayf("  vertices: adm %s -> %s, lb %s -> %s, aquifers %s -> %s (%.0f%% of original)\n",
       format(v0["adm"], big.mark = ","), format(v1["adm"], big.mark = ","),
       format(v0["lb"],  big.mark = ","), format(v1["lb"],  big.mark = ","),
       format(v0["aq"],  big.mark = ","), format(v1["aq"],  big.mark = ","),
       100 * sum(v1) / sum(v0))
  GEO <- list(aq_full = aq_full, aq_w = aq_w, adm_w = adm_w, lb_w = lb_w,
              adm_i = adm_i, lb_i = lb_i, boxes = boxes,
              vertices_before = v0, vertices_after = v1, settings = GEOM)
  cache_write(GEO, geom_path, geom_st)
} else {
  sayf("  [cached] simplified basemap (%.0f%% of original vertices)\n",
       100 * sum(GEO$vertices_after) / sum(GEO$vertices_before))
}
aq_sf <- GEO$aq_full          # full detail: the subject of each inset
aq_s  <- GEO$aq_w             # simplified: world map
adm_s <- GEO$adm_w
lb_s  <- GEO$lb_w
lab_pts <- suppressWarnings(st_point_on_surface(
  aq_sf %>% filter(unit_id %in% sel$.gid))) %>%
  left_join(sel %>% select(.gid, tag), by = c(unit_id = ".gid"))
lab_xy <- bind_cols(st_drop_geometry(lab_pts)["tag"],
                    as_tibble(st_coordinates(lab_pts)))
pmap <- ggplot() +
  geom_sf(data = adm_s, fill = "grey96", colour = "grey75", linewidth = 0.08) +
  geom_sf(data = lb_s, colour = "grey55", linewidth = 0.14) +
  geom_sf(data = aq_s %>% filter(TBn == 0), aes(fill = log10(dens)),
          colour = "grey35", linewidth = 0.06, alpha = 0.95) +
  geom_sf(data = aq_s %>% filter(TBn == 1), aes(fill = log10(dens)),
          colour = "black", linewidth = 0.45, alpha = 0.98) +
  ggrepel::geom_label_repel(
    data = lab_xy, aes(X, Y, label = tag), size = 2.4,
    label.padding = 0.12, box.padding = 0.28, min.segment.length = 0,
    segment.size = 0.25, max.overlaps = Inf, seed = 1) +
  scale_fill_viridis_c(name = expression(log[10]~"wells km"^-2),
                       option = "magma", direction = -1) +
  coord_sf(xlim = c(-170, 180), ylim = c(-58, 78), expand = FALSE) +
  labs(title = "Monitored aquifer-country segments: observation-well density",
       subtitle = paste0("Aquifer polygons from Jasechko et al. (2024), split ",
                         "at international borders; grey lines: land borders ",
                         "(geoBoundaries CGAZ). Thick black outlines: ",
                         "transboundary segments. Labels (a)-(",
                         letters[nrow(sel)],
                         "): segments selected by the coded cross-analysis rule.")) +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.3, colour = "grey30"),
        legend.position = "bottom",
        legend.key.height = unit(0.25, "cm"), legend.key.width = unit(1, "cm"))
save_pdf(pmap, "figure_main_wells_map", 9.5, 5.6)

# ---- per-aquifer insets ------------------------------------------------------
# Wells coloured by the sign of their fitted 2000-2022 trend. The well file
# carries the Theil-Sen point estimate only (no per-well standard error), so
# "no trend" is defined as |GWSlp| <= 10 mm/yr; this band is stated in the
# caption rather than presented as a significance test.
trend_band <- 10
inset_of <- function(gid, tag) {
  poly <- aq_sf %>% filter(unit_id == gid)   # full detail: the inset's subject
  w <- wells %>% filter(unit_id == gid)
  info <- seg_info %>% filter(unit_id == gid)
  # Box computed once when the geometry cache was built, so the inset extent
  # here and the crop used to build GEO$adm_i / GEO$lb_i cannot drift apart.
  bx <- GEO$boxes[[gid]]
  xl <- unname(c(bx["xmin"], bx["xmax"]))
  yl <- unname(c(bx["ymin"], bx["ymax"]))
  box <- st_as_sfc(bx)
  # Cropping the pre-cropped, lightly simplified layers: cheap, and planar by
  # design (s2 is off), so the routine lon/lat notice is silenced.
  adm_c <- suppressMessages(suppressWarnings(st_crop(GEO$adm_i, box)))
  lb_c  <- suppressMessages(suppressWarnings(st_crop(GEO$lb_i,  box)))
  w <- w %>% mutate(cls = case_when(GWSlp >  trend_band ~ "depleting",
                                    GWSlp < -trend_band ~ "replenishing",
                                    TRUE ~ "no trend"))
  ggplot() +
    geom_sf(data = adm_c, fill = "grey94", colour = "grey65",
            linewidth = 0.15) +
    geom_sf(data = lb_c, colour = "grey20", linewidth = 0.45) +
    geom_sf(data = poly, fill = "grey85", alpha = 0.35, colour = "black",
            linewidth = 0.55) +
    geom_point(data = w, aes(lon, lat, fill = cls), shape = 21, size = 1.1,
               stroke = 0.15, colour = "grey30") +
    scale_fill_manual(values = c(depleting = "#B2182B", `no trend` = "white",
                                 replenishing = "#2166AC"), name = NULL,
                      drop = FALSE) +
    coord_sf(xlim = xl, ylim = yl, expand = FALSE) +
    labs(title = stringr::str_wrap(paste0(tag, " ", info$Aquifer, " (",
                                          info$CC, "), ", info$n, " wells"),
                                   34)) +
    theme_void(base_size = 7) +
    theme(plot.title = element_text(size = 5.6, face = "bold",
                                    lineheight = 0.9),
          panel.border = element_rect(colour = "grey40", fill = NA,
                                      linewidth = 0.3),
          legend.position = "none")
}
insets <- map2(sel$.gid, sel$tag, inset_of)
legend_stub <- ggplot(tibble(x = 1:3, cls = c("depleting", "no trend",
                                              "replenishing")),
                      aes(x, 1, fill = cls)) +
  geom_point(shape = 21, size = 2.4, colour = "grey30") +
  scale_fill_manual(values = c(depleting = "#B2182B", `no trend` = "white",
                               replenishing = "#2166AC"), name = NULL) +
  theme_void() + theme(legend.position = "bottom",
                       legend.text = element_text(size = 7))
leg <- function(p) {
  g <- ggplotGrob(p)
  i <- grep("guide-box", sapply(g$grobs, function(x) x$name))[1]
  if (is.na(i)) return(grid::nullGrob())
  g$grobs[[i]]
}
pins <- wrap_plots(insets, ncol = 4) /
  patchwork::wrap_elements(leg(legend_stub)) +
  plot_layout(heights = c(8, 1))
save_pdf(pins, "figure_main_wells_insets", 10, 6.2)


# ---- SI Fig. S2: contrasts by hydraulic connectivity ------------------------

conn_tab <- rob$conn_tab
if (is.null(conn_tab)) {
  say("  connectivity subsets absent (IGRAC properties not found) -- skipping Fig. S2\n")
} else {

# subset order: ROI bins first (short/intermediate/long), then lithology
ord <- c(grep("ROI", names(subruns), value = TRUE)[
           order(match(grep("ROI", names(subruns), value = TRUE),
                       c("short ROI (highest connectivity attenuation)",
                         "intermediate ROI", "long ROI")))],
         sort(grep("ROI", names(subruns), value = TRUE, invert = TRUE)))
ct <- conn_tab %>%
  mutate(subset = factor(subset, levels = rev(ord)),
         grp = if_else(grepl("ROI", subset), "by radius of influence (Tr, S)",
                       "by lithology"))
mk_co <- function(est, lo, hi, pone, pref_v, xlab, title) {
  ctp <- ct %>% mutate(sig = is.finite(.data[[pone]]) & .data[[pone]] < cfg$alpha_one)
  ggplot(ctp, aes(.data[[est]], subset)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
    geom_vline(xintercept = pref_v, linetype = 2, linewidth = 0.35,
               colour = "grey15") +
    geom_errorbar(aes(xmin = .data[[lo]], xmax = .data[[hi]]),
                  orientation = "y", width = 0, linewidth = 0.55,
                  colour = pal[["teal"]], alpha = 0.7) +
    # point (diamond) size flags one-sided significance at alpha_one (95%),
    # the same test p1 reports
    geom_point(aes(size = sig), shape = 18, colour = pal[["teal_dark"]]) +
    geom_text(aes(label = sprintf("p1 = %s, k = %d", pstar(.data[[pone]]),
                                  n_treated)),
              vjust = -1.2, size = 1.9, colour = "grey30") +
    scale_size_manual(values = c(`TRUE` = 4.2, `FALSE` = 2.2),
                      labels = c(`TRUE` = sprintf("one-sided p < %.2f", cfg$alpha_one),
                                 `FALSE` = "not rejected"),
                      breaks = c("FALSE", "TRUE"), drop = FALSE,
                      limits = c("FALSE", "TRUE"), name = NULL) +
    facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = xlab, y = NULL, title = title) +
    theme_conn +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 0, size = 6.5),
          axis.text.y = element_text(size = 6.5),
          plot.margin = margin(4, 10, 4, 4))
}
pc1 <- mk_co("h1_estimate", "h1_lo", "h1_hi", "h1_p_one", pref$h1$estimate,
             "H1 contrast, mm/yr\n(90% wild cluster bootstrap CI)", "H1 by TBA subset")
pc2 <- mk_co("h2_estimate", "h2_lo", "h2_hi", "h2_p_one", pref$h2$estimate,
             "H2 contrast, Fisher z\n(90% wild cluster bootstrap CI)", "H2 by TBA subset")
pconn <- (pc1 | pc2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
pconn <- pconn +
  plot_annotation(
    title = "Transboundary contrasts by hydraulic-property subsets of treated segments",
    subtitle = paste(strwrap(paste0(
      "Each subset keeps the named treated segments and the full domestic ",
      "control pool; the complete pipeline (matching, weights, second stage, ",
      "inference) is re-estimated. Dashed line: all-TBA preferred estimate. ",
      "Error bars: 90% wild cluster bootstrap interval (test inversion). ",
      "p1 = one-sided wild cluster bootstrap p (pre-stated direction), the ",
      "same route the interval is built from; k = treated segments. Larger ",
      "diamond: one-sided p1 < 0.05."), width = 140), collapse = "\n"),
    theme = theme(plot.title = element_text(face = "bold", size = 9),
                  plot.subtitle = element_text(size = 6.3, colour = "grey30")))
save_pdf(pconn, "figure_si_wells_connectivity", 10.5, 5.6)
}
