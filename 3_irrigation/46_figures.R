# =============================================================================
# 46_figures.R -- CANONICAL FIGURES (PDF ONLY; figure/)
#
#   figure/figure_main_irrigation.pdf   Fig. 3: A intensity, B Gini, C Lorenz
#   figure/figure_si_robustness.pdf     SI: Fig. S8 counterparts + GWf>0
#   figure/figure_si_river_exclusion.pdf SI: river-border exclusion block
#   figure/figure_si_balance.pdf        SI: ensemble balance diagnostics (R3.2)
#
# Design (revised 2026-08; supersedes the v3 mirrored-density encoding):
#   - ONE quantity per panel. The v3 panel drew two mirrored half-densities on
#     a shared axis -- the matched-control baseline (a LEVEL, e.g. 2.1% of
#     area) and the transboundary contrast (a DIFFERENCE, e.g. +1.3 pp). Equal
#     visual weight invited a two-sample reading, when the claim being made is
#     that each difference distribution sits away from ZERO; and zero was drawn
#     across both, though it is meaningful only for the difference. The panel
#     now plots the difference alone, as a full violin anchored on a solid zero
#     line, so the axis carries one measure and the intended comparison is the
#     only one the geometry supports.
#   - The baseline is NOT dropped, it is demoted: mean and 5-95% range appear
#     in the panel subtitle, where they give the difference its scale without
#     competing with it for the reader's eye.
#   - Plotting two LEVELS instead (control vs transboundary) would make the
#     reader's instinct correct but is statistically worse: the difference is
#     computed WITHIN each control realization, so it is paired and far more
#     precise than two marginal level distributions would imply.
#   - mean point, IQR (thick) and 5-95% (thin) intervals on the difference;
#   - one labelled significance bar per panel (share of realizations with
#     p < 0.1: blue positive / grey n.s. / vermilion negative), with legend;
#   - intensity outcomes in percentage points;
#   - descriptive, non-causal panel titles (reviewer R1.3);
#   - panel C: Lorenz curves for three named aquifers, Gini computed from the
#     current curves; NO map insets (reviewer R2.3).
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!exists("CFG")) source(file.path(.root, "41_config.R"))
if (!exists("run_spec")) source(file.path(CFG$root, "42_core.R"))
if (is.null(.log$path)) log_open(CFG$log_file)
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
cfg <- CFG

main <- readRDS(file.path(cfg$cache_dir, "main_objects.rds"))
rob  <- readRDS(file.path(cfg$cache_dir, "robustness_objects.rds"))
dir.create(cfg$fig_dir, showWarnings = FALSE, recursive = TRUE)

# PDF device. cairo_pdf embeds fonts (which journals require) but needs a
# working cairo build -- on macOS, XQuartz. Without it grDevices::cairo_pdf
# warns "failed to load cairo DLL" and writes NOTHING, so the figures on disk
# silently stay at whatever the last successful run produced, while
# 48_run_all.R's assertion 1 -- which only tests file.exists() -- still passes.
# capabilities("cairo") is NOT a reliable check: it can report TRUE on a
# machine where the library fails to load. So probe by writing a real file,
# fall back to base pdf(), and verify every output after writing it.
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
}
save_pdf <- function(p, name, w, h) {
  path <- file.path(cfg$fig_dir, paste0(name, ".pdf"))
  if (file.exists(path)) unlink(path)  # a failed write must not leave a stale file
  if (.cairo_ok) {
    ggsave(path, p, width = w, height = h, device = grDevices::cairo_pdf)
  } else {
    ggsave(path, p, width = w, height = h, device = grDevices::pdf,
           encoding = "WinAnsi.enc")
  }
  if (!file.exists(path) || file.size(path) < 1000)
    stop("save_pdf: ", basename(path), " was not written (device failure). ",
         "Check the PDF device message above.")
  say("  wrote ", basename(path), "\n")
}

# ---- shared aesthetics ------------------------------------------------------
COL_EFF <- "#E69F00"                                  # Okabe-Ito (CVD-safe)
COL_SIG  <- c(positive = "#0072B2", `n.s.` = "#BDBDBD", negative = "#D55E00")
PP_OUTCOMES <- c("Ir", "OverIR3", "OverIR6", "IrNeed3", "GW")

theme_paper <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line  = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
        plot.subtitle = element_text(size = 6, colour = "grey30", hjust = 0.5),
        plot.tag   = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        legend.key.size = unit(0.32, "cm"),
        legend.title = element_text(size = 7),
        legend.text  = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

sum_stats <- function(v) data.frame(
  m = mean(v), q05 = quantile(v, .05), q25 = quantile(v, .25),
  md = median(v), q75 = quantile(v, .75), q95 = quantile(v, .95))

# Full (two-sided) violin for the difference. Built as two polygons -- the mass
# above zero and the mass below -- rather than one, so that the boundary point
# at y = 0 is interpolated exactly onto the zero line and the outline meets it
# cleanly. Both parts take the same fill here; keeping them separate leaves the
# option of colouring them by sign without touching the geometry.
violin_parts <- function(v, x0 = 1, w = 0.34) {
  v <- v[is.finite(v)]
  dd <- stats::density(v, n = 512)
  d  <- dd$y / max(dd$y)
  poly_for <- function(keep, lab) {
    if (!any(keep)) return(NULL)
    yy <- dd$x[keep]; ddv <- d[keep]
    if (min(dd$x) < 0 && max(dd$x) > 0) {          # add the exact zero crossing
      d0 <- stats::approx(dd$x, d, xout = 0)$y
      yy <- c(yy, 0); ddv <- c(ddv, d0)
      o <- order(yy); yy <- yy[o]; ddv <- ddv[o]
    }
    data.frame(x = c(x0 - w * ddv, rev(x0 + w * ddv)),
               y = c(yy, rev(yy)), part = lab)
  }
  rbind(poly_for(dd$x >= 0, "positive"), poly_for(dd$x <= 0, "negative"))
}

# The matched-control baseline, demoted to the subtitle: it sets the scale the
# difference should be read against, without sharing an axis with it.
fmt_base <- function(b, sc) {
  u <- if (sc == 100) "%" else ""
  sprintf("matched-control baseline %.2f%s  (5-95%%: %.2f-%.2f)",
          mean(b), u, quantile(b, .05), quantile(b, .95))
}

# One specification panel: the transboundary difference alone, on zero.
spec_panel <- function(res, ylab = NULL, sig_level = cfg$sig_level) {
  df <- res$summary_df
  sc <- if (res$spec$outcome %in% PP_OUTCOMES) 100 else 1
  b  <- df$int_eff * sc; e <- df$treat_eff * sc

  ply <- violin_parts(e)
  pts <- cbind(sum_stats(e), x = 1)

  sig <- ifelse(df$treat_p >= sig_level, "n.s.",
                ifelse(df$treat_eff > 0, "positive", "negative"))
  tb  <- as.data.frame(table(factor(sig, c("positive", "n.s.", "negative"))))
  tb$p <- tb$Freq / sum(tb$Freq)
  # The axis carries ONE measure, so its range is set by the difference alone,
  # and zero is forced into view -- it is the comparison being made.
  ymin <- min(c(e, 0)); ymax <- max(c(e, 0))
  yspan <- ymax - ymin; if (yspan == 0) yspan <- 1
  y_bar <- ymin - 0.15 * yspan
  x0 <- 1 - 0.42; ww <- 0.84; xc <- x0; bars <- NULL; labs <- NULL
  for (i in seq_len(nrow(tb))) {
    p <- tb$p[i]; if (p == 0) next
    bars <- rbind(bars, data.frame(x = xc, xend = xc + ww * p, y = y_bar,
                                   sig = as.character(tb$Var1[i])))
    if (p >= 0.12 && tb$Var1[i] != "n.s.")
      labs <- rbind(labs, data.frame(x = xc + ww * p / 2, y = y_bar,
                                     lab = sprintf("%d%%", round(100 * p))))
    xc <- xc + ww * p
  }

  g <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.45, colour = "grey25") +
    geom_polygon(data = ply, aes(x, y, group = part), fill = COL_EFF,
                 alpha = 0.9, colour = NA) +
    geom_linerange(data = pts, aes(x = x, ymin = q05, ymax = q95),
                   linewidth = 0.35, colour = "grey35") +
    geom_linerange(data = pts, aes(x = x, ymin = q25, ymax = q75),
                   linewidth = 1.05, colour = "grey15") +
    geom_point(data = pts, aes(x, m), size = 1.4, colour = "grey5") +
    geom_segment(data = bars, aes(x = x, xend = xend, y = y, yend = y,
                                  colour = sig), linewidth = 2.6, lineend = "butt") +
    scale_colour_manual(sprintf("Share of %d realizations (p < %.1f)",
                                nrow(df), sig_level),
                        values = COL_SIG, drop = FALSE) +
    scale_x_continuous(limits = c(0.5, 1.5), breaks = NULL) +
    coord_cartesian(ylim = c(ymin - 0.22 * yspan, ymax + 0.05 * yspan)) +
    labs(x = NULL, y = ylab, title = res$spec$title,
         subtitle = fmt_base(b, sc)) +
    theme_paper + theme(legend.position = "none")
  if (!is.null(labs))
    g <- g + geom_text(data = labs, aes(x, y, label = lab), size = 2.3,
                       colour = "white", fontface = "bold")
  g
}

strip_y <- function(g) g + labs(y = NULL)

# One shared legend row for the difference panels. The baseline no longer has a
# key: it is not drawn, it is reported in each panel's subtitle.
make_legend_row <- function(n = cfg$n_realizations, sig = cfg$sig_level) {
  # Two rows PER group, not one. This plot is never drawn -- only its legend is
  # extracted below -- but with a single point per colour group geom_line has no
  # segment to draw and ggplot warns "Each group consists of only one
  # observation" once per figure that uses this legend. The key glyphs are
  # synthesised by the guide system from the scale, so duplicating the rows
  # changes nothing about the legend and silences the warning at source.
  lab_eff <- "Transboundary difference across realizations"
  dd <- data.frame(
    x = rep(1:2, times = 2),
    y = rep(1:2, times = 2),
    f = lab_eff,
    s = factor(rep(c("positive", "negative"), each = 2),
               levels = c("positive", "n.s.", "negative")))
  gl <- ggplot(dd) +
    geom_col(aes(x, y, fill = f)) +
    geom_line(aes(x, y, colour = s), linewidth = 2) +
    scale_fill_manual(NULL, values = setNames(COL_EFF, lab_eff)) +
    scale_colour_manual(sprintf("Share of %d realizations (p < %.1f):", n, sig),
                        values = COL_SIG, drop = FALSE) +
    guides(fill = guide_legend(order = 1),
           colour = guide_legend(order = 2, nrow = 1)) +
    theme_paper + theme(legend.position = "bottom",
                        legend.box = "horizontal",
                        legend.spacing.x = unit(0.5, "cm"))
  patchwork::wrap_elements(full = cowplot::get_legend(gl))
}

# =============================================================================
# FIG 3 -- MAIN FIGURE
# =============================================================================
say("figure_main_irrigation:\n")
pA1 <- spec_panel(main$IrIntens,
                  "Difference in irrigation intensity (pp of area)") +
  labs(tag = "A")
pA2 <- strip_y(spec_panel(main$Overdraft))
pA3 <- strip_y(spec_panel(main$IrRivsInt))
pB1 <- spec_panel(main$IrGini, "Difference in borderward Gini") + labs(tag = "B")
pB2 <- strip_y(spec_panel(main$IrRivsGini))
pB3 <- strip_y(spec_panel(main$GWGini))
pB4 <- strip_y(spec_panel(main$GWRivsGini))

# ---- panel C: Lorenz curves, no maps (reviewer R2.3) ------------------------
lx <- lorenz_examples(cfg)
short_name <- function(nm) {
  nm <- sub(" Aquifer System.*", "", nm)
  nm <- sub(" Basin alluvial aquifer system", "", nm)
  sub("Saq-Ram", "Saq-Ram / Disi", nm)
}
lz <- do.call(rbind, lapply(names(lx), function(nm) {
  cc <- lx[[nm]]$curve
  data.frame(x = cc$LandPct, y = cc$IrrPct,
             lab = sprintf("%s (%s)\nG = %+.2f", short_name(lx[[nm]]$name),
                           lx[[nm]]$country, lx[[nm]]$gini))
}))
lz_cols <- setNames(c("#D55E00", "#E69F00", "#0072B2")[seq_along(unique(lz$lab))],
                    unique(lz$lab))
lab_pos <- data.frame(x = c(26, 63, 58), y = c(86, 34, 5),
                      lab = unique(lz$lab))
pC <- ggplot(lz, aes(x, y, colour = lab)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3,
              colour = "grey55") +
  geom_line(linewidth = 0.8) +
  geom_text(data = lab_pos, aes(label = lab), hjust = 0, size = 2.0,
            lineheight = 0.95, show.legend = FALSE) +
  scale_colour_manual(values = lz_cols, guide = "none") +
  coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(x = "Cumulative land area by distance to border (%)",
       y = "Cumulative irrigated area (%)",
       title = "Borderward Gini construction", tag = "C") +
  theme_paper + theme(axis.text.x = element_text(size = 6.5),
                      axis.text.y = element_text(size = 6.5),
                      axis.ticks.x = element_line(linewidth = 0.3))

fig3 <- ((pA1 | pA2 | pA3 | pC) + patchwork::plot_layout(widths = c(1, 1, 1, 1.9))) /
        (pB1 | pB2 | pB3 | pB4) /
        make_legend_row() +
  patchwork::plot_layout(heights = c(1, 1, 0.14))
save_pdf(fig3, "figure_main_irrigation", w = 9.0, h = 6.6)

# =============================================================================
# SI -- ROBUSTNESS COUNTERPARTS (Fig. S8 + attribution sensitivity)
# =============================================================================
say("figure_si_robustness:\n")
s1 <- spec_panel(rob$CropSuitInt, "Transboundary difference") + labs(tag = "A")
s2 <- strip_y(spec_panel(rob$IrNeedInt)) + labs(tag = "B")
s3 <- strip_y(spec_panel(rob$Overdraft_Irrig)) + labs(tag = "C")
s4 <- strip_y(spec_panel(rob$Overdraft6)) + labs(tag = "D")
s5 <- spec_panel(rob$IrIntDoubleRobust, "Transboundary difference") + labs(tag = "E")
s6 <- strip_y(spec_panel(rob$IrGiniDoubleRobust)) + labs(tag = "F")
s7 <- strip_y(spec_panel(rob$GWShareGt0)) + labs(tag = "G")
fig_si <- (s1 | s2 | s3 | s4) / (s5 | s6 | s7 | patchwork::plot_spacer()) /
  make_legend_row() + patchwork::plot_layout(heights = c(1, 1, 0.14))
save_pdf(fig_si, "figure_si_robustness", w = 9.0, h = 6.8)

# =============================================================================
# SI -- RIVER-BORDER EXCLUSION (reviewer R5.4)
# =============================================================================
say("figure_si_river_exclusion:\n")
e1 <- spec_panel(rob$IrIntens_exRiv, "Transboundary difference") + labs(tag = "A")
e2 <- strip_y(spec_panel(rob$IrGini_exRiv)) + labs(tag = "B")
e3 <- strip_y(spec_panel(rob$GWGini_exRiv)) + labs(tag = "C")
fig_ex <- (e1 | e2 | e3) / make_legend_row() +
  patchwork::plot_layout(heights = c(1, 0.2))
save_pdf(fig_ex, "figure_si_river_exclusion", w = 8.2, h = 4.1)

# =============================================================================
# SI -- BALANCE DIAGNOSTICS (reviewer R3.2)
# =============================================================================
say("figure_si_balance:\n")
all_res <- c(main, rob)
disp_label <- function(r)   # titles collide across intensity/Gini specs
  paste0(r$spec$title,
         ifelse(r$spec$group == "main_gini" & !grepl("Gini", r$spec$title),
                " [Gini]", ""))
smd_df <- do.call(rbind, lapply(all_res, function(r)
  data.frame(spec = disp_label(r), group = r$spec$group,
             mean_smd = r$summary_df$mean_smd, max_smd = r$summary_df$max_smd)))
smd_long <- rbind(
  data.frame(spec = smd_df$spec, stat = "Mean |SMD|", v = smd_df$mean_smd),
  data.frame(spec = smd_df$spec, stat = "Max |SMD|",  v = smd_df$max_smd))
smd_long$spec <- factor(smd_long$spec, levels = rev(unique(smd_df$spec)))

pBal1 <- ggplot(smd_long, aes(v, spec, colour = stat)) +
  geom_vline(xintercept = cfg$balance_tol_substantive, linetype = 2,
             linewidth = 0.35, colour = "grey40") +
  geom_boxplot(outlier.size = 0.35, linewidth = 0.35,
               position = position_dodge(width = 0.65), width = 0.55) +
  scale_colour_manual(NULL, values = c(`Mean |SMD|` = "#0072B2",
                                       `Max |SMD|`  = "#D55E00")) +
  labs(x = sprintf("Post-matching standardized mean difference across %d realizations",
                   cfg$n_realizations),
       y = NULL, title = "Covariate balance across the control ensemble") +
  theme_paper + theme(axis.text.x = element_text(size = 6.5),
                      axis.text.y = element_text(size = 6.5),
                      axis.ticks.x = element_line(linewidth = 0.3),
                      panel.grid.major.x = element_line(colour = "grey92"))

main_names <- names(main)
love <- do.call(rbind, lapply(main_names, function(nm) {
  bt <- main[[nm]]$best$balance
  bt <- bt[bt$covariate != "distance", ]
  lab <- disp_label(main[[nm]])
  rbind(data.frame(spec = lab, covariate = bt$covariate,
                   when = "Before matching", smd = bt$smd_pre),
        data.frame(spec = lab, covariate = bt$covariate,
                   when = "After matching", smd = bt$smd_post))
}))
love$when <- factor(love$when, c("Before matching", "After matching"))
pBal2 <- ggplot(love, aes(smd, covariate, colour = when, shape = when)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_vline(xintercept = c(-1, 1) * cfg$balance_tol_substantive,
             linetype = 2, linewidth = 0.3, colour = "grey40") +
  geom_point(size = 1.5) +
  scale_colour_manual(NULL, values = c(`Before matching` = "#BDBDBD",
                                       `After matching` = "#0072B2")) +
  scale_shape_manual(NULL, values = c(`Before matching` = 1,
                                      `After matching` = 16)) +
  facet_wrap(~spec, nrow = 2, scales = "free_x") +
  labs(x = "Standardized mean difference (best-balanced realization)", y = NULL,
       title = "Covariate balance, best-balanced ensemble member") +
  theme_paper + theme(axis.text.x = element_text(size = 6.5),
                      axis.text.y = element_text(size = 6.5),
                      axis.ticks.x = element_line(linewidth = 0.3),
                      panel.grid.major.x = element_line(colour = "grey92"),
                      strip.text = element_text(size = 6.5, face = "bold"))
fig_bal <- pBal1 / pBal2 + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "A")
save_pdf(fig_bal, "figure_si_balance", w = 7.5, h = 8.2)

say("46_figures done.\n")
