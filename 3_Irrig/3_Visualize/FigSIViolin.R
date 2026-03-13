#MAP: https://code.earthengine.google.com/508d15eb659bdb3f4e44a4bdc3ab93cd

combine_res_mirrored_density_plots <- function(
    res_list,
    sig_level = 0.1,
    ylb = "Gini",
    capt = "",
    adjust = 1,
    n = 256,
    max_half_width = 0.30,
    bar_halfwidth = 0.16,        # longer bars
    bar_lwd = 2,               # thicker bars
    bar_label_minwidth = 0.03,
    bar_gap = 0.0               # vertical separation in units of y_span
) {
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  
  labels <- sapply(res_list, `[[`, "label")
  
  q_safe <- function(x, probs) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(rep(NA_real_, length(probs)))
    as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
  }
  
  dens_df <- function(x, side, x_base, max_half_width = 0.30, adjust = 1, n = 256) {
    x <- x[is.finite(x)]
    if (length(x) < 2 || isTRUE(all.equal(stats::sd(x), 0))) {
      y <- if (length(x) == 0) 0 else x[1]
      return(tibble::tibble(
        side = side,
        y = c(y - 1e-6, y, y + 1e-6),
        x_left = x_base,
        x_right = x_base
      ))
    }
    d <- stats::density(x, adjust = adjust, n = n, na.rm = TRUE)
    dd <- tibble::tibble(y = d$x, dens = d$y) %>%
      mutate(w = max_half_width * dens / max(dens, na.rm = TRUE))
    
    if (side == "Control") {
      dd %>% mutate(side = side, x_left = x_base - w, x_right = x_base)
    } else {
      dd %>% mutate(side = side, x_left = x_base, x_right = x_base + w)
    }
  }
  
  # Build stacked (pos / ns / neg) bar at a given x_center & y_bar
  make_sig_bar <- function(eff, p, x_center, y_bar, y_lab, sig_level,
                           bar_halfwidth, label_minwidth,direction) {
    
    ok_p <- is.finite(p)
    ok_e <- is.finite(eff)
    
    sig <- ok_p & (p < sig_level) & ok_e
    ns  <- ok_p & (p >= sig_level)
    
    p_pos <- mean(sig & (eff > 0), na.rm = TRUE)
    p_neg <- mean(sig & (eff < 0), na.rm = TRUE)
    p_ns  <- mean(ns, na.rm = TRUE)
    
    # exact zero among significant -> treat visually as ns
    p_zero_sig <- mean(sig & (eff == 0), na.rm = TRUE)
    if (!is.na(p_zero_sig) && p_zero_sig > 0) p_ns <- p_ns + p_zero_sig
    
    tot <- p_pos + p_ns + p_neg
    if (!is.finite(tot) || tot <= 0) {
      p_pos <- 0; p_ns <- 1; p_neg <- 0
      tot <- 1
    }
    p_pos <- p_pos / tot
    p_ns  <- p_ns  / tot
    p_neg <- p_neg / tot
    
    if (direction == "right") {
      x0 <- x_center
      x1 <- x_center + 2 * bar_halfwidth
    } else {
      x0 <- x_center - 2 * bar_halfwidth
      x1 <- x_center
    }
    
    x_pos_end <- x0 + (x1 - x0) * p_pos
    x_ns_end  <- x_pos_end + (x1 - x0) * p_ns
    x_neg_end <- x1
    
    seg <- tibble::tibble(
      part = factor(c("pos", "ns", "neg"), levels = c("pos", "ns", "neg")),
      xs = c(x0, x_pos_end, x_ns_end),
      xe = c(x_pos_end, x_ns_end, x_neg_end),
      y = y_bar,
      yend = y_bar
    )
    
    lab <- seg %>%
      mutate(
        xm = (xs + xe) / 2,
        w = (xe - xs),
        lbl = sprintf("%d%%", round(100 * c(p_pos, p_ns, p_neg)))
      ) %>%
      filter(w > label_minwidth & part != "ns")
    
    list(seg = seg, lab = lab)
  }
  
  plots <- imap(res_list, function(res, i) {
    df <- res$summary_df
    stopifnot(all(c("int_eff", "treat_eff", "int_p", "treat_p") %in% names(df)))
    
    control_vals <- df$int_eff
    effect_vals  <- df$treat_eff
    
    # x positions (aligned bases)
    x_center <- 1
    x_nudge <- c(Control = -0.10, Effect = 0.10)
    x_base_control <- x_center + x_nudge["Control"]
    x_base_effect  <- x_center + x_nudge["Effect"]
    
    # densities
    dd_control <- dens_df(control_vals, "Control", x_base_control,
                          max_half_width = max_half_width, adjust = adjust, n = n)
    dd_effect  <- dens_df(effect_vals,  "Effect",  x_base_effect,
                          max_half_width = max_half_width, adjust = adjust, n = n)
    dd <- bind_rows(dd_control, dd_effect) %>%
      mutate(side = factor(side, levels = c("Control", "Effect")))
    
    # summaries
    summ <- tibble::tibble(
      side = factor(c("Control", "Effect"), levels = c("Control", "Effect")),
      mean = c(mean(control_vals, na.rm = TRUE), mean(effect_vals, na.rm = TRUE)),
      q05  = c(q_safe(control_vals, 0.05)[1], q_safe(effect_vals, 0.05)[1]),
      q25  = c(q_safe(control_vals, 0.25)[1], q_safe(effect_vals, 0.25)[1]),
      q50  = c(q_safe(control_vals, 0.50)[1], q_safe(effect_vals, 0.50)[1]),
      q75  = c(q_safe(control_vals, 0.75)[1], q_safe(effect_vals, 0.75)[1]),
      q95  = c(q_safe(control_vals, 0.95)[1], q_safe(effect_vals, 0.95)[1])
    ) %>%
      mutate(xc = ifelse(side == "Control", x_base_control, x_base_effect))
    
    # y-range for bar placement
    y_min <- min(c(control_vals, effect_vals), na.rm = TRUE)
    y_max <- max(c(control_vals, effect_vals), na.rm = TRUE)
    y_span <- y_max - y_min
    if (!is.finite(y_span) || y_span == 0) y_span <- 1
    
    # ---- separate bars vertically ----
    y_bar_eff  <- y_min - 0.14 * y_span
    y_lab_eff  <- y_min - 0.21 * y_span
    
    y_bar_ctrl <- y_min - (0.14 + bar_gap) * y_span
    y_lab_ctrl <- y_min - (0.21 + bar_gap) * y_span
    
    bar_control <- make_sig_bar(
      eff = df$int_eff, p = df$int_p,
      x_center = x_base_control,
      y_bar = y_bar_ctrl, y_lab = y_lab_ctrl,
      sig_level = sig_level,
      bar_halfwidth = bar_halfwidth,
      label_minwidth = bar_label_minwidth,direction='left'
    )
    bar_effect <- make_sig_bar(
      eff = df$treat_eff, p = df$treat_p,
      x_center = x_base_effect,
      y_bar = y_bar_eff, y_lab = y_lab_eff,
      sig_level = sig_level,
      bar_halfwidth = bar_halfwidth,
      label_minwidth = bar_label_minwidth,direction='right'
    )
    
    bar_seg <- bind_rows(
      bar_control$seg %>% mutate(which = "Control"),
      bar_effect$seg  %>% mutate(which = "Effect")
    )
    bar_lab <- bind_rows(
      bar_control$lab %>%
        mutate(which = "Control",
               y = y_lab_ctrl),          # below control bar
      
      bar_effect$lab %>%
        mutate(which = "Effect",y = y_lab_ctrl)
      # y = y_bar_eff + 0.07 * y_span)  # slightly above effect bar
    )
    
    # aesthetics
    fill_vals <- c("Control" = "#FADADD", "Effect" = "#D9A600")
    line_vals <- c("Control" = "#C97C8C", "Effect" = "#A67800")
    
    bar_cols <- c(
      pos = "#2C7FB8", # blue
      ns  = "#BDBDBD", # grey
      neg = "#D7301F"  # red
    )
    
    p <- ggplot() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      
      # densities
      geom_ribbon(
        data = dd,
        aes(y = y, xmin = x_left, xmax = x_right, fill = side),
        alpha = 0.9, color = NA
      ) +
      geom_path(data = dd, aes(x = x_left,  y = y, color = side), linewidth = 0.35) +
      geom_path(data = dd, aes(x = x_right, y = y, color = side), linewidth = 0.35) +
      scale_fill_manual(values = fill_vals, guide = "none") +
      scale_color_manual(values = line_vals, guide = "none") +
      
      # summaries
      geom_linerange(data = summ, aes(x = xc, ymin = q05, ymax = q95, color = side), linewidth = 0.5) +
      geom_linerange(data = summ, aes(x = xc, ymin = q25, ymax = q75, color = side), linewidth = 1.3) +
      geom_point(data = summ, aes(x = xc, y = mean, color = side), size = 2.1) +
      geom_point(data = summ, aes(x = xc, y = q50), shape = 95, size = 8, color = "gray10") +
      
      # two stacked bars (control + effect) with pos/ns/neg segments
      geom_segment(
        data = bar_seg,
        aes(x = xs, xend = xe, y = y, yend = yend, color = part),
        linewidth = bar_lwd,
        lineend = "butt"
      ) +
      scale_color_manual(values = c(line_vals, bar_cols), guide = "none") +
      
      geom_text(
        data = bar_lab,
        aes(x = xm, y = y, label = lbl),
        size = 3.7,
        color = "gray20"
      ) +
      
      coord_cartesian(ylim = c(y_min - (0.28 + bar_gap + 0.06) * y_span, y_max + 0.06 * y_span)) +
      scale_x_continuous(limits = c(0.5, 1.5), breaks = NULL) +
      labs(x = NULL, y = ylb) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", size = 0.6),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    if (!is.null(labels)) p <- p + ggtitle(labels[i])
    p
  })
  
  combined <- wrap_plots(plots, ncol = length(plots)) +
    plot_annotation(caption = capt) &
    theme(plot.caption = element_text(hjust = 0, size = 9, face = "italic"))
  
  print(combined)
  return(combined)
}


# source('Fct.R')
library(qs)
#load existing analyses
res_files <- list.files("resultOut", pattern = "\\.qs$", full.names = TRUE)
res_objects <- lapply(res_files, qs::qread)
names(res_objects) <- gsub("\\.qs$", "", basename(res_files))
list2env(res_objects, envir = .GlobalEnv)


#SI
CropSuitInt$label="Crop Suitability"
IrNeedInt$label="Irrigation Need"
Overdraft_Irrig$label="Overdraft | Irrig"
SI_A=combine_res_mirrored_density_plots(list(CropSuitInt,IrNeedInt,Overdraft_Irrig),ylb='Irrig. Intens.')

###OTHER SIs
IrGiniDoubleRobust$label="Double Robust"
ItIntdoublerobust$label="Double Robust"
Overdraft6$label="Overdraft\n(>6 months)"
SI_B1=combine_res_mirrored_density_plots(list(ItIntdoublerobust,Overdraft6),ylb='Irrig. Intens.')
SI_B2=combine_res_mirrored_density_plots(list(IrGiniDoubleRobust),ylb='Gini')
SI_B=SI_B1/SI_B2



