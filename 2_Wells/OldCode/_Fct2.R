
library(splines)
library(MASS)
library(tidyverse)
library(sandwich)
library(MatchIt)
library(cobalt)
library(metafor)
library(furrr)
library(clubSandwich)
library(lme4)
library(performance)
library(stringr)
library(knitr)
library(kableExtra)
library(patchwork)

# 
# analyze_one <- function(dat,Aquifer,outdir = "aqf_plt",trendVars = c(character(0)),alpha = 0.1,plt = FALSE,nM=1,knn=F) {
#   # --- Ensure required packages are available in all workers ---
#   requireNamespace("dplyr")
#   requireNamespace("sandwich")
#   requireNamespace("splines")
#   requireNamespace("MASS")
#   requireNamespace("ggplot2")
#   
#   ##Housekeeping
#   #Drop Nas
#   dat <- dat %>%
#     dplyr::filter(!is.na(GWSlp), !is.na(dist_LB_km))
#   #Drop missing covariates
#   trendVars <- intersect(trendVars, names(dat))
#   #Drop covariates with no variations
#   trendVars <- trendVars[vapply(trendVars, function(v) {
#     x <- dat[[v]]
#     if (is.null(x)) return(FALSE)
#     if (all(is.na(x))) return(FALSE)
#     dplyr::n_distinct(x, na.rm = TRUE) > 1
#   }, logical(1))]
#   
#   #Get covariates as factors
#   dat <- dat %>%
#     dplyr::mutate(across(
#       where(is.logical),
#       ~ factor(.x, levels = c(FALSE, TRUE))
#     ))
#   
#   
#   dat$w_clust=1
#   if(knn){
#   dat$w_clust=knn_reg_weights(dat)
#   }
#   
#   
#   ###########MODELS##############
#   # --- Linear model (no dist_LB_km) ---
#   f_lin_z <- reformulate(trendVars, "GWSlp")
#   m_lin_z <- stats::lm(f_lin_z, data = dat,weights = w_clust)
#   se_hc3_z <- sqrt(diag(get_vcov_safe(m_lin_z)))
#   beta_z_lin <- stats::coef(m_lin_z)["(Intercept)"]
#   se_z_lin <- se_hc3_z["(Intercept)"]
#   t_z_lin <- beta_z_lin / se_z_lin
#   df_z <- stats::nobs(m_lin_z) - length(stats::coef(m_lin_z))
#   p_z_lin <- 2 * stats::pt(-abs(t_z_lin), df_z)
#   
#   # --- Linear model ---
#   f_lin <- reformulate(c("dist_LB_km", trendVars), "GWSlp")
#   m_lin <- stats::lm(f_lin, data = dat,weights = w_clust)
#   se_hc3 <- sqrt(diag(get_vcov_safe(m_lin)))
#   beta_lin <- stats::coef(m_lin)["dist_LB_km"]
#   se_lin <- se_hc3["dist_LB_km"]
#   t_lin <- beta_lin / se_lin
#   df <- stats::nobs(m_lin) - length(stats::coef(m_lin))
#   p_lin <- 2 * stats::pt(-abs(t_lin), df)
#   n_unique_dist <- length(unique(stats::na.omit(dat$dist_LB_km)))
#   if (n_unique_dist <= 3) {
#     p_lin <- NA_real_
#     shape <- "degenerate (too few x points)"
#     direction <- NA_character_
#   }
#   
# 
#   ######################################################### --- Plot (optional) --- NO COVARIATE ADJUSTMENT.
#   if (plt && nrow(dat) > nM) {
#     
#     dir.create(outdir, showWarnings = FALSE)
#     
#     dat$y <- stats::model.response(stats::model.frame(m_lin))
#     grid <- data.frame(
#       dist_LB_km = seq(min(dat$dist_LB_km, na.rm = TRUE),
#                        max(dat$dist_LB_km, na.rm = TRUE),
#                        length.out = 100)
#     )
#     grid$pred_lin <- predict(m_lin, newdata = grid)
#     
#     p <- ggplot2::ggplot(dat, ggplot2::aes(dist_LB_km, y)) +
#       ggplot2::geom_point(alpha = 0.3, size = 1) +
#       ggplot2::geom_line(data = grid, ggplot2::aes(y = pred_lin),
#                          color = "blue", linetype = ifelse(p_lin < alpha, "solid", "dashed")) +
#       ggplot2::labs(
#         x = "Distance to border (km)", y = "GWSlp",
#         title = paste0("Aquifer ", Aquifer, ", TB: ", unique(dat$TB)),
#         subtitle = sprintf("b=%.2f | p=%.3f | n=%d", beta_lin, p_lin, nrow(dat))
#       ) +
#       ggplot2::theme_classic(base_size = 10)
#     
#     ggplot2::ggsave(
#       file.path(outdir, paste0("aquifer_", Aquifer, "_TB", unique(dat$TB), ".pdf")),
#       p, width = 6, height = 4, dpi = 150
#     )
#   }
#   
#   # --- Output summary row ---
#   dplyr::tibble(
#     n = nrow(dat),
#     beta_0 = beta_z_lin, se_0 = se_z_lin, p_0 = p_z_lin,
#     beta_1 = beta_lin, se_1 = se_lin, p_1 = p_lin
#   )
#   
# }



analyze_one_robust <- function(dat,trendVars = c(character(0))) {
    # --- Ensure required packages are available in all workers ---
    requireNamespace("dplyr")
    requireNamespace("sandwich")
    requireNamespace("splines")
    requireNamespace("MASS")
    requireNamespace("ggplot2")

    ##Housekeeping
    #Drop Nas
    dat <- dat %>%
      dplyr::filter(!is.na(GWSlp), !is.na(dist_LB_km))
    #Drop missing covariates
    trendVars <- intersect(trendVars, names(dat))
    #Drop covariates with no variations
    trendVars <- trendVars[vapply(trendVars, function(v) {
      x <- dat[[v]]
      if (is.null(x)) return(FALSE)
      if (all(is.na(x))) return(FALSE)
      dplyr::n_distinct(x, na.rm = TRUE) > 1
    }, logical(1))]
    
    #Get covariates as factors
    dat <- dat %>%
      dplyr::mutate(across(
        where(is.logical),
        ~ factor(.x, levels = c(FALSE, TRUE))
      ))
    
    
  ###########MODELS##############
  # --- Linear model (no dist_LB_km) ---
  f_lin_z <- reformulate(trendVars, "GWSlp")
  # --- Robust fit (MASS::rlm, no dist_LB_km) ---
  X_z <- stats::model.matrix(f_lin_z, dat)
  df_z <- nrow(X_z)-qr(X_z)$rank
 
  # Remove constant covariates
  const_vars_z <- setdiff(colnames(X_z)[apply(X_z, 2, sd, na.rm = TRUE) == 0], "(Intercept)")
  if (length(const_vars_z) > 0) {
    vars_z <- setdiff(all.vars(f_lin_z)[-1], const_vars_z)
    f_lin_z <- reformulate(vars_z, "GWSlp")
    X_z <- stats::model.matrix(f_lin_z, dat)
  }
  
  # Address rank deficiencies / collinearities (remove the relevant rows.)
  qr_z <- qr(X_z)
  if (qr_z$rank < ncol(X_z)) {
    aliased_z <- setdiff(colnames(X_z),
                         colnames(X_z)[qr_z$pivot[seq_len(qr_z$rank)]])
    aliased_z <- setdiff(aliased_z, "(Intercept)")
    if (length(aliased_z) > 0) {
      vars_z <- setdiff(all.vars(f_lin_z)[-1], aliased_z)
      f_lin_z <- reformulate(vars_z, "GWSlp")
      X_z <- stats::model.matrix(f_lin_z, dat)
    }
  }
  
  # Fit robust model
  if (qr(X_z)$rank < ncol(X_z)) {
    m_z_rlm <- NULL
  } else {
    m_z_rlm <- tryCatch(MASS::rlm(f_lin_z, data = dat, maxit = 50),
                        error = function(e) NULL)
  }
  
  # Extract coefficients
  beta_0 <- se_0 <- p_0 <- NA_real_
  if (!is.null(m_z_rlm)) {
    sr_z <- summary(m_z_rlm)$coef
    if ("(Intercept)" %in% rownames(sr_z)) {
      beta_0 <- sr_z["(Intercept)", "Value"]
      se_0   <- sr_z["(Intercept)", "Std. Error"]
      p_0    <- 2 * stats::pt(-abs(sr_z["(Intercept)", "t value"]), df_z)
    }
  }
  
  
  # --- Linear model ---
  f_lin <- reformulate(c("dist_LB_km", trendVars), "GWSlp")
  # --- Robust fit (MASS::rlm) ---
  X <- stats::model.matrix(f_lin, dat)
  df_1 <- nrow(X)-qr(X)$rank
  ##Addresses constant covariates, which are automatically dealt with in lm but not rlm
  const_vars <- setdiff(colnames(X)[apply(X, 2, sd, na.rm = TRUE) == 0], "(Intercept)")
  if (length(const_vars) > 0) {
    vars <- setdiff(all.vars(f_lin)[-1], const_vars)
    f_lin <- reformulate(vars, "GWSlp")
    X <- stats::model.matrix(f_lin, dat)
  }
  ##Addresses rank deficiencies /colinearities (again automatically dealt with in lm but no rlm)
  qr_obj <- qr(X)
  if (qr_obj$rank < ncol(X)) {
    aliased <- setdiff(colnames(X),
                       colnames(X)[qr_obj$pivot[seq_len(qr_obj$rank)]])
    aliased <- setdiff(aliased, "(Intercept)")
    if (length(aliased) > 0) {
      vars <- setdiff(all.vars(f_lin)[-1], aliased)
      f_lin <- reformulate(vars, "GWSlp")
      X <- stats::model.matrix(f_lin, dat)
    }
  }
  ##Finally fits the model
  if (qr(X)$rank < ncol(X)) {
    m_rlm <- NULL
  } else {
    m_rlm <- tryCatch(MASS::rlm(f_lin, data = dat, maxit = 50),
                      error = function(e) NULL)
  }
  #Assign estimated values
  beta_1 <- se_1 <- p_1 <- NA_real_
  if (!is.null(m_rlm)) {
    sr <- summary(m_rlm)$coef
    if ("dist_LB_km" %in% rownames(sr)) {
      beta_1 <- sr["dist_LB_km", "Value"]
      se_1 <- sr["dist_LB_km", "Std. Error"]
      p_1 <- 2 * stats::pt(-abs(sr["dist_LB_km", "t value"]), df_1)
    }
  }
  
  
  # --- Output summary row ---
  dplyr::tibble(
    n = nrow(dat),
    beta_0,se_0,p_0,beta_1,se_1,p_1
  )
}

# dclust = NULL,         # NEW: radius in km, or NULL
# kde = FALSE,           # NEW: if TRUE use KDE weights; default fixed-radius
# lon_col = "lon",       # NEW: coords for weights
run=function(plot=F,outdir='aqf_plt',tVar=c(character(0)),alph=0.1,match='full',cov_match =c("lat_c","lon_c","urbkHaKm2",'CS_max','LB_river'),rob=F,matchWg=T,nMin=20,shp=NULL, loo=NULL,r=1,dropTopSE=F,noRivWells=F,noAqfWells=F,dclust=NULL,kde=F){
  
  #parallel settings  
  plan(multisession, workers = parallel::detectCores() - 1)

  
  ##############################
  #Fetch data
  ############################## 
  df=read.csv('wellsData.csv')%>%
    mutate(LB_river=(LB_river==1))%>%filter(!CC%in%c("AUS", "NZL","TWN","JPN", "IRL", "GBR"))#%>% #remove island nations
    # mutate(arind=pet_mm/prec_mm)%>%filter(pet_mm>0,prec_mm>0) #aridity index
  
  if(noRivWells){
    df=df%>%filter(!LB_river)
  }
  
  if(noAqfWells){
    df <- df %>%
      group_by(Aquifer) %>%
      filter(!any(LB_river)) %>%
      ungroup()
  }
  #######################################
  #aquifer level trend analysis
  #######################################
  if(rob){
    if(plot) warning('no plotting possible for robust lm')
    results <- df %>%
      group_split(Aquifer, CC,.keep = TRUE) %>%
      future_map_dfr(~ {
        y <- dplyr::slice_head(.x, n = 1)
        out <- analyze_one_robust(
          .x,
          trendVars =tVar
        )
        # Add grouping info to result
        out$Aquifer <- y$Aquifer
        out$CC <- y$CC
        out
      }, .progress = TRUE)
  }else{
    results <- df %>%
      group_split(Aquifer, CC,.keep = TRUE) %>%
      future_map_dfr(~ {
        y <- dplyr::slice_head(.x, n = 1)
        out <- analyze_one(
          .x,
          Aquifer = paste0(y$Aquifer, "_", y$CC),
          outdir=outdir,
          plt=plot,
          trendVars =tVar,
          nM=nMin,dclust=dclust,kde=kde
        )
        # Add grouping info to result
        out$Aquifer <- y$Aquifer
        out$CC <- y$CC
        out
      }, .progress = TRUE)
  }
  
  
  aqf=results%>%
    left_join(read.csv('aqfData.csv'))%>%
    mutate(unit_id = paste0(Aquifer, "_", CC))%>%
    left_join(df%>%group_by(Aquifer,CC)%>%dplyr::summarize(LB_river=mean(LB_river))) #propensity that closest border is river.
  
  
  aqf=aqf%>%filter(n>nMin)%>%
    filter(!is.na(beta_0), !is.na(se_0), se_0 > 0,
           !is.na(beta_1), !is.na(se_1), se_1 > 0)%>%
    filter(!is.na(TB)) %>%
    filter(!if_any(all_of(cov_match), ~ is.na(.x))) 
  

  
  ##############################
  #MISC ROB CHECKS (remove if irrrelevant)
  ##############################
  if(!is.null(shp)){
    aqf=aqf%>%filter(shape%in%shp)
  }
  if(!is.null(loo)){ aqf=aqf%>%filter(!Aquifer%in%loo)}
  #drop highest weighted aqfs (ROBUSTNMESS ONLY)
  if(dropTopSE){
  aqf=aqf%>%filter(!Aquifer%in%c("Biscayne Aquifer","Central Minnesota Surficial and Buried Sand and Gravel Aquifers","East Bay Plain","Santa Clara Valley", "Orange County Coastal Plain" ,"Northern High Plains","Dakota Aquifer System","Delmarva Peninsula","Biscayne Aquifer", "Central Mississippi Embayment"))
  }

  ##############################
  #Matching on covariates
  ##############################
  form_match <- reformulate(cov_match, response = "TB")
  m.out <- matchit(form_match, data = aqf, method = match,ratio=r)
  aqf_full <- match.data(m.out)                  # adds `weights` and restrict to matches
  divisor <- if (matchWg) {
    aqf_full$weights
  } else {
    1
  }
  aqf_full <- aqf_full %>%
    dplyr::mutate(
      # Combine weights
      vi_0 = winsor((se_0)^2 / divisor),
      vi_1 = winsor((se_1)^2 / divisor)
    )


  ##############################
  #Fit metamodel
  ##############################
  fit_mv_0 <- tryCatch(
      metafor::rma.mv(
        yi = beta_0, V = vi_0, mods = ~ 1 + TB,
        random = ~ 1 | CC,
        data = aqf_full, method = "REML",
        control = list(optimizer = "optim", rel.tol = 1e-8)
      ),
      error = function(e) {
        warning("⚠️ MV REML model failed to converge — retrying with ML.")
        metafor::rma.mv(
          yi = beta_0, V = vi_0, mods = ~ 1 + TB,
          random = ~ 1 | CC,
          data = aqf_full, method = "ML"
        )
      })
  
  fit_mv_1 <- tryCatch(
          metafor::rma.mv(
            yi = beta_1, V = vi_1, mods = ~ 1 + TB,
            random = ~ 1 | CC,
            data = aqf_full, method = "REML",
            control = list(optimizer = "optim", rel.tol = 1e-8)
          ),
          error = function(e) {
            warning("⚠️ MV REML model failed to converge — retrying with ML.")
            metafor::rma.mv(
              yi = beta_1, V = vi_1, mods = ~ 1 + TB,
              random = ~ 1 | CC,
              data = aqf_full, method = "ML"
            )
          })
    
  
    out=dplyr::tibble(
    b_0  = as.numeric(coef(fit_mv_0)),
    se_0 = sqrt(diag(vcov(fit_mv_0))),
    p_0  = 2 * pnorm(abs(b_0 / se_0), lower.tail = FALSE),
    nm_0 = names(coef(fit_mv_0)),
    st0= getstars(p_0),
    tau2_CC_0  = fit_mv_0$sigma2[1],
    
    b_1  = as.numeric(coef(fit_mv_1)),
    se_1 = sqrt(diag(vcov(fit_mv_1))),
    p_1  = 2 * pnorm(abs(b_1 / se_1), lower.tail = FALSE),
    nm_1 = names(coef(fit_mv_1)),
    st1= getstars(p_1),
    tau2_CC_1  = fit_mv_1$sigma2[1],
    n         = nrow(aqf_full),
    k_cluster = dplyr::n_distinct(aqf_full$CC))

    
  
  future::plan(sequential)
  
  
  
  # --- Diagnosis plots ---------------------------------------------
  if (plot) {
    message("Saving matching diagnostics to ", outdir, " ...")
    
    f_love   = file.path(outdir, "__MATCH_loveplot.pdf")
    f_heat   <- file.path(outdir, "__MATCH_balance_heatmap.pdf")
    f_pshist <- file.path(outdir, "__MATCH_ps_histogram.pdf")
    
    # Love plot (ggplot2 object)
    p1 <- love.plot(m.out, abs = TRUE, thresholds = c(m = .1))
    ggsave(f_love, plot = p1, width = 7, height = 5)
    
    # Balance heatmap (ggplot2 object)
    suppressWarnings({
      bt <- bal.tab(m.out, un = TRUE, disp.means = TRUE)
      p2 <- plot(bt, type = "heat", var.order = "unadjusted")
      ggsave(f_heat, plot = p2, width = 7, height = 5)
    })
    
    # Propensity score histogram (base R plot)
    pdf(f_pshist, width = 7, height = 5)
    plot(m.out, type = "hist")
    dev.off()
    
    # --- Proportion summaries with color intensity ---
    prop_tbl <- aqf_full %>%
      mutate(
        # classify sign and significance
        sign = case_when(
          beta_1 > 0 ~ 1,
          beta_1 < 0 ~ -1,
          TRUE ~ 0
        ),
        # compute color value: direction × significance strength
        # p_lin capped at 0.3 so nonsignificant → gray
        color_val = sign * (1 - pmin(p_0 / 0.3, 1)),
        Res = case_when(
          p_1 > alph &beta_1<0 ~ "neg_ns",
          p_1 > alph  ~ "pos_ns",
          beta_1 < 0 ~ "neg",
          beta_1 > 0 ~ "pos",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(Res)) %>%
      group_by(TB, Res) %>%
      summarise(
        n = n(),
        w_sum = sum(weights, na.rm = TRUE),
        mean_color = mean(color_val, na.rm = TRUE), # <- average color intensity
        .groups = "drop"
      ) %>%
      group_by(TB) %>%
      mutate(
        prop_unweighted = n / sum(n),
        prop_weighted   = w_sum / sum(w_sum),
        total_n         = sum(n)
      ) %>%
      ungroup()
    
    # --- Build unified plotting table ---
    plot_tbl <- bind_rows(
      prop_tbl %>%
        filter(!TB) %>%
        transmute(Group = "Domestic", Res, Proportion = prop_unweighted, total_n, mean_color),
      prop_tbl %>%
        filter(!TB) %>%
        transmute(Group = "Domestic (matched)", Res, Proportion = prop_weighted, total_n, mean_color),
      prop_tbl %>%
        filter(TB) %>%
        transmute(Group = "Transboundary", Res, Proportion = prop_weighted, total_n, mean_color)
    ) %>%
      mutate(
        Res = factor(Res, levels = c("neg", "neg_ns","pos_ns", "pos")),
        Group = factor(Group, levels = c("Domestic", "Domestic (matched)", "Transboundary"))
      )%>%filter(Group!='Domestic')
    cols <- c(
      "pos"    = "#5D7DA9",  # strong blue (positive, significant)
      "pos_ns" = "#AEBBD3",  # light blue (positive, non-significant)
      "neg"    = "#CC513C",  # strong red (negative, significant)
      "neg_ns" = "#E6A59C"   # light red (negative, non-significant)
    )
    # cols <- c(
    #   "pos"    = "#5D7DA9",
    #   "pos_ns" = "#C0C9DA",
    #   "neg"    = "#CC513C",
    #   "neg_ns" = "#F1BBB3"
    # )
    # --- Plot with gradient fill ---
    p1 <- ggplot(plot_tbl, aes(x = Group, y = Proportion, fill = Res)) +
      geom_col(width = 0.7) +
      scale_fill_manual(values = cols, name = "Response") +
      labs(
        x = NULL,
        y = "Proportion",
        title = "Distance to border depletion trend"
      ) +
      theme_classic(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "top"
      )
    
    
    prop_tbl <- aqf_full %>%
      mutate(
        # classify sign and significance
        sign = case_when(
          beta_0 > 0 ~ 1,
          beta_0 < 0 ~ -1,
          TRUE ~ 0
        ),
        # compute color value: direction × significance strength
        # p_lin capped at 0.3 so nonsignificant → gray
        color_val = sign * (1 - pmin(p_0 / 0.3, 1)),
        Res = case_when(
          p_0 > alph &beta_0<0 ~ "neg_ns",
          p_0 > alph  ~ "pos_ns",
          beta_0 < 0 ~ "neg",
          beta_0 > 0 ~ "pos",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(Res)) %>%
      group_by(TB, Res) %>%
      summarise(
        n = n(),
        w_sum = sum(weights, na.rm = TRUE),
        mean_color = mean(color_val, na.rm = TRUE), # <- average color intensity
        .groups = "drop"
      ) %>%
      group_by(TB) %>%
      mutate(
        prop_unweighted = n / sum(n),
        prop_weighted   = w_sum / sum(w_sum),
        total_n         = sum(n)
      ) %>%
      ungroup()
    
    # --- Build unified plotting table ---
    plot_tbl <- bind_rows(
      prop_tbl %>%
        filter(!TB) %>%
        transmute(Group = "Domestic", Res, Proportion = prop_unweighted, total_n, mean_color),
      prop_tbl %>%
        filter(!TB) %>%
        transmute(Group = "Domestic (matched)", Res, Proportion = prop_weighted, total_n, mean_color),
      prop_tbl %>%
        filter(TB) %>%
        transmute(Group = "Transboundary", Res, Proportion = prop_weighted, total_n, mean_color)
    ) %>%
      mutate(
        Res = factor(Res, levels = c("neg", "neg_ns","pos_ns", "pos")),
        Group = factor(Group, levels = c("Domestic", "Domestic (matched)", "Transboundary"))
      )%>%filter(Group!='Domestic')
    cols <- c(
      "pos"    = "#5D7DA9",  # strong blue (positive, significant)
      "pos_ns" = "#AEBBD3",  # light blue (positive, non-significant)
      "neg"    = "#CC513C",  # strong red (negative, significant)
      "neg_ns" = "#E6A59C"   # light red (negative, non-significant)
    )
    # cols <- c(
    #   "pos"    = "#5D7DA9",
    #   "pos_ns" = "#C0C9DA",
    #   "neg"    = "#CC513C",
    #   "neg_ns" = "#F1BBB3"
    # )
    # --- Plot with gradient fill ---
    p0 <- ggplot(plot_tbl, aes(x = Group, y = Proportion, fill = Res)) +
      geom_col(width = 0.7) +
      scale_fill_manual(values = cols, name = "Response") +
      labs(
        x = NULL,
        y = "Proportion",
        title = "Average depletion"
      ) +
      theme_classic(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "top"
      )
    # --- Optional save ---
    ggsave(
      file.path(outdir, "__MATCH_proportions.pdf"),
      plot = p0+p1,
      width = 5, height = 4, dpi = 300
    )
    
    
  }
  return(list(aqf_full,out))
}


#######################################Helper Fct
winsor <- function(x, p = c(0.01, 0.99)) {
  qs <- quantile(x, p, na.rm = TRUE)
  pmin(pmax(x, qs[1]), qs[2])
}


getstars=function(p){
  
  dplyr::case_when(
    p < 0.01 ~ "***",
    p < 0.05  ~ "**",
    p < 0.1  ~ "*",
    TRUE      ~ ""
  )
}

get_vcov_safe <- function(model) {
  out <- tryCatch(sandwich::vcovHC(model, type = "HC3"), error = function(e) NULL)
  if (is.null(out) || any(!is.finite(diag(out)))) {
    out <- sandwich::vcovHC(model, type = "HC1")
  }
  out
}


set.seed(123)

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
  ggplot(plot_df, aes(x = component, y = est, fill = component)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  width = 0.15, linewidth = 0.6) +
    geom_point(size = 3, shape = 21, color = "black") +
    geom_text(aes(label = stars, y = ci_high * 1.08),
              vjust = 0, size = 5) +
    scale_fill_manual(values = c("Intercept" = "white", "TB effect" = "black")) +
    facet_wrap(~ outcome, scales = "free_y") +
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
    )
}

library(dplyr)
library(ggplot2)


plot_binned_depletion <- function(df, x, bin_size = 20, xmx = 500, 
                                  w = "w_final", ci = TRUE, mnLn = TRUE) {
  library(dplyr)
  library(ggplot2)
  # Count wells and aquifers per TB
  dflab <- x[[1]]%>%
    group_by(TB) %>%
    summarise(
      n_wells = sum(n),
      n_aq = n(),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0("Wells: ", n_wells, "\nSegments: ", n_aq)
    )
  

  # Extract tau2_CC
  tauCC <- x[[2]] %>%
    summarise(tau2_CC = mean(c(tau2_CC_0, tau2_CC_1), na.rm = TRUE)) %>%
    pull(tau2_CC)
  
  # Compute composite weights and bins
  df <- df %>%filter(dist_LB_km < xmx)%>%
    mutate(no = 1) %>%
    left_join(x[[1]] %>% dplyr::select(Aquifer, CC, vi_1, weights), by = c("Aquifer", "CC")) %>%
    mutate(
      vi_1 = ifelse(is.na(vi_1) | vi_1 <= 0, NA_real_, vi_1),
      weights = ifelse(is.na(weights) | weights < 0, 0, weights),
      w_final = weights / (vi_1 + tauCC)
    ) %>%
    mutate(bin = floor(dist_LB_km / bin_size) * bin_size) %>%
    filter(!is.na(w_final))
  
  # Binned summary
  dfx <- df %>%
    group_by(TB, bin) %>%
    summarise(
      n = sum(!is.na(GWSlp)),
      w_sum = sum(.data[[w]], na.rm = TRUE),
      mean_GWSlp = weighted.mean(GWSlp, .data[[w]], na.rm = TRUE),
      se = sqrt(sum(.data[[w]] * (GWSlp - mean_GWSlp)^2, na.rm = TRUE)) / sqrt(n),
      ci_low = mean_GWSlp - 1.645 * se,
      ci_high = mean_GWSlp + 1.645 * se,
      .groups = "drop"
    ) %>%
    group_by(TB) %>%
    mutate(
      w_sum_m = (w_sum - min(w_sum, na.rm = TRUE)) /
        (max(w_sum, na.rm = TRUE) - min(w_sum, na.rm = TRUE))
    ) %>%
    ungroup()
  
  # Mean per TB
  dfw <- df %>%
    group_by(TB) %>%
    summarise(mn = weighted.mean(GWSlp, .data[[w]], na.rm = TRUE), .groups = "drop")
  
  
  # Order facets (TRUE first)
  dfx$TB <- factor(dfx$TB, levels = c(TRUE, FALSE))
  
  # Custom facet labels
  tb_labels <- c(
    `TRUE` = "Transboundary",
    `FALSE` = "Non-Transboundary"
  )
  
  # Build plot
  p <- ggplot(data = dfx, aes(x = bin + bin_size / 2, y = mean_GWSlp, color = TB, fill = TB))
  
  if (ci) {
    p <- p + geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2, color = NA)
  }
  
  p <- p +
    geom_point(aes(alpha = w_sum_m)) +
    geom_line(aes(alpha = w_sum_m)) +
    scale_color_manual(values = c("TRUE" = "#1F77B4", "FALSE" = "#E64B35")) +
    scale_fill_manual(values = c("TRUE" = "#1F77B4", "FALSE" = "#E64B35"))
  
  if (mnLn) {
    p <- p +
      geom_hline(data = dfw, aes(yintercept = mn, color = TB),
                 linewidth = 0.5, linetype = "dashed") +
      geom_text(
        data = dfw,
        aes(x = 0, y = mn, label = sprintf("%.1f", mn), color = TB),
        hjust = -12, vjust = -0.2, size = 4, show.legend = FALSE
      )
  }
  
  # Add text labels for wells and aquifers
  p <- p +
    geom_text(
      data = dflab,
      aes(x = Inf, y = Inf, label = label),
      color = "black",
      hjust = 1.1, vjust = 1.1,
      size = 3.5,
      inherit.aes = FALSE
    ) +
    labs(
      x = "Distance to border (km)",
      y = "Mean groundwater depletion (GWSlp)"
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none") +
    facet_wrap(~TB, ncol = 1, scales = "free", labeller = as_labeller(tb_labels))
  
  p
}
make_regtable <- function(model_list, reg = 0, caption = "Regression results", fmt = "latex") {
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(knitr)
  library(kableExtra)
  
  b_col   <- paste0("b_",  reg)
  se_col  <- paste0("se_", reg)
  p_col   <- paste0("p_",  reg)
  nm_col  <- paste0("nm_", reg)
  st_col  <- paste0("st",  reg)
  tau_col <- paste0("tau2_CC_", reg)
  
  need_cols <- c(nm_col, b_col, se_col, p_col, st_col, tau_col, "n", "k_cluster", "lab")
  
  # Combine all list elements (force dplyr::select to avoid masking issues)
  df <- purrr::map_dfr(model_list, ~ .x %>%
                         dplyr::select(dplyr::all_of(need_cols)))
  
  # Rename for simplicity
  names(df) <- c("term", "b", "se", "p", "stars", "tau2", "n", "k_cluster", "label")
  
  fmt_num <- function(x) ifelse(is.na(x), "",
                                ifelse(abs(x) < 0.001, "$<$0.001", sprintf("%.3f", x)))
  fmt_tau <- function(x) ifelse(is.na(x), "", sprintf("%.1f", x))
  
  df_main <- df %>%
    dplyr::filter(term %in% c("intrcpt", "TBTRUE")) %>%
    dplyr::mutate(
      term = dplyr::recode(term, intrcpt = "Intercept", TBTRUE = "Transboundary"),
      entry = sprintf("%s (%s) [%s]%s", fmt_num(b), fmt_num(se), fmt_num(p), stars)
    ) %>%
    dplyr::select(term, label, entry) %>%
    tidyr::pivot_wider(names_from = label, values_from = entry)
  
  df_stats <- df %>%
    dplyr::distinct(label, tau2, n, k_cluster) %>%
    dplyr::mutate(
      tau2 = fmt_tau(tau2),
      n = as.character(n),
      k_cluster = as.character(k_cluster)
    ) %>%
    tidyr::pivot_longer(cols = c(tau2, n, k_cluster), names_to = "term", values_to = "entry") %>%
    dplyr::mutate(
      term = dplyr::recode(term, tau2 = "Tau^2", n = "Observations", k_cluster = "Clusters")
    ) %>%
    tidyr::pivot_wider(names_from = label, values_from = entry)
  
  df_main  <- df_main  %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  df_stats <- df_stats %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  
  df_out <- dplyr::bind_rows(df_main, df_stats)
  
  model_labels <- names(df_out)[-1]
  spec_nums <- paste0("(", seq_along(model_labels), ")")
  
  knitr::kable(df_out, format = fmt, booktabs = TRUE, escape = FALSE,
               caption = caption, align = "lcccc",
               col.names = c("", model_labels)) %>%
    kableExtra::add_header_above(c(" " = 1, stats::setNames(rep(1, length(model_labels)), spec_nums))) %>%
    kableExtra::kable_styling(latex_options = c("hold_position", "scale_down")) %>%
    kableExtra::row_spec(nrow(df_main) + 1, extra_latex_after = "\\midrule")
}

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
# ------------------------------------------------------------
# Declustering weights with a fixed radius of influence (r_km)
# - default: fixed-radius count weights (kde = FALSE)
# - optional: kernel density declustering with bandwidth = r_km (kde = TRUE)
# ------------------------------------------------------------
radius_reg_weights <- function(dat,
                               r_km,
                               kde = FALSE,
                               lon = "lon",
                               lat = "lat",
                               crs_in = 4326,
                               crs_proj = 3857,
                               include_self = TRUE,   # only for fixed-radius count
                               kernel = c("gaussian", "epanechnikov"),
                               w_floor = 1e-12,
                               w_cap_q = 0.99,
                               normalize = TRUE) {
  
  stopifnot(is.data.frame(dat))
  if (!is.numeric(r_km) || length(r_km) != 1 || !is.finite(r_km) || r_km <= 0) {
    stop("r_km must be a single finite number > 0.")
  }
  
  requireNamespace("sf")
  
  n <- nrow(dat)
  if (n == 0) return(numeric(0))
  if (n < 3) return(rep(1, n))
  
  if (!all(c(lon, lat) %in% names(dat))) stop("lon/lat columns not found in dat.")
  if (any(is.na(dat[[lon]]) | is.na(dat[[lat]]))) stop("Missing lon/lat values.")
  
  # project to meters
  pts <- sf::st_as_sf(dat, coords = c(lon, lat), crs = crs_in, remove = FALSE)
  pts_m <- sf::st_transform(pts, crs_proj)
  
  r_m <- r_km * 1000
  
  if (!kde) {
    # ---- Fixed-radius count declustering: w_i ∝ 1 / n_i(r) ----
    neigh <- sf::st_is_within_distance(pts_m, pts_m, dist = r_m)
    n_in <- lengths(neigh)
    if (!include_self) n_in <- pmax(n_in - 1L, 1L)
    
    w <- 1 / n_in
    
  } else {
    # ---- Kernel density declustering: w_i ∝ 1 / fhat_i (bandwidth = r_km) ----
    kernel <- match.arg(kernel)
    
    xy <- sf::st_coordinates(pts_m)
    dx <- outer(xy[, 1], xy[, 1], "-")
    dy <- outer(xy[, 2], xy[, 2], "-")
    d  <- sqrt(dx^2 + dy^2)
    
    if (kernel == "gaussian") {
      # Gaussian kernel: exp(-0.5 (d/h)^2)
      K <- exp(-0.5 * (d / r_m)^2)
      # density up to constant is fine for weights; include constant for stability if you want
      fhat <- rowSums(K) / (2 * pi * r_m^2)
    } else {
      # Epanechnikov (compact support): K(u)=0.75*(1-u^2) for u<=1 else 0
      u <- d / r_m
      K <- 0.75 * (1 - u^2)
      K[u > 1] <- 0
      # constant differs; irrelevant for weights, but keep positive scaling
      fhat <- rowSums(K) / (pi * r_m^2)
    }
    
    fhat <- pmax(fhat, w_floor)
    w <- 1 / fhat
  }
  
  # cap extremes + normalize
  wcap <- stats::quantile(w, probs = w_cap_q, na.rm = TRUE, type = 7)
  w <- pmin(w, as.numeric(wcap))
  
  if (normalize) w <- w / mean(w, na.rm = TRUE)
  
  w
}


# ------------------------------------------------------------
# analyze_one: replace knn with dclust (radius in km)
# - dclust = NULL (default) => no declustering weights (w=1)
# - dclust = number => apply radius_reg_weights(dat, r_km = dclust, kde = kde)
# ------------------------------------------------------------
analyze_one <- function(dat,
                        Aquifer,
                        outdir = "aqf_plt",
                        trendVars = c(character(0)),
                        alpha = 0.1,
                        plt = FALSE,
                        nM = 1,
                        dclust = NULL,         # NEW: radius in km, or NULL
                        kde = FALSE,           # NEW: if TRUE use KDE weights; default fixed-radius
                        lon_col = "lon",       # NEW: coords for weights
                        lat_col = "lat") {
  
  requireNamespace("dplyr")
  requireNamespace("sandwich")
  requireNamespace("splines")
  requireNamespace("MASS")
  requireNamespace("ggplot2")
  
  # Drop NAs
  dat <- dat %>%
    dplyr::filter(!is.na(GWSlp), !is.na(dist_LB_km))
  
  # Drop missing covariates
  trendVars <- intersect(trendVars, names(dat))
  trendVars <- trendVars[vapply(trendVars, function(v) {
    x <- dat[[v]]
    if (is.null(x)) return(FALSE)
    if (all(is.na(x))) return(FALSE)
    dplyr::n_distinct(x, na.rm = TRUE) > 1
  }, logical(1))]
  
  # Get covariates as factors
  dat <- dat %>%
    dplyr::mutate(across(where(is.logical), ~ factor(.x, levels = c(FALSE, TRUE))))
  
  # Declustering weights
  dat$w_clust <- 1
  if (!is.null(dclust)) {
    dat$w_clust <- radius_reg_weights(dat, r_km = dclust, kde = kde, lon = lon_col, lat = lat_col)
  }
  
  ########### MODELS ##############
  f_lin_z <- reformulate(trendVars, "GWSlp")
  m_lin_z <- stats::lm(f_lin_z, data = dat, weights = w_clust)
  se_hc3_z <- sqrt(diag(get_vcov_safe(m_lin_z)))
  beta_z_lin <- stats::coef(m_lin_z)["(Intercept)"]
  se_z_lin <- se_hc3_z["(Intercept)"]
  t_z_lin <- beta_z_lin / se_z_lin
  df_z <- stats::nobs(m_lin_z) - length(stats::coef(m_lin_z))
  p_z_lin <- 2 * stats::pt(-abs(t_z_lin), df_z)
  
  f_lin <- reformulate(c("dist_LB_km", trendVars), "GWSlp")
  m_lin <- stats::lm(f_lin, data = dat, weights = w_clust)
  se_hc3 <- sqrt(diag(get_vcov_safe(m_lin)))
  beta_lin <- stats::coef(m_lin)["dist_LB_km"]
  se_lin <- se_hc3["dist_LB_km"]
  t_lin <- beta_lin / se_lin
  df <- stats::nobs(m_lin) - length(stats::coef(m_lin))
  p_lin <- 2 * stats::pt(-abs(t_lin), df)
  
  n_unique_dist <- length(unique(stats::na.omit(dat$dist_LB_km)))
  if (n_unique_dist <= 3) {
    p_lin <- NA_real_
    shape <- "degenerate (too few x points)"
    direction <- NA_character_
  }
  
  # Plot (optional)
  if (plt && nrow(dat) > nM) {
    dir.create(outdir, showWarnings = FALSE)
    dat$y <- stats::model.response(stats::model.frame(m_lin))
    grid <- data.frame(
      dist_LB_km = seq(min(dat$dist_LB_km, na.rm = TRUE),
                       max(dat$dist_LB_km, na.rm = TRUE),
                       length.out = 100)
    )
    grid$pred_lin <- predict(m_lin, newdata = grid)
    
    p <- ggplot2::ggplot(dat, ggplot2::aes(dist_LB_km, y)) +
      ggplot2::geom_point(alpha = 0.3, size = 1) +
      ggplot2::geom_line(data = grid, ggplot2::aes(y = pred_lin),
                         color = "blue",
                         linetype = ifelse(p_lin < alpha, "solid", "dashed")) +
      ggplot2::labs(
        x = "Distance to border (km)", y = "GWSlp",
        title = paste0("Aquifer ", Aquifer, ", TB: ", unique(dat$TB)),
        subtitle = sprintf("b=%.2f | p=%.3f | n=%d%s",
                           beta_lin, p_lin, nrow(dat),
                           ifelse(is.null(dclust), "", paste0(" | dclust=", dclust, "km", ifelse(kde, " (kde)", ""))))
      ) +
      ggplot2::theme_classic(base_size = 10)
    
    ggplot2::ggsave(
      file.path(outdir, paste0("aquifer_", Aquifer, "_TB", unique(dat$TB),
                               ifelse(is.null(dclust), "", paste0("_dclust", dclust, ifelse(kde, "_kde", ""))),
                               ".pdf")),
      p, width = 6, height = 4, dpi = 150
    )
  }
  
  dplyr::tibble(
    n = nrow(dat),
    beta_0 = beta_z_lin, se_0 = se_z_lin, p_0 = p_z_lin,
    beta_1 = beta_lin,   se_1 = se_lin,   p_1 = p_lin
  )
}