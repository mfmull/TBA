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
library(sf)
# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script estimates, for each Aquifer × Country (CC) unit, the relationship
# between groundwater depletion slope (GWSlp) and distance to the nearest land
# border (dist_LB_km) using first-stage regressions on well-level data, then
# exports unit-level results under three specifications:
#   (1) baseline weighted OLS with heteroskedasticity-robust (HC3/HC1 fallback) SEs
#   (2) robust regression (MASS::rlm) as an outlier-resistant alternative
#   (3) spatial declustering sensitivity where observations are reweighted using
#       fixed-radius declustering weights across a range of radii
#
# Inputs
# - Well-level dataset: ../1_data/wellsData.csv
#     Contains GWSlp, dist_LB_km, lon/lat, Aquifer, CC, and LB_river indicator.
# - Aquifer-by-country dataset: ../1_data/aqfData.csv
#     Used to append unit-level descriptors after estimating each unit’s model.
#
# Key helper functions
# - get_vcov_safe(model):
#     Computes heteroskedasticity-robust vcov (HC3), with fallback to HC1 if HC3
#     fails or yields non-finite variances.
# - radius_reg_weights(dat, r_km, kde=FALSE, ...):
#     Computes spatial declustering weights from well coordinates in two modes:
#       * Fixed-radius neighbor count: w_i ∝ 1 / n_i(within r)
#       * Kernel density declustering (optional): w_i ∝ 1 / fhat_i (bandwidth r)
#     Weights are optionally capped at a high quantile and normalized to mean 1.
# - analyze_one(dat, trendVars=..., dclust=..., kde=...):
#     Runs two weighted OLS models per unit with robust SEs and diagnostic flags:
#       * Intercept-only (plus optional covariates): GWSlp ~ (trendVars)
#       * Border-distance model: GWSlp ~ dist_LB_km + (trendVars)
#     Returns coefficients, SEs, p-values, and flags for HC instability
#     (high leverage / singular covariance warnings).
# - analyze_one_robust(dat, trendVars=...):
#     Runs the same two-model structure using MASS::rlm, including preprocessing
#     to remove constant and aliased columns, and returns coefficient, SE, p,
#     and convergence warning flags.
#
# Main workflow
# 1) Parallel setup:
#    Uses furrr multisession with (detectCores() - 1) workers to run unit models
#    in parallel.
# 2) Load data and basic cleaning:
#    - Read wellsData.csv and aqfData.csv
#    - Convert LB_river to logical
#    - Drop selected island nations from CC (to focus on land-border contexts)
# 3) Baseline first stage (firstStageMain.csv):
#    - Split wells by Aquifer × CC
#    - For each unit, run analyze_one() (default covariates and no declustering)
#    - Append aquifer-by-country attributes (aqfData.csv)
#    - Create unit_id = "Aquifer_CC"
#    - Add unit-level mean(LB_river) as the propensity that the nearest border is a river
# 4) Robust first stage (firstStageMainRobust.csv):
#    - Same unit split, but run analyze_one_robust()
#    - Append the same unit metadata and LB_river share
# 5) Spatial declustering sensitivity (firstStageDC.csv):
#    - Loop over declustering radii dclust_seq = {0, 10, 20, …, 100} km
#    - For each radius:
#        * if 0: no declustering weights
#        * else: compute declustering weights via radius_reg_weights()
#               estimate unit-level OLS via analyze_one() and store dclust in output
#    - Append the same unit metadata and LB_river share
#
# Outputs
# - firstStageMain.csv:
#     Unit-level OLS/HC results (intercept model + dist model) + diagnostics + metadata
# - firstStageMainRobust.csv:
#     Unit-level robust regression (rlm) results + convergence flags + metadata
# - firstStageDC.csv:
#     Unit-level OLS/HC results across declustering radii, with dclust column + metadata
# ------------------------------------------------------------------------------

######Helper Functions################
get_vcov_safe=function(model) {
  out <- tryCatch(sandwich::vcovHC(model, type = "HC3"), error = function(e) NULL)
  if (is.null(out) || any(!is.finite(diag(out)))) {
    out <- sandwich::vcovHC(model, type = "HC1")
  }
  out
}
radius_reg_weights=function(dat,
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

analyze_one <- function(dat,
                        trendVars = character(0),
                        dclust = NULL,
                        kde = FALSE,
                        lon_col = "lon",
                        lat_col = "lat") {
  
  requireNamespace("dplyr")
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
  
  # Logical -> factor
  dat <- dat %>%
    dplyr::mutate(across(where(is.logical), ~ factor(.x, levels = c(FALSE, TRUE))))
  
  # Declustering weights
  dat$w_clust <- 1
  if (!is.null(dclust)) {
    dat$w_clust <- radius_reg_weights(dat, r_km = dclust, kde = kde, lon = lon_col, lat = lat_col)
  }
  
  # Warning flags we will fill
  warn_hc_hat_0 <- FALSE
  warn_hc_hat_1 <- FALSE
  warn_hc_singular_0 <- FALSE
  warn_hc_singular_1 <- FALSE
  
  # Helper: compute HC vcov, capture the "meatHC" warnings, and (optionally) flag high leverage
  get_vcov_flagged <- function(m) {
    warn_hat <- FALSE
    warn_singular <- FALSE
    
    w_handler <- function(w) {
      msg <- conditionMessage(w)
      if (grepl("hat values close to 1", msg, ignore.case = TRUE)) warn_hat <<- TRUE
      if (grepl("singular", msg, ignore.case = TRUE)) warn_singular <<- TRUE
      invokeRestart("muffleWarning")
    }
    
    V <- tryCatch(
      withCallingHandlers(get_vcov_safe(m), warning = w_handler),
      error = function(e) NULL
    )
    
    # Extra deterministic leverage flag (in case warnings do not trigger)
    # Threshold chosen to match warning spirit; adjust if you want.
    h <- tryCatch(stats::hatvalues(m), error = function(e) NULL)
    if (!is.null(h) && any(h >= 0.999999, na.rm = TRUE)) warn_hat <- TRUE
    
    list(V = V, warn_hat = warn_hat, warn_singular = warn_singular)
  }
  
  ########### MODELS ##############
  f_lin_z <- reformulate(trendVars, "GWSlp")
  m_lin_z <- tryCatch(stats::lm(f_lin_z, data = dat, weights = dat$w_clust), error = function(e) NULL)
  
  beta_z_lin <- se_z_lin <- p_z_lin <- NA_real_
  if (!is.null(m_lin_z)) {
    vc0 <- get_vcov_flagged(m_lin_z)
    warn_hc_hat_0 <- vc0$warn_hat
    warn_hc_singular_0 <- vc0$warn_singular
    
    if (!is.null(vc0$V)) {
      se_hc3_z <- sqrt(diag(vc0$V))
      if ("(Intercept)" %in% names(stats::coef(m_lin_z)) && "(Intercept)" %in% names(se_hc3_z)) {
        beta_z_lin <- stats::coef(m_lin_z)["(Intercept)"]
        se_z_lin <- se_hc3_z["(Intercept)"]
        t_z_lin <- beta_z_lin / se_z_lin
        df_z <- stats::nobs(m_lin_z) - length(stats::coef(m_lin_z))
        p_z_lin <- 2 * stats::pt(-abs(t_z_lin), df_z)
      }
    }
  }
  
  f_lin <- reformulate(c("dist_LB_km", trendVars), "GWSlp")
  m_lin <- tryCatch(stats::lm(f_lin, data = dat, weights = dat$w_clust), error = function(e) NULL)
  
  beta_lin <- se_lin <- p_lin <- NA_real_
  if (!is.null(m_lin)) {
    vc1 <- get_vcov_flagged(m_lin)
    warn_hc_hat_1 <- vc1$warn_hat
    warn_hc_singular_1 <- vc1$warn_singular
    
    if (!is.null(vc1$V)) {
      se_hc3 <- sqrt(diag(vc1$V))
      if ("dist_LB_km" %in% names(stats::coef(m_lin)) && "dist_LB_km" %in% names(se_hc3)) {
        beta_lin <- stats::coef(m_lin)["dist_LB_km"]
        se_lin <- se_hc3["dist_LB_km"]
        t_lin <- beta_lin / se_lin
        df <- stats::nobs(m_lin) - length(stats::coef(m_lin))
        p_lin <- 2 * stats::pt(-abs(t_lin), df)
      }
    }
  }
  
  # Degeneracy safeguard
  n_unique_dist <- length(unique(stats::na.omit(dat$dist_LB_km)))
  if (n_unique_dist <= 3) {
    p_lin <- NA_real_
  }
  
  dplyr::tibble(
    n = nrow(dat),
    beta_0 = beta_z_lin, se_0 = se_z_lin, p_0 = p_z_lin,
    beta_1 = beta_lin,   se_1 = se_lin,   p_1 = p_lin,
    
    # flags (HC instability)
    warn_hc_hatTooFewObsFi_0 = warn_hc_hat_0,
    warn_hc_singularCoVar_0 = warn_hc_singular_0,
    warn_hc_hatTooFewObsFit_1 = warn_hc_hat_1,
    warn_hc_singularCoVar_1 = warn_hc_singular_1
  )
}


# Robust-lm version of analyze_one:
# Adds robust (MASS::rlm) models WITH the same collinearity/constant-column checks
#     you used before for rlm, while keeping your lm + HC3 (get_vcov_safe) workflow. No spatial declustering
#
# Returns: lm (HC3) + rlm (robust) for intercept-only and with dist_LB_km

analyze_one_robust <- function(dat,
                               trendVars = character(0),
                               maxit = 50) {
  requireNamespace("dplyr")
  requireNamespace("MASS")
  
  # Drop NAs
  dat <- dat %>%
    dplyr::filter(!is.na(GWSlp), !is.na(dist_LB_km))
  
  # Keep covariate filtering logic
  trendVars <- intersect(trendVars, names(dat))
  trendVars <- trendVars[vapply(trendVars, function(v) {
    x <- dat[[v]]
    if (is.null(x)) return(FALSE)
    if (all(is.na(x))) return(FALSE)
    dplyr::n_distinct(x, na.rm = TRUE) > 1
  }, logical(1))]
  
  # Logical -> factor
  dat <- dat %>%
    dplyr::mutate(across(where(is.logical), ~ factor(.x, levels = c(FALSE, TRUE))))
  
  # --- robust-prep: remove constants + aliased cols (needed for rlm) ---
  prep_rlm <- function(f) {
    X <- stats::model.matrix(f, dat)
    df <- nrow(X) - qr(X)$rank
    
    sds <- apply(X, 2, stats::sd, na.rm = TRUE)
    const_vars <- setdiff(colnames(X)[sds == 0], "(Intercept)")
    if (length(const_vars) > 0) {
      vars <- setdiff(all.vars(f)[-1], const_vars)
      f <- stats::reformulate(vars, response = all.vars(f)[1])
      X <- stats::model.matrix(f, dat)
      df <- nrow(X) - qr(X)$rank
    }
    
    qr_obj <- qr(X)
    if (qr_obj$rank < ncol(X)) {
      keep <- colnames(X)[qr_obj$pivot[seq_len(qr_obj$rank)]]
      aliased <- setdiff(colnames(X), keep)
      aliased <- setdiff(aliased, "(Intercept)")
      if (length(aliased) > 0) {
        vars <- setdiff(all.vars(f)[-1], aliased)
        f <- stats::reformulate(vars, response = all.vars(f)[1])
        X <- stats::model.matrix(f, dat)
        df <- nrow(X) - qr(X)$rank
      }
    }
    
    list(f = f, X = X, df = df)
  }
  
  safe_rlm <- function(f) {
    pp <- prep_rlm(f)
    if (qr(pp$X)$rank < ncol(pp$X)) {
      return(list(m = NULL, df = pp$df, converged = NA, warn_converge = FALSE))
    }
    
    warn_converge <- FALSE
    w_handler <- function(w) {
      msg <- conditionMessage(w)
      if (grepl("failed to converge", msg, ignore.case = TRUE)) warn_converge <<- TRUE
      invokeRestart("muffleWarning")
    }
    
    m <- tryCatch(
      withCallingHandlers(
        MASS::rlm(pp$f, data = dat, maxit = maxit),
        warning = w_handler
      ),
      error = function(e) NULL
    )
    
    converged <- NA
    if (!is.null(m) && !is.null(m$converged)) converged <- isTRUE(m$converged)
    
    list(m = m, df = pp$df, converged = converged, warn_converge = warn_converge)
  }
  
  rlm_term <- function(rlm_obj, term, df) {
    beta <- se <- p <- NA_real_
    if (!is.null(rlm_obj)) {
      sr <- tryCatch(summary(rlm_obj)$coef, error = function(e) NULL)
      if (!is.null(sr) && term %in% rownames(sr)) {
        beta <- sr[term, "Value"]
        se   <- sr[term, "Std. Error"]
        p    <- 2 * stats::pt(-abs(sr[term, "t value"]), df)
      }
    }
    list(beta = beta, se = se, p = p)
  }
  
  # ---------------- RLM MODELS ----------------
  f_lin_z <- reformulate(trendVars, "GWSlp")
  f_lin   <- reformulate(c("dist_LB_km", trendVars), "GWSlp")
  
  r0 <- safe_rlm(f_lin_z)
  o0 <- rlm_term(r0$m, "(Intercept)", df = r0$df)
  
  r1 <- safe_rlm(f_lin)
  o1 <- rlm_term(r1$m, "dist_LB_km", df = r1$df)
  
  # Degeneracy safeguard
  n_unique_dist <- length(unique(stats::na.omit(dat$dist_LB_km)))
  if (n_unique_dist <= 3) {
    o1$p <- NA_real_
  }
  
  dplyr::tibble(
    n = nrow(dat),
    
    beta_0 = o0$beta, se_0 = o0$se, p_0 = o0$p,
    beta_1 = o1$beta, se_1 = o1$se, p_1 = o1$p,
    
    # convergence flags
    rlm0_converged = r0$converged,
    rlm1_converged = r1$converged,
    rlm0_warn_failed_to_converge = r0$warn_converge,
    rlm1_warn_failed_to_converge = r1$warn_converge
  )
}

#parallel settings  
plan(multisession, workers = parallel::detectCores() - 1)


##############################
#Fetch data
############################## 
df=read.csv('../1_data/wellsData.csv')%>%
  mutate(LB_river=(LB_river==1))%>%filter(!CC%in%c("AUS", "NZL","TWN","JPN", "IRL", "GBR"))#remove island nations
aqf=read.csv('../1_data/aqfData.csv')

###############
#Main first stage
##############
firstStage <- df %>%
  group_split(Aquifer, CC,.keep = TRUE) %>%
  future_map_dfr(~ {
    y <- dplyr::slice_head(.x, n = 1)
    out <- analyze_one(.x )
    # Add grouping info to result
    out$Aquifer <- y$Aquifer
    out$CC <- y$CC
    out
  }, .progress = TRUE)%>%
  left_join(aqf, by = c("Aquifer","CC")) %>%
  mutate(unit_id = paste0(Aquifer, "_", CC))%>%
  left_join(df%>%group_by(Aquifer,CC)%>%dplyr::summarize(LB_river=mean(LB_river))) #propensity that closest border is river.

write.csv(firstStage,'firstStageMain.csv')

#######################
#Robust first stage
######################
firstStageR <- df %>%
  group_split(Aquifer, CC,.keep = TRUE) %>%
  future_map_dfr(~ {
    y <- dplyr::slice_head(.x, n = 1)
    out <- analyze_one_robust(.x )
    # Add grouping info to result
    out$Aquifer <- y$Aquifer
    out$CC <- y$CC
    out
  }, .progress = TRUE)%>%
  left_join(aqf, by = c("Aquifer","CC")) %>%
  mutate(unit_id = paste0(Aquifer, "_", CC))%>%
  left_join(df%>%group_by(Aquifer,CC)%>%dplyr::summarize(LB_river=mean(LB_river))) #propensity that closest border is river.

write.csv(firstStageR,'firstStageMainRobust.csv')



########################
#Spatial Declust KDE
#######################
dclust_seq <- c(0, seq(10, 100, 10))

firstStageDC <- purrr::map_dfr(dclust_seq, function(dc) {
  message("Running dclust = ", dc)
  
  df %>%
    group_split(Aquifer, CC, .keep = TRUE) %>%
    future_map_dfr(~{
      y <- dplyr::slice_head(.x, n = 1)
      out <- analyze_one(
        .x,
        # kde=T,
        kde=F,
        dclust = if (dc == 0) NULL else dc
        # optionally: kde = TRUE/FALSE, lon_col = "lon", lat_col = "lat"
      )
      out$Aquifer <- y$Aquifer
      out$CC <- y$CC
      out$dclust <- dc
      out
    }, .progress = TRUE)
}) %>%
  left_join(aqf, by = c("Aquifer","CC")) %>%
  mutate(unit_id = paste0(Aquifer, "_", CC)) %>%
  left_join(df %>% group_by(Aquifer, CC) %>% summarize(LB_river = mean(LB_river), .groups="drop"),
            by = c("Aquifer","CC"))

write.csv(firstStageDC, "firstStageDC.csv", row.names = FALSE)
