# ==============================================================================
# Build analysis table (JJ split aquifers): covariates, outcomes, gini metrics
# Keys are Aquifer × country (shapeName / shapeGroup). Treatment = transboundary
# (Aquifer appears in >1 distinct country code), control otherwise.
# ==============================================================================
#https://code.earthengine.google.com/8ec4e37cf05eb8376a7416e2f098fc5d?noload=1

library(tidyverse)
library(dplyr)
library(ggplot2)
library(future)
library(future.apply)

# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------
dir_gee <- "geeOut"      # change if needed
alpha  <- 0.10           # not used here but kept for consistency

# ------------------------------------------------------------------------------
# Read “means” tables (one row per Aquifer×Country segment)
# ------------------------------------------------------------------------------
outCovB <- read.csv(file.path(dir_gee, "covariatesVct_JJ.csv")) %>%
  as_tibble() %>%
  filter(!is.na(Aquifer)) %>%
  select(-system.index, -Broader, -Details, -dummy, -.geo, -shapeType) %>%
  rename(CntrNm = shapeName, CC = shapeGroup) %>%
  mutate(
    n_countries = ave(CC, Aquifer, FUN = function(x) dplyr::n_distinct(x)),
    type = if_else(n_countries > 1, "treat", "ctrl")
  ) %>%
  filter(!CC %in% c(120, 124)) %>%
  select(-n_countries)

outLHSB <- read.csv(file.path(dir_gee, "outVct_JJ.csv")) %>%
  as_tibble() %>%
  filter(!is.na(Aquifer)) %>%
  select(-system.index, -Broader, -Details, -dummy, -.geo, -shapeType) %>%
  rename(CntrNm = shapeName, CC = shapeGroup) %>%
  filter(!CC %in% c(120, 124))

# ------------------------------------------------------------------------------
# Percentile-curve reader (keys: Aquifer + country)
# ------------------------------------------------------------------------------
read_pct_JJ <- function(path) {
  read.csv(path) %>%
    as_tibble() %>%
    select(Aquifer, shapeName, shapeGroup, starts_with("p")) %>%
    rename(CntrNme = shapeName, CC = shapeGroup) %>%
    filter(!CC %in% c(120, 124))
}

# ------------------------------------------------------------------------------
# Gini from two percentile curves (robust to missing rows via key-join)
# ------------------------------------------------------------------------------
compute_gini_from_xy_parallel <- function(x, y,
                                          keys = c("Aquifer", "CntrNme"),
                                          save_plots = FALSE,
                                          save_dir = "plots",
                                          n_cores = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' required.")
  if (!requireNamespace("future.apply", quietly = TRUE)) stop("Package 'future.apply' required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' required.")
  
  stopifnot(all(keys %in% names(x)), all(keys %in% names(y)))
  
  options(progressr.enabled = FALSE)
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  
  cols_x <- grep("^p", names(x), value = TRUE)
  cols_y <- grep("^p", names(y), value = TRUE)
  
  xy <- x %>%
    dplyr::select(dplyr::all_of(keys), dplyr::all_of(cols_x)) %>%
    dplyr::inner_join(
      y %>% dplyr::select(dplyr::all_of(keys), dplyr::all_of(cols_y)),
      by = keys,
      suffix = c("_x", "_y")
    )
  
  cols_x2 <- paste0(cols_x, "_x")
  cols_y2 <- paste0(cols_y, "_y")
  
  if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
  
  if (.Platform$OS.type == "unix") {
    plan_type <- "multicore"
    future::plan(future::multicore, workers = n_cores)
  } else {
    plan_type <- "multisession"
    future::plan(future::multisession, workers = n_cores)
  }
  message("Running on ", n_cores, " cores (", plan_type, " plan)...")
  
  compute_row <- function(row_id) {
    l <- as.numeric(xy[row_id, cols_x2])
    i <- as.numeric(xy[row_id, cols_y2])
    p <- as.numeric(sub("p", "", cols_x))
    
    aq  <- xy$Aquifer[row_id]
    ctn <- xy$CntrNme[row_id]
    
    if (all(is.na(l)) || all(is.na(i))) {
      return(data.frame(Aquifer = aq, CntrNme = ctn, gini = NA_real_))
    }
    
    ord <- order(p)
    p <- p[ord]; l <- l[ord]; i <- i[ord]
    keep <- !(is.na(l) | is.na(i))
    p <- p[keep]; l <- l[keep]; i <- i[keep]
    if (length(l) < 3) {
      return(data.frame(Aquifer = aq, CntrNme = ctn, gini = NA_real_))
    }
    
    l <- cummax(l); i <- cummax(i)
    
    IrrPct_for_Land <- approx(x = i, y = p, xout = l, rule = 2, ties = "ordered")$y
    LandPct_for_Irr <- approx(x = l, y = p, xout = i, rule = 2, ties = "ordered")$y
    
    df <- rbind(
      data.frame(LandPct = p, IrrPct = IrrPct_for_Land),
      data.frame(LandPct = LandPct_for_Irr, IrrPct = p)
    )
    
    df <- df[order(df$LandPct, df$IrrPct), ]
    df$LandPct <- pmax(pmin(df$LandPct, 100), 0)
    df$IrrPct  <- pmax(pmin(df$IrrPct,  100), 0)
    df <- df[!duplicated(round(df[, c("LandPct", "IrrPct")], 6)), ]
    
    xvals <- df$LandPct / 100
    yvals <- df$IrrPct  / 100
    area <- sum(diff(xvals) * (head(yvals, -1) + tail(yvals, -1))) / 2
    gini <- 1 - 2 * area
    
    if (save_plots) {
      plot_title <- paste0("Aquifer: ", aq, " | ", ctn, " | Gini = ", round(gini, 3))
      fname <- file.path(
        save_dir,
        paste0("Lorenz_", gsub("[^A-Za-z0-9]+", "_", aq), "_", gsub("[^A-Za-z0-9]+", "_", ctn), ".pdf")
      )
      
      g <- ggplot2::ggplot(df, ggplot2::aes(x = LandPct, y = IrrPct)) +
        ggplot2::geom_line(linewidth = 1, color = "#0072B2") +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
        ggplot2::coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
        ggplot2::labs(title = plot_title, x = "Land percentile", y = "Irrigation percentile") +
        ggplot2::theme_minimal()
      
      ggplot2::ggsave(fname, g, width = 3, height = 3)
    }
    
    data.frame(Aquifer = aq, CntrNme = ctn, gini = gini)
  }
  
  results_list <- tryCatch(
    future.apply::future_lapply(
      seq_len(nrow(xy)),
      compute_row,
      future.seed = TRUE,
      future.scheduling = 1
    ),
    error = function(e) {
      message("⚠ Parallel run interrupted: ", e$message)
      message("Switching to sequential mode ...")
      future::plan(future::sequential)
      lapply(seq_len(nrow(xy)), compute_row)
    }
  )
  
  results <- dplyr::bind_rows(results_list) %>%
    dplyr::mutate(gini = as.numeric(gini))
  
  future::plan(future::sequential)
  message("✓ Completed (", plan_type, " mode).")
  results
}

# ------------------------------------------------------------------------------
# Read percentile curves (JJ)
# ------------------------------------------------------------------------------
lpctB     <- read_pct_JJ(file.path(dir_gee, "landpct_JJ.csv"))
IrpctB    <- read_pct_JJ(file.path(dir_gee, "irpct_JJ.csv"))
GwpctB    <- read_pct_JJ(file.path(dir_gee, "gwpct_JJ.csv"))
crppctB   <- read_pct_JJ(file.path(dir_gee, "crppct_JJ.csv"))
GWBWSpctB <- read_pct_JJ(file.path(dir_gee, "GWBWSpct_JJ.csv"))
irBWSpctB <- read_pct_JJ(file.path(dir_gee, "irBWSpct_JJ.csv"))
CSIB <- read_pct_JJ(file.path(dir_gee, "crpCSIpct_JJ.csv"))


# ------------------------------------------------------------------------------
# Compute Ginis (JJ)
# ------------------------------------------------------------------------------
giniOutB <- compute_gini_from_xy_parallel(lpctB, IrpctB, save_plots = FALSE) %>%
  rename(giniIr = gini) %>%
  left_join(compute_gini_from_xy_parallel(lpctB, GwpctB) %>% rename(giniGw = gini),
            by = c("Aquifer", "CntrNme")) %>%
  left_join(compute_gini_from_xy_parallel(lpctB, crppctB) %>% rename(giniCr = gini),
            by = c("Aquifer", "CntrNme")) %>%
  left_join(compute_gini_from_xy_parallel(lpctB, irBWSpctB) %>% rename(giniUnsustIr = gini),
            by = c("Aquifer", "CntrNme")) %>%
  left_join(compute_gini_from_xy_parallel(lpctB, GWBWSpctB) %>% rename(giniUnsustGw = gini),
            by = c("Aquifer", "CntrNme"))%>%
  left_join(compute_gini_from_xy_parallel(lpctB, CSIB) %>% rename(giniCSI = gini),
            by = c("Aquifer", "CntrNme"))


# ------------------------------------------------------------------------------
# Merge and write master table (JJ)
# ------------------------------------------------------------------------------
out0 <- outCovB %>%
  full_join(outLHSB, by = c("Aquifer", "CntrNm", "CC")) %>%
  left_join(giniOutB, by = c("Aquifer" = "Aquifer", "CntrNm" = "CntrNme")) %>%
  mutate(aq_id = row_number())

# countries that have at least one treated row
treat_countries <- out0 %>%
  filter(type == "treat") %>%
  distinct(CntrNm) %>%
  pull(CntrNm)

# keep: all treat rows + only control rows in those countries
out0_keep <- out0 %>%
  filter(type == "treat" | (type != "treat" & CntrNm %in% treat_countries))



##Tag river borders

LBr <- st_read("../../../RawData/GSRB/GSRB_Level0.shp", quiet = TRUE) %>%
  st_transform(3857)

LBr_buf <- st_buffer(LBr, dist = 50000) # 50 km

aq_sf <- st_read("../../../RawData/jasechko_aquifs/jasechko_CountrySplit.shp", quiet = TRUE) %>%
  rename(aqName = Aquifer) %>%
  st_transform(3857) %>%
  st_make_valid()

ints <- st_intersects(aq_sf, LBr_buf)
aq_sf <- aq_sf %>%
  mutate(RiverBorder = lengths(ints) > 0)

out0_keep <- out0_keep %>%
  left_join(
    aq_sf %>%
      st_drop_geometry() %>%
      transmute(Aquifer = aqName, CntrNm = shapeName, RiverBorder),
    by = c("Aquifer", "CntrNm")
  )




write.csv(out0_keep, "_dataMain_JJ.csv", row.names = FALSE)
