library(tidyverse)
library(sf)

# ==============================================================================
# Script purpose (preamble)
# ==============================================================================
# This script builds the master analysis table by combining polygon-level summaries
# exported from Google Earth Engine (GEE) for:
#   (1) IGRAC transboundary aquifers (treatment polygons; type = "treat")
#   (2) HydroBASINS control basins (control polygons; type = "bas_cntr")
#
# It assumes you already ran the two GEE scripts and downloaded their CSV exports
# into `geeOut/`:
#   - Control basins:https://code.earthengine.google.com/4ea83dc49991c19cfdaea804b7f4de35?noload=1
#   - Aquifers (TBA): https://code.earthengine.google.com/74de9d6141ba8012ed146a1fe5615c64?noload=1
#
# What the script does
#   A) Ingest “means” tables (one row per polygon)
#      - covariatesVct*.csv : polygon means of climate, land cover, suitability, roads, etc.
#      - outVct*.csv        : polygon means of irrigation and (un)sustainability overlays
#
#   B) Ingest “distance-to-border percentile curve” tables (p0..p100 by 10)
#      - landpct*.csv       : distance percentiles for all land pixels (reference curve)
#      - irpct*.csv         : irrigated pixels
#      - gwpct*.csv         : groundwater-irrigated pixels
#      - crppct*.csv        : cropland pixels
#      - GWBWSpct*.csv      : GW irrigation under BWS≥3
#      - irBWSpct*.csv      : irrigation under BWS≥3
#      - crpCSIpct*.csv     : CSI-predicted cropland (threshold calibrated per polygon)
#      - robustness variants (e.g., *_gt0, *_gt6, *_GWS_GT0, *_GWS_GT6)
#
#   C) Compute distribution metrics (Lorenz/Gini-style) for each polygon
#      - For each mask curve (Y) vs land curve (X), compute a scalar gini
#        describing whether activity is concentrated near borders or in the interior
#        (sign convention is yours).
#
#   D) Merge everything into a single table and export
#      - Full-join covariates + outcomes + gini metrics within each type
#      - Bind treatment + controls into one dataset
#      - Add an ID column
#      - Tag polygons with “RiverBorder” using a buffered GSRB river-border layer
#      - Filter to valid treatment polygons and to the final control basin set
#      - Write `_dataMain.csv`
#
# IMPORTANT (variable naming for downstream analyses)
#   The analyses expect these gini columns to exist:
#     - giniIr_Need
#     - gini_crpGWSpct_gt0
#     - gini_gwpct_gt0
#   This script creates them explicitly (with giniIr_Need and gini_crpGWSpct_gt0
#   intentionally identical, since both refer to the same underlying curve here).
#
# Inputs
#   1) GEE exports (in `geeOut/`)
#      Treatment (aquifers):
#        covariatesVct.csv, outVct.csv,
#        landpct.csv, irpct.csv, gwpct.csv, crppct.csv,
#        GWBWSpct.csv, irBWSpct.csv, crpCSIpct.csv,
#        gwpct_gt0.csv, crppct_gt0.csv, irBWSpct_gt6.csv,
#        crpGWSpctcGT0.csv, crpGWSpctc_GWS_GT0.csv, crpGWSpct_GWS_GT6.csv
#      Controls (basins): same filenames with suffix *_B.csv
#
#   2) Spatial layers used locally in R
#      - Aquifers: ../../RawData/IGRAC/TBA_Split_2025.shp
#      - Final control basins shapefile: ../0_prepControls/matched_basins_out/controlBas2025.shp
#      - River-border layer: ../../RawData/GSRB/GSRB_Level0.shp
#
# Outputs
#   - `_dataMain.csv` : master dataset for matching + analysis
#   - `Data_S2.csv`   : treatment-only table for supplement (requires aqx in memory)
#
# CRS / geometry notes
#   - RiverBorder tagging uses EPSG:3857 + 50 km buffer.
#   - Everything else is tabular from GEE exports.
# ==============================================================================

# -------------------------------
# Settings
# -------------------------------
dir_gee <- "geeOut"

# Helper: read percentile-curve CSV and standardize names
read_pct <- function(path, id_col, cn_col) {
  read.csv(path) %>%
    as_tibble() %>%
    select({{ id_col }}, {{ cn_col }}, starts_with("p")) %>%
    rename(aq_id = {{ id_col }}, CntrName = {{ cn_col }})
}

# ==============================================================================
# 1) Covariates + mean outcomes
# ==============================================================================

# Aquifers
outCov <- read.csv(file.path(dir_gee, "covariatesVct.csv")) %>%
  as_tibble() %>%
  rename(CntrName = CntryNm) %>%
  mutate(type = "treat")

outLHS <- read.csv(file.path(dir_gee, "outVct.csv")) %>%
  as_tibble() %>%
  select(
    aq_id, CntryNm, crp, GW, Ir,
    IrNeed3, OverGW3, OverIR3,
    IrNeed0, OverGW0, OverIR0,
    IrNeed6, OverGW6, OverIR6
  ) %>%
  rename(CntrName = CntryNm) %>%
  mutate(type = "treat")

# Basins (controls)
outCovB <- read.csv(file.path(dir_gee, "covariatesVct_B.csv")) %>%
  as_tibble() %>%
  rename(aq_id = HYBAS_ID) %>%
  mutate(type = "bas_cntr")

outLHSB <- read.csv(file.path(dir_gee, "outVct_B.csv")) %>%
  as_tibble() %>%
  select(
    HYBAS_ID, CntrName, crp, GW, Ir,
    IrNeed3, OverGW3, OverIR3,
    IrNeed0, OverGW0, OverIR0,
    IrNeed6, OverGW6, OverIR6
  ) %>%
  rename(aq_id = HYBAS_ID) %>%
  mutate(type = "bas_cntr")

# ==============================================================================
# 2) Gini from two percentile curves (parallelized)
# ==============================================================================

compute_gini_from_xy_parallel <- function(x, y,
                                          save_plots = FALSE,
                                          save_dir = "plots",
                                          n_cores = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' required.")
  if (!requireNamespace("future.apply", quietly = TRUE)) stop("Package 'future.apply' required.")
  
  options(progressr.enabled = FALSE)
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  
  cols_x <- grep("^p", names(x), value = TRUE)
  cols_y <- grep("^p", names(y), value = TRUE)
  
  if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
  
  if (.Platform$OS.type == "unix") {
    plan_type <- "multicore"
    future::plan(future::multicore, workers = n_cores)
  } else {
    plan_type <- "multisession"
    future::plan(future::multisession, workers = n_cores)
  }
  message("Running on ", n_cores, " cores (", plan_type, " plan)...")
  
  compute_row <- function(row_id, x, y, cols_x, cols_y, save_plots, save_dir) {
    l <- as.numeric(x[row_id, cols_x])
    i <- as.numeric(y[row_id, cols_y])
    p <- as.numeric(sub("p", "", cols_x))
    
    if (all(is.na(l)) || all(is.na(i))) {
      return(data.frame(aq_id = x$aq_id[row_id], CntrName = x$CntrName[row_id], gini = NA_real_))
    }
    
    ord <- order(p)
    p <- p[ord]; l <- l[ord]; i <- i[ord]
    keep <- !(is.na(l) | is.na(i))
    p <- p[keep]; l <- l[keep]; i <- i[keep]
    if (length(l) < 3) {
      return(data.frame(aq_id = x$aq_id[row_id], CntrName = x$CntrName[row_id], gini = NA_real_))
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
      plot_title <- paste0("aq_id ", x$aq_id[row_id], " | ", x$CntrName[row_id],
                           " | Gini = ", round(gini, 3))
      fname <- file.path(save_dir, paste0("Lorenz_aq", x$aq_id[row_id], "_", x$CntrName[row_id], ".pdf"))
      
      g <- ggplot2::ggplot(df, ggplot2::aes(x = LandPct, y = IrrPct)) +
        ggplot2::geom_line(linewidth = 1, color = "#0072B2") +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
        ggplot2::coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
        ggplot2::labs(title = plot_title, x = "Land percentile", y = "Irrigation percentile") +
        ggplot2::theme_minimal()
      
      ggplot2::ggsave(fname, g, width = 3, height = 3)
    }
    
    data.frame(aq_id = x$aq_id[row_id], CntrName = x$CntrName[row_id], gini = gini)
  }
  
  results_list <- tryCatch(
    future.apply::future_lapply(
      seq_len(nrow(x)),
      compute_row,
      x = x, y = y, cols_x = cols_x, cols_y = cols_y,
      save_plots = save_plots, save_dir = save_dir,
      future.seed = TRUE,
      future.scheduling = 1
    ),
    error = function(e) {
      message("⚠ Parallel run interrupted: ", e$message)
      message("Switching to sequential mode ...")
      future::plan(future::sequential)
      lapply(
        seq_len(nrow(x)),
        compute_row,
        x = x, y = y, cols_x = cols_x, cols_y = cols_y,
        save_plots = save_plots, save_dir = save_dir
      )
    }
  )
  
  results <- bind_rows(results_list) %>%
    mutate(aq_id = as.numeric(aq_id), gini = as.numeric(gini))
  
  future::plan(future::sequential)
  message("✓ Completed (", plan_type, " mode).")
  results
}

# ==============================================================================
# 3) Percentile curves (Aquifers) and Ginis
# ==============================================================================

lpct <- read_pct(file.path(dir_gee, "landpct.csv"), aq_id, CntryNm)

Irpct     <- read_pct(file.path(dir_gee, "irpct.csv"),     aq_id, CntryNm)
Gwpct     <- read_pct(file.path(dir_gee, "gwpct.csv"),     aq_id, CntryNm)
crppct    <- read_pct(file.path(dir_gee, "crppct.csv"),    aq_id, CntryNm)
crpCSIpct <- read_pct(file.path(dir_gee, "crpCSIpct.csv"), aq_id, CntryNm)
GWBWSpct  <- read_pct(file.path(dir_gee, "GWBWSpct.csv"),  aq_id, CntryNm)
irBWSpct  <- read_pct(file.path(dir_gee, "irBWSpct.csv"),  aq_id, CntryNm)

gwpct_gt0    <- read_pct(file.path(dir_gee, "gwpct_gt0.csv"),    aq_id, CntryNm)
crppct_gt0   <- read_pct(file.path(dir_gee, "crppct_gt0.csv"),   aq_id, CntryNm)
irBWSpct_gt6 <- read_pct(file.path(dir_gee, "irBWSpct_gt6.csv"), aq_id, CntryNm)

# GWS>=3 cropland curve (available file in current exports)
crpGWSpct_gt0     <- read_pct(file.path(dir_gee, "crpGWSpctcGT0.csv"),      aq_id, CntryNm)
crpGWSpct_GWS_gt0 <- read_pct(file.path(dir_gee, "crpGWSpctc_GWS_GT0.csv"), aq_id, CntryNm)
crpGWSpct_GWS_gt6 <- read_pct(file.path(dir_gee, "crpGWSpct_GWS_GT6.csv"),  aq_id, CntryNm)

# --- Aquifers: compute Ginis (final names used in analyses) ---
giniOut <- compute_gini_from_xy_parallel(lpct, Irpct, save_plots = FALSE) %>%
  rename(giniIr = gini) %>%
  left_join(compute_gini_from_xy_parallel(lpct, Gwpct) %>% rename(giniGw = gini),
            by = c("aq_id", "CntrName")) %>%
  left_join(compute_gini_from_xy_parallel(lpct, crppct) %>% rename(giniCr = gini),
            by = c("aq_id", "CntrName")) %>%
  # Need / scarcity cropland curve:
  left_join(compute_gini_from_xy_parallel(lpct, crpGWSpct_gt0) %>% rename(giniIr_Need = gini),
            by = c("aq_id", "CntrName")) %>%
  mutate(gini_crpGWSpct_gt0 = giniIr_Need) %>%
  left_join(compute_gini_from_xy_parallel(lpct, irBWSpct) %>% rename(giniUnsustIr = gini),
            by = c("aq_id", "CntrName")) %>%
  left_join(compute_gini_from_xy_parallel(lpct, GWBWSpct) %>% rename(giniUnsustGw = gini),
            by = c("aq_id", "CntrName")) %>%
  left_join(compute_gini_from_xy_parallel(lpct, crpCSIpct) %>% rename(giniCSI = gini),
            by = c("aq_id", "CntrName")) %>%
  # GWf>0 curve gini:
  left_join(compute_gini_from_xy_parallel(lpct, gwpct_gt0) %>% rename(gini_gwpct_gt0 = gini),
            by = c("aq_id", "CntrName")) %>%
  # Optional robustness:
  left_join(compute_gini_from_xy_parallel(lpct, crppct_gt0) %>% rename(gini_crppct_gt0 = gini),
            by = c("aq_id", "CntrName")) %>%
  left_join(compute_gini_from_xy_parallel(lpct, irBWSpct_gt6) %>% rename(gini_irBWSpct_gt6 = gini),
            by = c("aq_id", "CntrName")) %>%
  left_join(compute_gini_from_xy_parallel(lpct, crpGWSpct_GWS_gt0) %>% rename(giniCrpGWS0_crpGT05 = gini),
            by = c("aq_id", "CntrName")) %>%
  left_join(compute_gini_from_xy_parallel(lpct, crpGWSpct_GWS_gt6) %>% rename(giniCrpGWS6_crpGT05 = gini),
            by = c("aq_id", "CntrName")) %>%
  mutate(type = "treat")

# ==============================================================================
# 4) Percentile curves (Basins) and Ginis
# ==============================================================================

lpctB <- read_pct(file.path(dir_gee, "landpct_B.csv"), HYBAS_ID, CntrName)

IrpctB     <- read_pct(file.path(dir_gee, "irpct_B.csv"),     HYBAS_ID, CntrName)
GwpctB     <- read_pct(file.path(dir_gee, "gwpct_B.csv"),     HYBAS_ID, CntrName)
crppctB    <- read_pct(file.path(dir_gee, "crppct_B.csv"),    HYBAS_ID, CntrName)
crpCSIpctB <- read_pct(file.path(dir_gee, "crpCSIpct_B.csv"), HYBAS_ID, CntrName)
GWBWSpctB  <- read_pct(file.path(dir_gee, "GWBWSpct_B.csv"),  HYBAS_ID, CntrName)
irBWSpctB  <- read_pct(file.path(dir_gee, "irBWSpct_B.csv"),  HYBAS_ID, CntrName)

gwpct_gt0B    <- read_pct(file.path(dir_gee, "gwpct_gt0_B.csv"),    HYBAS_ID, CntrName)
crppct_gt0B   <- read_pct(file.path(dir_gee, "crppct_gt0_B.csv"),   HYBAS_ID, CntrName)
irBWSpct_gt6B <- read_pct(file.path(dir_gee, "irBWSpct_gt6_B.csv"), HYBAS_ID, CntrName)

crpGWSpct_gt0B     <- read_pct(file.path(dir_gee, "crpGWSpctcGT0_B.csv"),      HYBAS_ID, CntrName)
crpGWSpct_GWS_gt0B <- read_pct(file.path(dir_gee, "crpGWSpctc_GWS_GT0_B.csv"), HYBAS_ID, CntrName)
crpGWSpct_GWS_gt6B <- read_pct(file.path(dir_gee, "crpGWSpct_GWS_GT6_B.csv"),  HYBAS_ID, CntrName)

# --- Basins: compute Ginis (final names used in analyses) ---
giniOutB <- compute_gini_from_xy_parallel(lpctB, IrpctB, save_plots = FALSE) %>%
  rename(giniIr = gini) %>%
  left_join(compute_gini_from_xy_parallel(lpctB, GwpctB) %>% rename(giniGw = gini),
            by = c("aq_id", "CntrName")) %>%
  left_join(compute_gini_from_xy_parallel(lpctB, crppctB) %>% rename(giniCr = gini),
            by = c("aq_id", "CntrName")) %>%
  left_join(compute_gini_from_xy_parallel(lpctB, crpGWSpct_gt0B) %>% rename(giniIr_Need = gini),
            by = c("aq_id", "CntrName")) %>%
  mutate(gini_crpGWSpct_gt0 = giniIr_Need) %>%
  left_join(compute_gini_from_xy_parallel(lpctB, irBWSpctB) %>% rename(giniUnsustIr = gini),
            by = c("aq_id", "CntrName")) %>%
  left_join(compute_gini_from_xy_parallel(lpctB, GWBWSpctB) %>% rename(giniUnsustGw = gini),
            by = c("aq_id", "CntrName")) %>%
  left_join(compute_gini_from_xy_parallel(lpctB, crpCSIpctB) %>% rename(giniCSI = gini),
            by = c("aq_id", "CntrName")) %>%
  left_join(compute_gini_from_xy_parallel(lpctB, gwpct_gt0B) %>% rename(gini_gwpct_gt0 = gini),
            by = c("aq_id", "CntrName")) %>%
  left_join(compute_gini_from_xy_parallel(lpctB, crppct_gt0B) %>% rename(gini_crppct_gt0 = gini),
            by = c("aq_id", "CntrName")) %>%
  left_join(compute_gini_from_xy_parallel(lpctB, irBWSpct_gt6B) %>% rename(gini_irBWSpct_gt6 = gini),
            by = c("aq_id", "CntrName")) %>%
  left_join(compute_gini_from_xy_parallel(lpctB, crpGWSpct_GWS_gt0B) %>% rename(giniCrpGWS0_crpGT05 = gini),
            by = c("aq_id", "CntrName")) %>%
  left_join(compute_gini_from_xy_parallel(lpctB, crpGWSpct_GWS_gt6B) %>% rename(giniCrpGWS6_crpGT05 = gini),
            by = c("aq_id", "CntrName")) %>%
  mutate(type = "bas_cntr")

# ==============================================================================
# 5) Merge everything and write master dataset
# ==============================================================================

out  <- outCov  %>% full_join(outLHS,  by = c("aq_id", "CntrName", "type")) %>%
  full_join(giniOut, by = c("aq_id", "CntrName", "type"))

outB <- outCovB %>% full_join(outLHSB,  by = c("aq_id", "CntrName", "type")) %>%
  full_join(giniOutB, by = c("aq_id", "CntrName", "type"))

if (exists("outTB")) {
  outB <- outB %>% filter(!aq_id %in% outTB$aq_id)
}

out0 <- bind_rows(out, outB) %>%
  mutate(ID = row_number())

# ==============================================================================
# 6) Tag River borders
# ==============================================================================

LBr <- st_read("../../RawData/GSRB/GSRB_Level0.shp", quiet = TRUE) %>%
  st_transform(3857)

LBr_buf <- st_buffer(LBr, dist = 50000) # 50 km

aq_sf <- st_read("../../RawData/IGRAC/TBA_Split_2025.shp", quiet = TRUE) %>%
  rename(aqName = name) %>%
  st_transform(3857) %>%
  st_make_valid()

ints <- st_intersects(aq_sf, LBr_buf)
aq_idsRivs <- aq_sf$aq_id[lengths(ints) > 0] %>% unique()

ctrlVctB <- st_read("../0_prepControls/matched_basins_out/controlBas2025.shp", quiet = TRUE) %>%
  st_transform(3857) %>%
  st_make_valid()

ints <- st_intersects(ctrlVctB, LBr_buf)
aq_idsRivs <- c(aq_idsRivs, ctrlVctB$HYBAS_ID[lengths(ints) > 0] %>% unique())

out0 <- out0 %>%
  mutate(RiverBorder = !aq_id %in% aq_idsRivs)

# ==============================================================================
# 7) Filter final control set and write _dataMain.csv
# ==============================================================================

ctrl_ids <- st_read("../0_prepControls/matched_basins_out/controlBas2025.shp", quiet = TRUE)$HYBAS_ID

d <- out0 %>%
  filter(!is.na(aq_id)) %>%
  filter(type == "treat" | aq_id %in% ctrl_ids)%>%
  select(
    aq_id, CntrName, type,
    area_km2, lat_c, lon_c,
    CS_max, RDS_mainDist,
    GW, Ir, IrNeed3, IrNeed0,
    OverIR3, OverIR6,
    giniIr, giniGw,
    giniIr_Need,
    giniCSI,
    gini_crpGWSpct_gt0,
    gini_gwpct_gt0,
    RiverBorder
  )

write.csv(d, "_dataMain.csv", row.names = FALSE)

# ==============================================================================
# 8) Build dout for supplement (requires aqx in workspace)
# ==============================================================================

# This block is only used if you have an aqx object loaded (with columns:
# aq_id, code, name, CC, countrs, and geometry). If aqx is not present, skip.
if (exists("aqx")) {
  dout <- d %>%
    filter(type == "treat") %>%
    left_join(aqx %>% rename(CntrName=CntryNm)%>%select(-geometry), by = c("aq_id",'CntrName')) %>%
    rename(
      CropSuit     = CS_max,
      IrNeed       = IrNeed3,
      Overdraft    = OverIR3,
      Gw           = GW,
      giniCropSuit = giniCSI,
      giniIrNeed   = giniIr_Need
    ) %>%
    rename(
      IGRAC_Code   = code,
      AquiferName  = name,
      CountryCode  = CC,
      CountryName  = CntrName,
      Countries    = countrs,
      Irrig        = Ir,
      GWIrrig      = Gw,
      IrrigNeed    = IrNeed,
      G_Irrig      = giniIr,
      G_IrrigNeed  = giniIrNeed,
      G_IrrigGW    = giniGw,
      CSI          = CropSuit,
      G_CSI        = giniCropSuit
    )
  
  write.csv(dout, "Data_S2.csv", row.names = FALSE)
}
