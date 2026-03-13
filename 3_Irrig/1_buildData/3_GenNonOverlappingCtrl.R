# ==============================================================================
# Control-basin preprocessing + non-overlapping random samples (HydroSHEDS)
# ==============================================================================
# This script:
#   1) Loads candidate control basins (HydroSHEDS; must include HYBAS_ID).
#   2) Generates N random *non-overlapping* basin sets by greedy selection:
#        - randomize basin order
#        - keep a basin if it does not intersect any already kept basin
#      The procedure is repeated N times and is restartable via an .rds checkpoint.
#
# Outputs:
#   - CtrlNoOverlapHYBAS_<suffix>.rds   (checkpoint / final list of HYBAS_ID sets)
#
# Notes:
#   - Sampling uses an equal-area CRS (EPSG:6933) for robust geometry operations.
#   - "Non-overlapping" here means geometries do not intersect (touching counts as intersect).
#   - For large N (e.g., 500), expect this to be compute-heavy.
# ==============================================================================
library(terra)
library(sf)
library(progress)

#' Generate N random non-overlapping basin samples (restartable with checkpoints)
#'
#' @param ctrlVct SpatVector of HydroSHEDS basins (must contain HYBAS_ID)
#' @param N Number of random realizations to generate
#' @param crs_out Optional target CRS (default = EPSG:6933, equal-area)
#' @param checkpoint_dir Directory to store or read checkpoint files
#' @param file_suffix Optional string appended to checkpoint filenames
#' @return A list of length N, each element a vector of HYBAS_IDs kept
#'
sample_nonoverlap_sets <- function(ctrlVct, N = 10, crs_out = "EPSG:6933",
                                   checkpoint_dir = ".",
                                   file_suffix = "") {
  # --- Validate input ---
  stopifnot(inherits(ctrlVct, "SpatVector"))
  if (!"HYBAS_ID" %in% names(ctrlVct))
    stop("HYBAS_ID field not found in ctrlVct")
  
  # --- Preprocess: validity + projection ---
  ctrlVct <- makeValid(ctrlVct)
  ctrlVct <- project(ctrlVct, crs_out)
  
  # --- Construct suffix string ---
  suffix_str <- if (nzchar(file_suffix)) paste0("_", file_suffix) else ""
  
  # --- Check for existing checkpoint ---
  pattern <- sprintf("^CtrlNoOverlapHYBAS%s.rds$", suffix_str)
  chk_files <- list.files(checkpoint_dir, pattern = pattern, full.names = TRUE)
  
  if (length(chk_files) > 0) {
    last_chk <- chk_files[which.max(file.info(chk_files)$mtime)]
    message("🟡 Found existing checkpoint: ", basename(last_chk))
    samples_list <- readRDS(last_chk)
    completed <- max(which(!sapply(samples_list, is.null)))
    message("Resuming from sample ", completed + 1, " of ", N)
  } else {
    samples_list <- vector("list", N)
    completed <- 0
    message("🟢 No checkpoint found. Starting fresh.")
  }
  
  # --- Internal helper for one draw ---
  one_sample <- function(v, k, N) {
    message(sprintf("\nRunning sample %d of %d...", k, N))
    v <- v[sample(nrow(v)), ]   # randomize order
    keep_ids <- c()
    sel_geom <- v[0]
    
    pb <- progress_bar$new(
      total = nrow(v),
      format = sprintf("  Sample %d [:bar] :percent eta: :eta", k),
      clear = FALSE, width = 60
    )
    
    for (i in seq_len(nrow(v))) {
      g <- v[i]
      if (nrow(sel_geom) == 0 || !any(relate(g, sel_geom, "intersects"))) {
        keep_ids <- c(keep_ids, g$HYBAS_ID)
        sel_geom <- rbind(sel_geom, g)
      }
      pb$tick()
    }
    keep_ids
  }
  
  # --- Run remaining replicates ---
  t0 <- Sys.time()
  # Skip if already complete or over-complete
  if (completed >= N) {
    message("\n✅ Checkpoint already contains ", completed, 
            " samples (≥ N = ", N, "). Nothing to do.")
    return(samples_list)
  }
  
  # Only run forward from the next sample to N
  for (k in seq.int(from = completed + 1, to = N, by = 1)) {
    samples_list[[k]] <- one_sample(ctrlVct, k, N)
    if (k %% 10 == 0 || k == N || k == 1) {
      message("  💾 Saving checkpoint...")
      chk_file <- file.path(checkpoint_dir,
                            sprintf("CtrlNoOverlapHYBAS%s.rds", suffix_str))
      saveRDS(samples_list, chk_file)
      message("  🔹 Saved checkpoint at sample ", k)
      gc()
    }
  }
  message("\n✅ Completed ", N, " samples in ",
          round(difftime(Sys.time(), t0, units = "mins"), 2), " min")
  
  # --- Return ---
  return(samples_list)
}



ctrlVctB <- st_read("../0_PrepControls/matched_basins_out/controlBas2025.shp")

x=sample_nonoverlap_sets(vect(ctrlVctB), N = 500,checkpoint_dir = ".",file_suffix = "B") 


