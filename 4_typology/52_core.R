# =============================================================================
# 52_core.R -- CORE LIBRARY (transboundary aquifer dyad typology)
#
# One file, no side effects: sourcing it (after 51_config.R) defines every
# function the run scripts need. It classifies nothing and writes nothing.
#
# What this replaces. The legacy folder implemented the same ordered
# classification SEVEN times -- Fig4_and_mosaic.R, Fig4_scatter.R and five
# FigSI_*.R scripts that differ from each other only in a buffer suffix, a
# threshold, or a variable name. Each copy carried its own hard-coded output
# paths, and three of them mislabelled a plotted variable in a way that
# produced an "NA" facet strip in the published SI figures. Here the rule
# exists once, in classify_dyads(), and a robustness variant is a list of
# configuration overrides -- never a copy of the rule.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(readr); library(ggplot2)
})
have_multcomp  <- requireNamespace("multcompView", quietly = TRUE)
have_alluvial  <- requireNamespace("ggalluvial",   quietly = TRUE)
have_vcd       <- requireNamespace("vcd",          quietly = TRUE)
have_ggrepel   <- requireNamespace("ggrepel",      quietly = TRUE)
options(warn = 1)

if (!exists("CFG")) stop("Source 51_config.R before 52_core.R.")

# =============================================================================
# 1. CONFIG OVERRIDES
# =============================================================================
cfg_with <- function(cfg = CFG, ...) {
  ov <- list(...)
  unknown <- setdiff(names(ov), names(CFG))
  if (length(unknown))
    stop("unknown configuration override(s): ", paste(unknown, collapse = ", "))
  for (nm in names(ov)) cfg[nm] <- list(ov[[nm]])
  cfg
}

# =============================================================================
# 2. UTILITIES AND LOGGING
# =============================================================================
`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

.log <- new.env(parent = emptyenv()); .log$path <- NULL
log_open <- function(path) {
  .log$path <- path
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(sprintf("# run log opened %s", format(Sys.time())), path)
}
say <- function(...) {
  txt <- paste0(...)
  cat(txt)
  if (!is.null(.log$path)) try(cat(txt, file = .log$path, append = TRUE),
                               silent = TRUE)
  invisible(NULL)
}
sayf <- function(...) say(sprintf(...))

# Standard error of the mean. FIX (audit 3.2): the legacy divisor was sqrt(n()),
# which counts rows whose value is NA -- three dyads have NA FGW from missing
# country totals, and aov() already drops them.
se_mean <- function(x, drop_na = CFG$se_drop_na) {
  n <- if (isTRUE(drop_na)) sum(!is.na(x)) else length(x)
  if (n < 2) return(NA_real_)
  stats::sd(x, na.rm = TRUE) / sqrt(n)
}

# =============================================================================
# 3. CACHE INTEGRITY
# =============================================================================
# Hash the object, never a printed rendering of it: str()-style digests
# truncate vectors and round numbers, so they can collide silently.
.md5_of <- function(obj) {
  f <- tempfile(); on.exit(unlink(f), add = TRUE)
  con <- file(f, "wb"); serialize(obj, con, ascii = FALSE); close(con)
  unname(tools::md5sum(f))
}
.file_stamp <- function(paths) {
  paths <- paths[!is.null(paths)]
  paths <- paths[file.exists(paths)]
  if (!length(paths)) return(character(0))
  info <- file.info(paths)
  paste0(basename(paths), ":", unname(tools::md5sum(paths)), ":", info$size)
}
.cfg_for_stamp <- function(cfg) {
  # Paths are dropped (they differ per machine and change no result); the
  # dyad_zones FLAG and the file NAMES are kept, because they do.
  drop <- c("root", "data_dir", "out_dir", "fig_dir", "cache_dir", "audit_dir",
            "log_file", "source_files", "cache_check", "gee_dir",
            "dyad_zone_dir")
  cfg[setdiff(names(cfg), drop)]
}
cache_stamp <- function(cfg = CFG, overrides = list(), extra = NULL) {
  list(cfg_md5   = .md5_of(.cfg_for_stamp(cfg)),
       ov_md5    = .md5_of(overrides[order(names(overrides))]),
       code      = .file_stamp(cfg$source_files),
       # The correction exports are analysis inputs: re-exporting them must
       # invalidate every cached variant. .file_stamp skips absent files, so
       # this is inert until they exist.
       data      = .file_stamp(c(cfg$dyad_file, dyad_zone_paths(cfg))),
       extra_md5 = .md5_of(extra),
       r_version = paste(R.version$major, R.version$minor, sep = "."),
       written_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}
stamp_matches <- function(a, b) {
  if (is.null(a) || is.null(b)) return(FALSE)
  k <- c("cfg_md5", "ov_md5", "code", "data", "extra_md5")
  all(vapply(k, function(n) identical(a[[n]], b[[n]]), logical(1)))
}
cache_read <- function(path, stamp, cfg = CFG, force = FALSE) {
  if (force || !file.exists(path)) return(NULL)
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(obj)) {
    say("  [cache] unreadable: ", basename(path), " -- recomputing\n")
    return(NULL)
  }
  if (!isTRUE(cfg$cache_check)) return(obj)
  if (!stamp_matches(attr(obj, "stamp"), stamp)) {
    say("  [cache] stale: ", basename(path), " -- recomputing\n")
    return(NULL)
  }
  obj
}
cache_write <- function(obj, path, stamp) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  attr(obj, "stamp") <- stamp
  saveRDS(obj, path)
  invisible(obj)
}

# =============================================================================
# 4. PDF DEVICE
# =============================================================================
# cairo_pdf embeds fonts (which journals require) but needs a working cairo
# build -- on macOS, XQuartz. Without it grDevices::cairo_pdf warns "failed to
# load cairo DLL" and writes NOTHING while the calling code reports success.
# capabilities("cairo") is NOT a reliable check: it returns TRUE on machines
# where the library fails to load. Probe by writing a real file, fall back to
# base pdf(), and verify every output afterwards.
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
pdf_device_note <- function() {
  if (.cairo_ok) say("  PDF device: cairo_pdf (fonts embedded).\n")
  else say("  PDF device: base pdf() -- cairo is unavailable in this R build.\n",
           "    Figures are still written, but fonts are NOT embedded; install\n",
           "    XQuartz (https://www.xquartz.org) and restart R before submission.\n")
}
save_pdf <- function(p, name, w, h, cfg = CFG) {
  path <- file.path(cfg$fig_dir, paste0(name, ".pdf"))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(path)) unlink(path)   # a failed write must not leave a stale file
  if (.cairo_ok) ggsave(path, p, width = w, height = h,
                        device = grDevices::cairo_pdf)
  else           ggsave(path, p, width = w, height = h,
                        device = grDevices::pdf, encoding = "WinAnsi.enc")
  if (!file.exists(path) || file.size(path) < 1000)
    stop("save_pdf: ", basename(path), " was not written (device failure).")
  sayf("  wrote %s.pdf\n", name)
  invisible(path)
}
# Base-graphics equivalent, for the vcd mosaics.
save_pdf_base <- function(expr, name, w, h, cfg = CFG) {
  path <- file.path(cfg$fig_dir, paste0(name, ".pdf"))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(path)) unlink(path)
  if (.cairo_ok) grDevices::cairo_pdf(path, width = w, height = h)
  else           grDevices::pdf(path, width = w, height = h,
                                 encoding = "WinAnsi.enc")
  err <- NULL
  ok <- tryCatch({ force(expr); TRUE },
                 error = function(e) { err <<- conditionMessage(e); FALSE })
  grDevices::dev.off()
  # Report the actual cause. "was not written" on its own is useless: the two
  # failure modes -- the drawing call errored, or the device produced nothing --
  # need different fixes, and the size tells them apart.
  if (!ok || !file.exists(path) || file.size(path) < 1000)
    stop("save_pdf_base: ", basename(path), " was not written.\n",
         if (!is.null(err)) paste0("  the drawing call failed: ", err, "\n")
         else if (!file.exists(path)) "  the device produced no file at all.\n"
         else sprintf("  the device produced only %d bytes.\n", file.size(path)))
  sayf("  wrote %s.pdf\n", name)
  invisible(path)
}

# =============================================================================
# 5. DATA LOADING
# =============================================================================
.cache <- new.env(parent = emptyenv())
load_dyads <- function(cfg = CFG) {
  if (!is.null(.cache$dyads)) return(.cache$dyads)
  if (!file.exists(cfg$dyad_file))
    stop("Dyad table not found: ", cfg$dyad_file,
         "\n  (1_data/FinalDiads.csv is the frozen analysis input; ",
         "build_dyads() can regenerate it if the geometry inputs are supplied.)")
  d <- read.csv(cfg$dyad_file, stringsAsFactors = FALSE) %>% as_tibble()
  need <- c("code", "name", "CC_1", "CC_2", "ArrangementScore",
            "has_river_border")
  miss <- setdiff(need, names(d))
  if (length(miss))
    stop("FinalDiads.csv is missing required column(s): ",
         paste(miss, collapse = ", "))
  # One row per dyad. The legacy scripts each applied this independently.
  n0 <- nrow(d)
  d <- d %>% distinct(code, CC_1, CC_2, .keep_all = TRUE)
  if (nrow(d) < n0)
    sayf("  dropped %d duplicate dyad row(s)\n", n0 - nrow(d))
  .cache$dyads <- d
  d
}
cache_clear <- function() rm(list = ls(.cache), envir = .cache)

# Column name for a variable / zone / side, e.g. zone_col("GW", "_B2", 1).
# Passing zone = "" gives the aquifer-country total; zone = "_C" the country
# total. Every classification and metric goes through this, so a variant that
# changes the variable or the buffer cannot leave a stale column behind.
zone_col <- function(var, zone, side) paste0(var, zone, "_", side)

require_cols <- function(d, cols, what) {
  miss <- setdiff(cols, names(d))
  if (length(miss))
    stop(what, ": column(s) not in the dyad table: ",
         paste(miss, collapse = ", "))
  invisible(TRUE)
}

# =============================================================================
# 5b. DYAD-SPECIFIC BORDER ZONES
# =============================================================================
# Every zone column in FinalDiads.csv is a GEE reduction over an aquifer-
# country SEGMENT with distance measured to `lb`, the global land-border
# layer. Two consequences, both wrong for a dyad:
#
#   - a segment facing several counterparts carries the SAME near-zone area in
#     all of them (the 3+-riparian case the audit flagged);
#   - a segment near a border that is not this dyad's border -- a third
#     country's, or a stretch of the same border running outside the aquifer
#     footprint -- carries area that this dyad should not get credit for.
#
# gee/dyad_border_corrections.js re-measures each zone against the stretch of
# the dyad's own border lying under the aquifer, one row per
# (code, CC, CC_other). apply_dyad_zones() overlays those values on the
# matching FinalDiads columns and leaves unmatched dyads at their segment
# values, so the overlay is a strict refinement that cannot introduce NA.
#
# Called from build_typology() and distance_class(), which is where cfg (and
# therefore the per-variant override) is in scope. load_dyads() is left alone
# deliberately: it caches one table for the whole session, and variants must
# be able to differ.
#
# The export arrives in PARTS -- the main script partitions by aquifer, and
# gee/finish_dyad_10k.js re-runs the ones that time out -- so each zone is
# named by a GLOB, not a filename, and every matching part is concatenated.
.dz <- new.env(parent = emptyenv())

# All parts on disk for one zone, in sorted order.
dyad_zone_parts <- function(cfg = CFG, src) {
  if (is.null(cfg$dyad_zone_files[[src]])) return(character(0))
  sort(Sys.glob(file.path(cfg$dyad_zone_dir, cfg$dyad_zone_files[[src]])))
}

dyad_zone_paths <- function(cfg = CFG) {
  if (is.null(cfg$dyad_zone_files)) return(character(0))
  unlist(lapply(names(cfg$dyad_zone_files), function(s)
    dyad_zone_parts(cfg, s)), use.names = FALSE)
}

# TRUE when every zone has at least one part AND the keys are complete. Used by
# 59_dyad_correction.R to skip cleanly, and by 51_config.R to decide whether to
# register the dyad_border variant at all.
dyad_zone_ready <- function(cfg = CFG) {
  if (is.null(cfg$dyad_zone_files)) return(FALSE)
  all(vapply(names(cfg$dyad_zone_files),
             function(s) length(dyad_zone_parts(cfg, s)) > 0L, logical(1)))
}

# Read and concatenate one zone's parts. Repeated keys arise where the finisher
# re-exported an aquifer the main export had already completed. They agree only
# to within floating-point noise -- Earth Engine re-serialises the reduction, so
# the last one or two ulps move -- so the collapse is TOLERANCE-based and takes
# one row rather than summing, which would silently double those keys.
read_dyad_zone <- function(cfg, src) {
  parts <- dyad_zone_parts(cfg, src)
  if (!length(parts))
    stop("no file matching ", cfg$dyad_zone_files[[src]], " in ",
         cfg$dyad_zone_dir,
         "\n  Run 4_typology/gee/dyad_border_corrections.js in the Earth ",
         "Engine Code Editor and put the CSVs in 1_data/geeOut/.")
  need <- c("code", "CC", "CC_other")
  vars <- cfg$dyad_zone_vars
  z <- bind_rows(lapply(parts, function(p) {
    x <- read.csv(p, stringsAsFactors = FALSE)
    miss <- setdiff(need, names(x))
    if (length(miss))
      stop(basename(p), " is missing column(s): ",
           paste(miss, collapse = ", "),
           "\n  Expected one row per (code, CC, CC_other).")
    v <- intersect(vars, names(x))
    if (!length(v))
      stop(basename(p), " carries none of the expected value columns (",
           paste(vars, collapse = ", "), ").")
    x[, c(need, v), drop = FALSE]
  }))
  v <- intersect(vars, names(z))
  rtol <- cfg$dyad_zone_dup_tol %||% 1e-6
  agg <- z %>%
    group_by(across(all_of(need))) %>%
    summarise(across(all_of(v),
                     list(lo = ~ min(.x, na.rm = TRUE),
                          hi = ~ max(.x, na.rm = TRUE)),
                     .names = "{.col}..{.fn}"), .groups = "drop")
  for (nm in v) {
    lo <- agg[[paste0(nm, "..lo")]]; hi <- agg[[paste0(nm, "..hi")]]
    den <- pmax(abs(hi), .Machine$double.eps)
    bad <- which(is.finite(lo) & is.finite(hi) & (hi - lo) / den > rtol)
    if (length(bad))
      stop("dyad zone ", src, ": duplicate key(s) DISAGREE on ", nm, " -- ",
           sprintf("%s %s-%s: %.10g vs %.10g", agg$code[bad[1]],
                   agg$CC[bad[1]], agg$CC_other[bad[1]], lo[bad[1]],
                   hi[bad[1]]),
           "\n  (", length(bad), " key(s)). Two exports of the same aquifer ",
           "returned different areas; re-export rather than choosing one.")
  }
  out <- agg[, need, drop = FALSE]
  for (nm in v) out[[nm]] <- agg[[paste0(nm, "..lo")]]
  out
}

apply_dyad_zones <- function(d, cfg = CFG, quiet = TRUE) {
  if (!isTRUE(cfg$dyad_zones)) return(d)

  paths <- dyad_zone_paths(cfg)
  if (!dyad_zone_ready(cfg))
    stop("dyad_zones = TRUE but the correction export(s) are missing:\n  ",
         paste(unlist(cfg$dyad_zone_files), collapse = "\n  "),
         "\n  Run 4_typology/gee/dyad_border_corrections.js in the Earth ",
         "Engine Code Editor and put the CSVs in 1_data/geeOut/.")

  key <- .md5_of(list(.file_stamp(paths), cfg$dyad_zone_files,
                      cfg$dyad_zone_vars, dim(d), names(d)))
  hit <- .dz[[key]]
  if (!is.null(hit)) return(hit)

  vars <- cfg$dyad_zone_vars
  tol  <- cfg$dyad_zone_tol %||% 1e-6
  n_repl <- 0L; viol <- list(); covered <- rep(FALSE, nrow(d))

  for (src in names(cfg$dyad_zone_files)) {
    z <- read_dyad_zone(cfg, src)
    need <- c("code", "CC", "CC_other")
    v <- intersect(vars, names(z))

    for (side in 1:2) {
      other <- if (side == 1L) 2L else 1L
      zz <- z %>% rename_with(~ paste0(.x, "..dz"), all_of(v))
      # (code, this side's country, the OTHER side's country). The export
      # carries both directions, so both sides find their own row.
      by <- c("code",
              setNames("CC",       paste0("CC_", side)),
              setNames("CC_other", paste0("CC_", other)))
      d <- left_join(d, zz, by = by)
      for (nm in v) {
        oldc <- zone_col(nm, src, side)
        newc <- paste0(nm, "..dz")
        if (!oldc %in% names(d)) { d[[newc]] <- NULL; next }
        h <- !is.na(d[[newc]])
        bad <- which(h & d[[newc]] > d[[oldc]] * (1 + tol) + 1e-12)
        if (length(bad))
          viol[[length(viol) + 1L]] <- data.frame(
            column = oldc, code = d$code[bad],
            CC_1 = d$CC_1[bad], CC_2 = d$CC_2[bad],
            segment = d[[oldc]][bad], dyad = d[[newc]][bad],
            stringsAsFactors = FALSE)
        d[[oldc]][h] <- d[[newc]][h]
        d[[newc]] <- NULL
        n_repl <- n_repl + sum(h)
        covered <- covered | h
      }
    }
  }

  # Corrected values that EXCEED the segment value they replace. These are not,
  # as first supposed, impossible. The two zones are not nested: `lb` is
  # ms_innerlines() of the country layer and exists only where two countries
  # share a LAND boundary, whereas the dyad zone is the counterpart segment
  # buffered, which can reach across water. AS140 QAT-SAU across the Gulf of
  # Salwa and AF081 AGO-COD at the Congo mouth are the clear cases, and four
  # values replace an exact zero. Below that sits a construction artefact of a
  # few per cent -- a raster distance field against a geodesic vector buffer
  # built from a counterpart simplified to 500 m - 2 km.
  #
  # So this is REPORTED, with the full list written out, and stops the run only
  # when the excess is large enough to indicate a join error rather than either
  # real effect. Adjust cfg$dyad_zone_excess_stop, not this comment, if the
  # export changes.
  if (length(viol)) {
    vv <- do.call(rbind, viol)
    vv$rel_excess <- ifelse(vv$segment > 0, vv$dyad / vv$segment - 1, Inf)
    vv <- vv[order(-vv$rel_excess), ]
    if (!is.null(cfg$out_dir) && dir.exists(cfg$out_dir))
      utils::write.csv(vv, file.path(cfg$out_dir,
                                     "table_dyad_border_excess.csv"),
                       row.names = FALSE)
    fin <- vv$rel_excess[is.finite(vv$rel_excess)]
    msg <- sprintf(
      paste0("%d corrected value(s) on %d dyad(s) exceed the segment value ",
             "they replace\n  (median %+.1f%%, %d at or above 10%%, %d ",
             "replacing an exact zero).\n  Expected: the dyad zone is the ",
             "counterpart segment buffered and can cross water, while `lb` ",
             "exists\n  only on shared LAND boundaries. See ",
             "table_dyad_border_excess.csv.\n"),
      nrow(vv), nrow(dplyr::distinct(vv, code, CC_1, CC_2)),
      100 * (if (length(fin)) stats::median(fin) else NA_real_),
      sum(fin >= 0.10), sum(vv$segment == 0))
    stop_at <- cfg$dyad_zone_excess_stop %||% 5
    if (isTRUE(cfg$dyad_zone_strict) || max(fin, 0) > stop_at)
      stop(msg, "  Largest relative excess is ", sprintf("%+.0f%%", 100 *
           max(fin, 0)), ", above the ", sprintf("%+.0f%%", 100 * stop_at),
           " tolerance: check the join keys.")
    if (!quiet) say(msg)
  }

  if (!quiet)
    sayf(paste0("  dyad zones: %d value(s) replaced across %d of %d dyads; ",
                "%d dyad(s) kept segment values.\n"),
         n_repl, sum(covered), nrow(d), sum(!covered))

  attr(d, "dyad_zone_covered") <- covered
  .dz[[key]] <- d
  d
}

# =============================================================================
# 6. THE CLASSIFICATION -- ONE IMPLEMENTATION
# =============================================================================
# Deliberately written in base vectorised R rather than case_when(): the rule
# is the scientific content of this pipeline and should be readable and
# testable without a data-frame grammar in the way. Order matters and is the
# published order; each branch is evaluated only on rows not already assigned.
classify_dyads <- function(d, cfg = CFG, future = FALSE) {
  v <- cfg$land_var; near <- cfg$near_zone; far <- cfg$far_zone

  cols <- c(zone_col(v, "", 1:2), zone_col(v, near, 1:2), zone_col(v, far, 1:2),
            zone_col("UR", far, 1:2))
  if (future)
    cols <- c(cols, zone_col(cfg$future_var, "", 1:2),
              zone_col(cfg$future_var, near, 1:2),
              zone_col(cfg$future_var, far, 1:2))
  require_cols(d, cols, "classify_dyads")

  g  <- lapply(1:2, function(s) d[[zone_col(v, "",   s)]])
  bn <- lapply(1:2, function(s) d[[zone_col(v, near, s)]])
  bf <- lapply(1:2, function(s) d[[zone_col(v, far,  s)]])
  ur <- lapply(1:2, function(s) d[[zone_col("UR", far, s)]])

  if (future) {
    # Per zone, the maximum of current extent and water-stressed cropland.
    fv <- cfg$future_var
    for (s in 1:2) {
      g[[s]]  <- pmax(g[[s]],  d[[zone_col(fv, "",   s)]])
      bn[[s]] <- pmax(bn[[s]], d[[zone_col(fv, near, s)]])
      bf[[s]] <- pmax(bf[[s]], d[[zone_col(fv, far,  s)]])
    }
  }

  eps <- cfg$eps; thr <- cfg$thresh
  # Absolute floor on the side that triggers BL. One unit of the cropland bands
  # is 1e8 m2 = 10,000 ha (see the unit note in 51_config.R), so bl_min_ha =
  # 1e-2 is 100 ha.
  # Zero by default, so the preferred specification is unchanged; the bl_floor
  # variant sets it. Rationale: BL is a SHARE, so it fires on a side holding a
  # trivial area. Of the 45 BL dyads in the published classification, 24 bind
  # on a side holding under 100 ha and 16 on under 50 ha; 34 bind on under
  # 1,000 ha. The floor asks whether the category survives when the
  # concentrated side has to hold a non-trivial amount of cropland. At 100 ha
  # BL falls to 21, at 500 ha to 16 and at 5,000 ha to 3.
  bl_floor <- cfg$bl_min_ha %||% 0
  # Share of a side's total that sits inside the near / far zone.
  r_near <- pmax(bn[[1]] / g[[1]], bn[[2]] / g[[2]])
  r_far  <- pmax(bf[[1]] / g[[1]], bf[[2]] / g[[2]])

  cls <- rep(NA_character_, nrow(d))
  assign_if <- function(cls, test, label) {
    hit <- is.na(cls) & !is.na(test) & test
    cls[hit] <- label
    cls
  }
  cls <- assign_if(cls, pmax(bf[[1]], bf[[2]]) < eps &
                        pmax(ur[[1]], ur[[2]]) < eps, "IU0")
  cls <- assign_if(cls, (bf[[1]] < eps & ur[[1]] < eps) |
                        (bf[[2]] < eps & ur[[2]] < eps), "IU1")
  cls <- assign_if(cls, pmin(bf[[1]], bf[[2]]) < eps, "UB")
  # BL requires a side that is BOTH concentrated above `thr` AND, if a floor is
  # set, larger than it. With bl_floor = 0 this is identical to r_near > thr.
  bl_hit <- ((bn[[1]] / g[[1]] > thr) & (g[[1]] >= bl_floor)) |
            ((bn[[2]] / g[[2]] > thr) & (g[[2]] >= bl_floor))
  cls <- assign_if(cls, bl_hit, "BL")
  # Everything that reaches here is bilateral and distributed. The legacy code
  # had a redundant `pmax(rG2) > thresh ~ "DA"` branch before an unreachable
  # `TRUE ~ "DA"` fallback; both resolve to DA, so this is the same rule.
  cls[is.na(cls)] <- "DA"

  # Manual reclassifications from the external audit (Section S5). NULL in the
  # preferred specification, so this changes no published number. Applied to the
  # CURRENT classification only: the audit tested present-day land use, and the
  # future scenario is a projection that no external source can corroborate.
  # Every entry states the class it expects to replace, and the run stops if the
  # rule no longer produces that class -- so the override cannot silently
  # survive a change upstream of it.
  # Applied to BOTH the current and the future classification. The audit tested
  # present-day land use, but `classF` takes pmax(current, future) per zone, so a
  # correction establishing presence on a side necessarily also establishes it in
  # the future scenario -- the propagation is monotone by construction. Each row
  # therefore carries two pairs: from/to for `class` and fromF/toF for `classF`.
  # Where the future scenario already reflected the presence (because the
  # water-stressed cropland proxy dominates that zone) toF equals fromF and the
  # row is a guarded no-op. A table without fromF/toF is treated as current-only,
  # which is the pre-audit behaviour.
  mc <- cfg$manual_class
  if (!is.null(mc) && nrow(mc)) {
    has_f <- all(c("fromF", "toF") %in% names(mc))
    if (!future || has_f) {
      key  <- paste(d$code, d$CC_1, d$CC_2)
      mkey <- paste(mc$code, mc$CC_1, mc$CC_2)
      hit  <- match(mkey, key)
      if (anyNA(hit))
        stop("manual_class: dyad not found in the dyad table: ",
             paste(mkey[is.na(hit)], collapse = "; "))
      fr <- if (future) mc$fromF else mc$from
      to <- if (future) mc$toF   else mc$to
      bad <- cls[hit] != fr
      if (any(bad))
        stop("manual_class: the rule no longer produces the expected ",
             if (future) "future class" else "class", " for ",
             paste(mkey[bad], collapse = "; "),
             ". Re-check the audit entry before overriding.")
      cls[hit] <- to
    }
  }

  factor(cls, levels = cfg$class_levels)
}

# Current and future class, plus the outcome metrics, in one table.
build_typology <- function(cfg = CFG, d = NULL) {
  if (is.null(d)) d <- load_dyads(cfg)
  # Overlay dyad-specific zones before anything reads a zone column. A no-op
  # unless cfg$dyad_zones is TRUE.
  d <- apply_dyad_zones(d, cfg)
  v <- cfg$land_var
  require_cols(d, c(zone_col(v, "", 1:2), zone_col(v, "_C", 1:2),
                    zone_col(paste0(v, "3"), "", 1:2)), "build_typology")

  tot <- d[[zone_col(v, "", 1)]] + d[[zone_col(v, "", 2)]]
  d %>%
    mutate(
      class  = classify_dyads(d, cfg, future = FALSE),
      classF = classify_dyads(d, cfg, future = TRUE),
      # FGW: the share of the two countries' COMBINED national extent that
      # overlies the dyad. NOTE (audit 2.1): this pooled ratio is what the
      # published figure plots; the Methods sentence describing "the larger of
      # the two values" is being corrected to match. The pooled version also
      # yields 3 NAs against 14 for the per-country maximum.
      FGW = tot / (d[[zone_col(v, "_C", 1)]] + d[[zone_col(v, "_C", 2)]]),
      # FBW: share of the dyad's extent under >=3-month water scarcity.
      FBW = (d[[zone_col(paste0(v, "3"), "", 1)]] +
             d[[zone_col(paste0(v, "3"), "", 2)]]) / tot,
      land_total = tot)
}

# Distance class from aquifer EXTENT within each buffer (audit 2.4b): whether
# inland exploitation is geometrically feasible, not where use is observed.
distance_class <- function(d, cfg = CFG) {
  # Area_B1 / Area_B2 are zone columns too, so the distance class has to be
  # computed from the same corrected geometry as the classification.
  d <- apply_dyad_zones(d, cfg)
  near <- cfg$near_zone; far <- cfg$far_zone
  require_cols(d, c(zone_col("Area", "", 1:2), zone_col("Area", near, 1:2),
                    zone_col("Area", far, 1:2)), "distance_class")
  r_near <- pmax(d[[zone_col("Area", near, 1)]] / d[[zone_col("Area", "", 1)]],
                 d[[zone_col("Area", near, 2)]] / d[[zone_col("Area", "", 2)]])
  r_far  <- pmax(d[[zone_col("Area", far, 1)]] / d[[zone_col("Area", "", 1)]],
                 d[[zone_col("Area", far, 2)]] / d[[zone_col("Area", "", 2)]])
  lv <- cfg$dist_class_levels
  factor(ifelse(r_near >= cfg$thresh, lv[1],
         ifelse(r_far  >= cfg$thresh, lv[2], lv[3])), levels = lv)
}

# =============================================================================
# 6b. THE DYAD-BORDER CORRECTION, AS A LEDGER
# =============================================================================
# Runs the same specification twice -- once on the segment zones, once on the
# dyad zones -- and returns everything the SI needs to describe the difference:
# a per-dyad ledger, the two correction matrices, the counts, and which dyads
# the audit also touches. Reads nothing that build_typology() does not, and
# writes nothing; 59_dyad_correction.R does the writing.
dyad_border_correction <- function(cfg = CFG, audit = NULL) {
  if (!dyad_zone_ready(cfg))
    stop("dyad_border_correction(): the correction exports are not present.")
  seg <- cfg; seg$dyad_zones <- FALSE; seg$manual_class <- NULL
  dyd <- cfg; dyd$dyad_zones <- TRUE;  dyd$manual_class <- NULL

  d  <- load_dyads(seg)
  a  <- build_typology(seg, d)
  dc <- apply_dyad_zones(d, dyd, quiet = FALSE)
  b  <- build_typology(dyd, d)
  covered <- attr(dc, "dyad_zone_covered")
  if (is.null(covered)) covered <- rep(NA, nrow(d))

  v <- seg$land_var
  share <- function(x, cf) {
    n <- pmax(x[[zone_col(v, cf$near_zone, 1)]] / x[[zone_col(v, "", 1)]],
              x[[zone_col(v, cf$near_zone, 2)]] / x[[zone_col(v, "", 2)]])
    ifelse(is.finite(n), n, NA_real_)
  }

  led <- tibble(
    code = d$code, CC_1 = d$CC_1, CC_2 = d$CC_2, name = d$name,
    measured   = covered,
    class_in   = as.character(a$class),  class_out  = as.character(b$class),
    classF_in  = as.character(a$classF), classF_out = as.character(b$classF),
    near_share_in  = share(d,  seg),
    near_share_out = share(dc, dyd))
  led$changed  <- led$class_in  != led$class_out
  led$changedF <- led$classF_in != led$classF_out

  # Which of these dyads the external audit also names, and whether the class
  # it recorded as the starting point still holds after this correction.
  #
  # `audit_from_stale` should now be FALSE everywhere, and that is the point:
  # AUDIT_FINAL carries `from` values re-derived ON the stage-2 frame, so this
  # column is a check on the re-derivation rather than a warning about it. A
  # TRUE here means the frozen table has drifted from the exports and the
  # manual_class guard in classify_dyads() is about to stop the run anyway.
  tag <- function(led, au, col) {
    led[[col]] <- FALSE
    if (is.null(au) || !nrow(au)) return(led)
    i <- match(paste(au$code, au$CC_1, au$CC_2),
               paste(led$code, led$CC_1, led$CC_2))
    ok <- !is.na(i)
    led[[col]][i[ok]] <- TRUE
    led
  }
  au <- audit %||% cfg$manual_class
  led <- tag(led, if (exists("AUDIT_FINAL_T1")) AUDIT_FINAL_T1 else au,
             "audit_t1")
  led <- tag(led, au, "audit_t12")
  led$audit_from <- NA_character_; led$audit_to <- NA_character_
  led$audit_from_stale <- NA
  if (!is.null(au) && nrow(au)) {
    i  <- match(paste(au$code, au$CC_1, au$CC_2),
                paste(led$code, led$CC_1, led$CC_2))
    ok <- !is.na(i)
    led$audit_from[i[ok]] <- au$from[ok]
    led$audit_to[i[ok]]   <- au$to[ok]
    led$audit_from_stale[i[ok]] <- au$from[ok] != led$class_out[i[ok]]
  }

  lv  <- cfg$class_levels
  mat <- bind_rows(
    as.data.frame(table(published = factor(led$class_in,  lv),
                        corrected = factor(led$class_out, lv)),
                  responseName = "n") %>% mutate(scenario = "current"),
    as.data.frame(table(published = factor(led$classF_in,  lv),
                        corrected = factor(led$classF_out, lv)),
                  responseName = "n") %>% mutate(scenario = "future"))

  cnt <- function(x) as.integer(table(factor(x, lv)))
  list(ledger = led, matrix = mat,
       counts = list(current_in  = cnt(led$class_in),
                     current_out = cnt(led$class_out),
                     future_in   = cnt(led$classF_in),
                     future_out  = cnt(led$classF_out)))
}

# =============================================================================
# 7. FIG. 4B/C STATISTICS
# =============================================================================
# ANOVA + Tukey HSD on the two outcome metrics, with compact-letter groupings.
#
# FIX (audit 3.1): the legacy SI scripts named the first metric FIR and then
# applied factor(Variable, c("FGW","FBW")), silently turning every FIR row into
# NA -- which is why S11-S13 carry an "NA" facet strip. Here the metric names
# are fixed and only the axis LABEL varies with cfg$land_var.
metric_labels <- function(cfg = CFG) {
  v <- if (identical(cfg$land_var, "IR")) "irrigated cropland"
       else "groundwater-irrigated cropland"
  c(FGW = paste0("National relevance\n(share of both countries' ", v, ")"),
    FBW = paste0("Overdraft exposure\n(share of dyad ", v, " under scarcity)"))
}

tukey_table <- function(tp, cfg = CFG) {
  if (!have_multcomp)
    stop("Package 'multcompView' is required for the Tukey letter groupings.")
  d_long <- tp %>%
    filter(.data$land_total > 0) %>%
    select(class, FGW, FBW) %>%
    pivot_longer(-class, names_to = "Variable", values_to = "Value")

  one <- function(dd) {
    fit <- stats::aov(Value ~ class, dd)
    tuk <- stats::TukeyHSD(fit, conf.level = cfg$tukey_conf)
    lets <- multcompView::multcompLetters4(fit, tuk,
                                           threshold = cfg$tukey_alpha)$class[[1]] %>%
      tibble::enframe("class", "letters")
    kw <- if (isTRUE(cfg$report_kruskal))
      suppressWarnings(stats::kruskal.test(Value ~ class, dd)$p.value) else NA_real_
    dd %>%
      group_by(class) %>%
      summarise(n = sum(!is.na(Value)),
                mean = mean(Value, na.rm = TRUE),
                se   = se_mean(Value, cfg$se_drop_na),
                ci   = list(tryCatch(stats::t.test(Value)$conf.int,
                                     error = function(e) c(NA_real_, NA_real_))),
                .groups = "drop") %>%
      mutate(ci_low = map_dbl(ci, 1), ci_high = map_dbl(ci, 2),
             kruskal_p = kw) %>%
      select(-ci) %>%
      left_join(lets, by = "class")
  }
  d_long %>%
    group_by(Variable) %>%
    tidyr::nest() %>%
    mutate(res = map(data, one)) %>%
    select(-data) %>%
    tidyr::unnest(res) %>%
    ungroup() %>%
    mutate(Variable = factor(Variable, levels = c("FGW", "FBW")),
           class    = factor(class, levels = cfg$class_levels))
}

# Country-clustered bootstrap of the class means (audit 4.2). A country appears
# in many dyads, so the dyads are not independent; this resamples COUNTRIES and
# recomputes the class means, giving an interval that respects that structure.
cluster_bootstrap <- function(tp, cfg = CFG) {
  if (!isTRUE(cfg$cluster_boot)) return(NULL)
  dd <- tp %>% filter(.data$land_total > 0) %>%
    select(class, FGW, FBW, CC_1, CC_2)
  ccs <- sort(unique(c(dd$CC_1, dd$CC_2)))
  idx <- lapply(ccs, function(c) which(dd$CC_1 == c | dd$CC_2 == c))
  names(idx) <- ccs
  stat <- function(z) {
    z %>% group_by(class) %>%
      summarise(FGW = mean(FGW, na.rm = TRUE),
                FBW = mean(FBW, na.rm = TRUE), .groups = "drop")
  }
  old <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()),
          add = TRUE)
  set.seed(cfg$cluster_boot_seed)
  reps <- map_dfr(seq_len(cfg$cluster_boot_reps), function(b) {
    draw <- sample(ccs, length(ccs), replace = TRUE)
    z <- dd[unlist(idx[draw], use.names = FALSE), , drop = FALSE]
    if (!nrow(z)) return(NULL)
    stat(z) %>% mutate(rep = b)
  })
  reps %>%
    pivot_longer(c(FGW, FBW), names_to = "Variable", values_to = "Value") %>%
    group_by(Variable, class) %>%
    summarise(boot_lo = quantile(Value, 0.025, na.rm = TRUE),
              boot_hi = quantile(Value, 0.975, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(Variable = factor(Variable, levels = c("FGW", "FBW")))
}

# Contingency test for the mosaics. FIX (audit 4.3): with 11 of 20 cells below
# an expected count of 5, the asymptotic chi-square is unreliable; a Monte
# Carlo p-value is reported instead and the cell sparsity is recorded.
contingency_test <- function(tab, cfg = CFG) {
  ex <- suppressWarnings(stats::chisq.test(tab)$expected)
  small <- sum(ex < 5); cells <- length(ex)
  tt <- if (isTRUE(cfg$chisq_simulate))
    suppressWarnings(stats::chisq.test(tab, simulate.p.value = TRUE,
                                       B = cfg$chisq_B))
  else suppressWarnings(stats::chisq.test(tab))
  tibble(statistic = unname(tt$statistic), p_value = tt$p.value,
         method = tt$method,
         cells_expected_lt5 = small, n_cells = cells,
         note = if (small > 0)
           sprintf("%d of %d cells have expected count < 5", small, cells)
         else "no sparse cells")
}

# =============================================================================
# 8. REBUILDING THE DYAD TABLE FROM GEOMETRY (optional)
# =============================================================================
# 1_data/FinalDiads.csv is the frozen analysis input and nothing in the
# pipeline calls this function. It exists so the provenance chain is not lost:
# it is a port of the legacy 1_makeDyad.R, with the paths resolved from CFG
# rather than hard-coded (the legacy script pointed at "0_PreporcessIGRAC",
# misspelled, and at a RawData/ tree that is not in the repository).
#
# Two of the three geometry inputs are not currently present anywhere, so this
# stops with a message naming what is missing rather than failing part-way.
build_dyads <- function(cfg = CFG, out_file = NULL) {
  need <- c(tba_shapefile = cfg$tba_shapefile,
            rivers_shapefile = cfg$rivers_shapefile,
            agreements_file = cfg$agreements_file)
  missing <- names(need)[vapply(need, function(p)
    is.null(p) || !nzchar(p) || !file.exists(p), logical(1))]
  if (length(missing))
    stop("build_dyads() needs geometry inputs that are not configured or not ",
         "found:\n  ", paste(sprintf("%s = %s", missing,
                                     vapply(need[missing], function(p)
                                       if (is.null(p)) "NULL" else p,
                                       character(1))), collapse = "\n  "),
         "\nSet these in 51_config.R. FinalDiads.csv ships as the frozen ",
         "analysis input, so this is only needed to regenerate it.")
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required to rebuild the dyad table.")
  out_file <- out_file %||% cfg$dyad_file
  say("Rebuilding the dyad table from geometry. This is slow.\n")

  tba <- sf::st_read(cfg$tba_shapefile, quiet = TRUE)
  codes <- unique(tba$code)

  # A dyad exists where two country polygons of the same aquifer share a
  # boundary of non-zero length.
  pairs_in_aquifer <- function(cd) {
    g <- tba[tba$code == cd, ]
    if (nrow(g) < 2) return(NULL)
    nb <- sf::st_touches(g)
    out <- list()
    for (i in seq_len(nrow(g))) for (j in nb[[i]]) {
      if (i >= j) next
      len <- sum(as.numeric(sf::st_length(suppressWarnings(sf::st_intersection(
        sf::st_boundary(g$geometry[i]), sf::st_boundary(g$geometry[j]))))))
      if (isTRUE(len > 0))
        out[[length(out) + 1]] <- data.frame(
          code = cd, CC_1 = min(g$CC[i], g$CC[j]), CC_2 = max(g$CC[i], g$CC[j]),
          stringsAsFactors = FALSE)
    }
    if (!length(out)) NULL else distinct(bind_rows(out))
  }
  dy <- bind_rows(lapply(codes, pairs_in_aquifer)) %>%
    left_join(tba %>% select(code, name) %>% sf::st_drop_geometry() %>%
                distinct(), by = "code")
  sayf("  %d dyads across %d aquifers.\n", nrow(dy), n_distinct(dy$code))

  # Cooperation score: the lower of the two sides, then the Burchi override.
  agr <- read.csv(cfg$agreements_file, stringsAsFactors = FALSE) %>%
    select(CC, code, ArrangementScore) %>%
    filter(!is.na(CC), !is.na(code), !is.na(ArrangementScore)) %>%
    group_by(CC, code) %>%
    summarise(score = max(ArrangementScore), .groups = "drop")
  dy <- dy %>%
    left_join(agr %>% rename(CC_1 = CC, s1 = score), by = c("code", "CC_1")) %>%
    left_join(agr %>% rename(CC_2 = CC, s2 = score), by = c("code", "CC_2")) %>%
    mutate(ArrangementScore = coalesce(pmin(s1, s2), 0)) %>%
    select(-s1, -s2) %>%
    mutate(ArrangementScore = if_else(code %in% cfg$burchi_codes,
                                      as.numeric(cfg$burchi_score),
                                      as.numeric(ArrangementScore)))

  # GEE summaries, joined per side. `src` is the column suffix the rest of the
  # pipeline addresses through zone_col().
  gee <- function(f) read.csv(file.path(cfg$gee_dir, f), stringsAsFactors = FALSE)
  vars <- c("CR3", "GW", "GW3", "UR", "Area", "IR", "IR3")
  join_side <- function(x, df, side, src, by_code = TRUE) {
    keys <- if (by_code) c("code", "CC") else "CC"
    if (!by_code && "GID_0" %in% names(df)) df <- rename(df, CC = GID_0)
    v <- intersect(vars, names(df))
    z <- df %>% select(any_of(keys), any_of(v)) %>%
      group_by(across(all_of(keys))) %>%
      summarise(across(all_of(v), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
      rename_with(~ paste0(.x, src, "_", side), all_of(v))
    by <- if (by_code) c("code", setNames("CC", paste0("CC_", side)))
          else setNames("CC", paste0("CC_", side))
    left_join(x, z, by = by)
  }
  srcs <- list(list(f = "outB10k_Typol.csv",  src = "_B1",  by_code = TRUE),
               list(f = "outB100k_Typol.csv", src = "_B2",  by_code = TRUE),
               list(f = "outTBA_Typol.csv",   src = "",     by_code = TRUE),
               list(f = "outCntr_Typol.csv",  src = "_C",   by_code = FALSE),
               list(f = "outB5k_Typol.csv",   src = "_B11", by_code = TRUE),
               list(f = "outB200k_Typol.csv", src = "_B22", by_code = TRUE))
  for (s in srcs) {
    df <- gee(s$f)
    for (side in 1:2) dy <- join_side(dy, df, side, s$src, s$by_code)
  }

  # River-border flag: does the shared boundary run within tolerance of a
  # mapped river? s2 is on, so the tolerance is in metres.
  sf::sf_use_s2(TRUE)
  riv <- sf::st_read(cfg$rivers_shapefile, quiet = TRUE) %>%
    sf::st_transform(sf::st_crs(tba)) %>% sf::st_make_valid()
  flag_one <- function(cd, c1, c2) {
    g1 <- tba[tba$code == cd & tba$CC == c1, ]
    g2 <- tba[tba$code == cd & tba$CC == c2, ]
    if (!nrow(g1) || !nrow(g2)) return(FALSE)
    shared <- suppressWarnings(sf::st_intersection(
      sf::st_boundary(sf::st_union(g1$geometry)),
      sf::st_boundary(sf::st_union(g2$geometry))))
    if (!length(shared) || all(sf::st_is_empty(shared))) return(FALSE)
    any(lengths(sf::st_intersects(
      sf::st_buffer(shared, cfg$river_border_tol_m), riv)) > 0)
  }
  dy$has_river_border <- vapply(seq_len(nrow(dy)), function(i)
    isTRUE(flag_one(dy$code[i], dy$CC_1[i], dy$CC_2[i])), logical(1))

  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  write.csv(dy, out_file, row.names = FALSE)
  sayf("  wrote %s (%d dyads, %d columns)\n", basename(out_file),
       nrow(dy), ncol(dy))
  invisible(dy)
}
