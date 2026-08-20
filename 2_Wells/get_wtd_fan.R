# =============================================================================
# get_wtd_fan.R -- FETCH FAN ET AL. (2013) WATER TABLE DEPTH, SAMPLE AT WELLS
#
# One-shot helper for the natural water-table-depth robustness check that
# 34_run_robustness.R reports as a row of SI Table S4.
# Downloads the Fan, Li & Miguez-Macho (2013, Science 339:940, doi:10.1126/
# science.1229881) ~30 arc-second equilibrium water-table-depth grids from the
# authors' THREDDS server, samples the nearest cell at every well in
# 1_data/wellsData.csv, and writes
#
#     1_data/wtd_by_well.csv     columns: id, lon, lat, wtd_fan_m
#
# which 34_run_robustness.R picks up directly (it aggregates to segment
# means internally, because the matched design needs a segment-level
# covariate; the well-level file is kept as the canonical exchange format and
# also allows well-level use later).
#
# WHY THIS LAYER. The Fan et al. grids are a NATURAL EQUILIBRIUM simulation
# (no pumping), so the sampled depth acts as a pre-treatment proxy for natural
# lift cost -- unlike observed or transient modeled depths, which are
# post-treatment (see response to R1.4/R1.1). If you prefer the transient
# GLOBGM layer as a sensitivity, use GEE Scripts/getWTD.txt instead/in
# addition; both can coexist as separate wtd_* columns in wtd_by_well.csv.
#
# NETWORK. The server (thredds-gfnl.usc.es, U. Santiago de Compostela) was
# unreachable from the cloud sandbox on 14 Aug 2026 but is usually up from
# European networks; run this from your machine. Downloads total a few GB and
# resume across runs (files are kept under 1_data/fan_wtd/). If the server
# stays down, fall back to the GEE route (GEE Scripts/getWTD.txt).
#
# Requires: ncdf4 (install.packages("ncdf4")). Uses only base R otherwise.
# =============================================================================

suppressPackageStartupMessages({
  ok_nc <- requireNamespace("ncdf4", quietly = TRUE)
})
if (!ok_nc) stop("Package 'ncdf4' is required: install.packages('ncdf4')")

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!file.exists(file.path(.root, "31_config.R")))
  stop("Run from 2_wells/ (31_config.R not found next to this script).")

BASE     <- "http://thredds-gfnl.usc.es/thredds"
CAT_TOP  <- paste0(BASE, "/catalog/GLOBALWTDFTP/catalog.xml")
DEST_DIR <- file.path(.root, "1_data", "fan_wtd")
OUT_CSV  <- file.path(.root, "1_data", "wtd_by_well.csv")
WELLS    <- file.path(.root, "1_data", "wellsData.csv")
dir.create(DEST_DIR, showWarnings = FALSE, recursive = TRUE)
options(timeout = max(3600, getOption("timeout")))

cat("Reading catalog:", CAT_TOP, "\n")
read_url <- function(u) tryCatch(
  paste(readLines(u, warn = FALSE), collapse = "\n"), error = function(e) NULL)

# ---- discover netCDF files by traversing the THREDDS catalog ---------------
# We parse the catalog XML with regexes (attributes urlPath= and catalogRef
# href=) to avoid an xml2 dependency. Preference order: files under an
# "annualmean"-like branch, else any *.nc whose name mentions wtd/ann.
crawl <- function(cat_url, depth = 0) {
  if (depth > 3) return(character(0))
  x <- read_url(cat_url)
  if (is.null(x)) { cat("  [unreachable]", cat_url, "\n"); return(character(0)) }
  paths <- unique(regmatches(x, gregexpr('urlPath="[^"]+\\.nc"', x))[[1]])
  paths <- sub('^urlPath="', "", sub('"$', "", paths))
  refs  <- unique(regmatches(x, gregexpr('href="[^"]+catalog\\.xml"', x))[[1]])
  refs  <- sub('^href="', "", sub('"$', "", refs))
  refs  <- refs[!grepl("^http", refs)]                   # relative refs only
  base_dir <- sub("catalog\\.xml.*$", "", cat_url)
  for (r in refs) paths <- c(paths, crawl(paste0(base_dir, r), depth + 1))
  unique(paths)
}
nc_paths <- crawl(CAT_TOP)
if (!length(nc_paths))
  stop("No netCDF files discovered. The THREDDS server may be down -- ",
       "try again later or use the GEE route (GEE Scripts/getWTD.txt).")
ann <- nc_paths[grepl("annual|_ann", nc_paths, ignore.case = TRUE)]
sel <- if (length(ann)) ann else
  nc_paths[grepl("wtd|ann", basename(nc_paths), ignore.case = TRUE)]
if (!length(sel)) sel <- nc_paths
# one file per continent-ish basename; drop obvious monthly series
sel <- sel[!grepl("month", sel, ignore.case = TRUE)]
cat("Selected", length(sel), "file(s):\n"); for (p in sel) cat("   ", p, "\n")

# ---- download (resumable) ---------------------------------------------------
local_files <- character(0)
for (p in sel) {
  dst <- file.path(DEST_DIR, basename(p))
  if (file.exists(dst) && file.size(dst) > 1e6) {
    cat("  keeping existing", basename(dst), "\n")
  } else {
    u <- paste0(BASE, "/fileServer/", p)
    cat("  downloading", u, "\n")
    ok <- tryCatch(download.file(u, dst, mode = "wb", quiet = FALSE) == 0,
                   error = function(e) FALSE)
    if (!ok) { cat("  FAILED:", u, "\n"); next }
  }
  local_files <- c(local_files, dst)
}
if (!length(local_files)) stop("No files downloaded.")

# ---- wells ------------------------------------------------------------------
w <- read.csv(WELLS, stringsAsFactors = FALSE)
stopifnot(all(c("id", "lon", "lat") %in% names(w)))
w <- w[, c("id", "lon", "lat")]
w$wtd_fan_m <- NA_real_
cat("Sampling", nrow(w), "wells...\n")

# ---- nearest-cell sampling, row-slab reads ----------------------------------
sample_nc <- function(nc_file, wells) {
  nc <- ncdf4::nc_open(nc_file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  vn <- names(nc$var)
  v  <- vn[grepl("wtd|water", vn, ignore.case = TRUE)]
  v  <- if (length(v)) v[1] else vn[which.max(sapply(nc$var, function(z) z$ndims))]
  dm <- nc$var[[v]]$dim
  dn <- tolower(sapply(dm, function(z) z$name))
  ilon <- which(grepl("lon|^x$", dn))[1]; ilat <- which(grepl("lat|^y$", dn))[1]
  if (is.na(ilon) || is.na(ilat)) { cat("  [skip] no lon/lat dims in", basename(nc_file), "\n"); return(wells) }
  lon <- dm[[ilon]]$vals; lat <- dm[[ilat]]$vals
  dlon <- median(diff(lon)); dlat <- median(diff(lat))
  inb <- which(is.na(wells$wtd_fan_m) &
               wells$lon >= min(lon) - abs(dlon) & wells$lon <= max(lon) + abs(dlon) &
               wells$lat >= min(lat) - abs(dlat) & wells$lat <= max(lat) + abs(dlat))
  if (!length(inb)) return(wells)
  jx <- pmin(pmax(round((wells$lon[inb] - lon[1]) / dlon) + 1L, 1L), length(lon))
  jy <- pmin(pmax(round((wells$lat[inb] - lat[1]) / dlat) + 1L, 1L), length(lat))
  cat("  ", basename(nc_file), ":", length(inb), "wells in extent; var =", v, "\n")
  miss <- nc$var[[v]]$missval %||% NA_real_
  nd <- nc$var[[v]]$ndims
  # read one latitude row at a time (rows are small; random cell reads are not)
  for (j in sort(unique(jy))) {
    sel <- inb[jy == j]
    st <- rep(1L, nd); ct <- rep(1L, nd)
    st[ilon] <- 1L; ct[ilon] <- -1L        # full longitude row
    st[ilat] <- j;  ct[ilat] <- 1L         # one latitude
    row <- as.numeric(ncdf4::ncvar_get(nc, v, start = st, count = ct))
    val <- row[jx[jy == j]]
    if (!is.na(miss)) val[val == miss] <- NA_real_
    val[!is.finite(val) | abs(val) > 1e4] <- NA_real_
    # Fan grids store the water table as negative depth below land surface in
    # some distributions and positive depth in others; normalise to positive m.
    wells$wtd_fan_m[sel] <- abs(val)
  }
  wells
}
`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a
for (f in local_files) w <- tryCatch(sample_nc(f, w), error = function(e) {
  cat("  [error]", basename(f), conditionMessage(e), "\n"); w })

cat(sprintf("Coverage: %d of %d wells (%.1f%%) have a WTD value.\n",
            sum(is.finite(w$wtd_fan_m)), nrow(w),
            100 * mean(is.finite(w$wtd_fan_m))))
write.csv(w, OUT_CSV, row.names = FALSE)
cat("Written:", OUT_CSV, "\nNow rerun 34_run_robustness.R; the water-table row will appear.\n")
