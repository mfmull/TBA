# =============================================================================
# 32_core.R -- CORE LIBRARY
#
# One file, no side effects: sourcing it (after 31_config.R) defines every
# function the run scripts need. It fits nothing and writes nothing.
#
# Refactored from 21_core.R. What is preserved: the central configuration
# object; run_specification() as the single entry point (a robustness check is
# a list of overrides, never a copy of fitting code); deterministic caching
# with stamps that refuse stale objects; shared functions for the first stage,
# matching/weighting, second stage and inference; explicit provenance for
# spatial uncertainty; the H3 shrinkage machinery with BLUP validation and the
# coded selection rule. What changed: the pipeline is self-contained on
# wellsData.csv (segment covariates and all first-stage quantities are derived
# from the wells; the former firststageMain.csv is gone), the Conley kernel no
# longer needs `sf`, and figure/table writers moved to 36/37.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(metafor)
})
have_cs    <- requireNamespace("clubSandwich", quietly = TRUE)
have_repel <- requireNamespace("ggrepel",      quietly = TRUE)
have_match <- requireNamespace("MatchIt",      quietly = TRUE)
have_mass  <- requireNamespace("MASS",         quietly = TRUE)
options(warn = 1)

if (!exists("CFG")) stop("Source 31_config.R before 32_core.R.")

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
# 2. UTILITIES
# =============================================================================
`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

winsor <- function(x, p = c(0.01, 0.99)) {
  if (is.null(p)) return(x)
  qs <- quantile(x, p, na.rm = TRUE)
  pmin(pmax(x, qs[1]), qs[2])
}

wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

ess <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (!length(w)) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

crit_val <- function(df = NA_real_, level = CFG$ci_level) {
  a <- 1 - (1 - level) / 2
  if (is.finite(df) && df > 0) stats::qt(a, df) else stats::qnorm(a)
}

pstar <- function(p) {
  if (!length(p)) return("NA")
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

smd <- function(x, tb, w) {
  ok <- is.finite(x) & !is.na(tb) & is.finite(w) & w >= 0
  x <- x[ok]; tb <- tb[ok]; w <- w[ok]
  if (!any(tb == 1) || !any(tb == 0)) return(NA_real_)
  s <- stats::sd(x[tb == 1])
  if (!is.finite(s) || s == 0) return(NA_real_)
  (wmean(x[tb == 1], w[tb == 1]) - wmean(x[tb == 0], w[tb == 0])) / s
}

# Warnings collected with context.
.notes <- new.env(parent = emptyenv()); .notes$msgs <- character(0)
quiet <- function(expr, tag) {
  withCallingHandlers(expr, warning = function(w) {
    .notes$msgs <- unique(c(.notes$msgs, sprintf("[%s] %s", tag, conditionMessage(w))))
    invokeRestart("muffleWarning")
  })
}
notes_get   <- function() .notes$msgs
notes_reset <- function() .notes$msgs <- character(0)
notes_add   <- function(m) .notes$msgs <- unique(c(.notes$msgs, m))

# Run log.
.log <- new.env(parent = emptyenv()); .log$path <- NULL; .log$in_worker <- FALSE
log_open <- function(path) {
  .log$path <- path
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(sprintf("# run log opened %s", format(Sys.time())), path)
}
say <- function(...) {
  if (isTRUE(.log$in_worker)) return(invisible(NULL))
  txt <- paste0(...)
  cat(txt)
  if (!is.null(.log$path))
    try(cat(txt, file = .log$path, append = TRUE), silent = TRUE)
  invisible(NULL)
}
sayf <- function(...) say(sprintf(...))

# =============================================================================
# 3. REPRODUCIBLE PARALLEL MAP (L'Ecuyer streams; mclapply on unix)
# =============================================================================
.rng_save <- function() list(
  kind = RNGkind(),
  seed = if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv()) else NULL)
.rng_restore <- function(st) {
  suppressWarnings(RNGkind(st$kind[1], st$kind[2], st$kind[3]))
  if (!is.null(st$seed)) assign(".Random.seed", st$seed, envir = globalenv())
  invisible(NULL)
}
with_seed <- function(seed, expr, kind = "Mersenne-Twister") {
  st <- .rng_save(); on.exit(.rng_restore(st), add = TRUE)
  set.seed(seed, kind = kind)
  force(expr)
}
.rng_streams <- function(seed, n) {
  if (n < 1) return(list())
  st <- .rng_save(); on.exit(.rng_restore(st), add = TRUE)
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  s <- get(".Random.seed", envir = globalenv())
  out <- vector("list", n)
  for (i in seq_len(n)) { out[[i]] <- s; s <- parallel::nextRNGStream(s) }
  out
}
par_workers <- function(cfg = CFG) {
  w <- cfg$workers %||% max(1L, parallel::detectCores())
  min(8L, max(1L, as.integer(w)))
}
# par_map(x, f, seed): job i draws from stream i whatever process runs it, so
# sequential and parallel runs agree bit for bit and worker count is irrelevant.
# `preschedule` controls how mclapply distributes jobs, and the right answer
# depends entirely on the shape of the workload:
#
#   TRUE  -- split the jobs into one block per worker up front. One fork per
#            worker, so it is the only sane choice for very many cheap jobs
#            (the 1500+ per-unit Conley fits). But each child then runs its
#            whole block SEQUENTIALLY in one process: memory accumulates across
#            the block, there is no load balancing, and one child dying takes
#            its entire block with it -- reported as "scheduled cores ... did
#            not deliver results".
#   FALSE -- one fork per job. The child exits after a single job, so peak
#            memory is one job rather than a block; slow and fast jobs balance
#            across workers; and a failure is isolated to that job and comes
#            back as a try-error carrying a real message instead of a silent
#            block-wide loss.
#
# Few, heavy, unequal jobs (a robustness specification is a full pipeline fit:
# rma.mv, cluster-robust SEs, two 999-replicate wild bootstraps) therefore want
# FALSE. The default below picks TRUE only when the job count is large enough
# that per-job fork overhead would matter. Reproducibility is unaffected either
# way: job i draws from stream i whatever process runs it.
par_map <- function(x, f, label = "working", cfg = CFG, seed = NULL,
                    preschedule = NULL) {
  n <- length(x)
  if (n == 0L) return(list())
  if (is.null(preschedule)) preschedule <- n >= 500L
  if (!is.null(seed)) { .st <- .rng_save(); on.exit(.rng_restore(.st), add = TRUE) }
  streams <- if (is.null(seed)) vector("list", n) else .rng_streams(seed, n)
  run_one <- function(i) {
    if (!is.null(streams[[i]]))
      assign(".Random.seed", streams[[i]], envir = globalenv())
    f(x[[i]])
  }
  nw <- par_workers(cfg)
  # NEVER fork from inside a forked worker. Grandchildren inherit the child's
  # result pipe to the master; when they exit they corrupt the serialised
  # stream the master reads, and mclapply then loses EVERY job, not just the
  # nested one ("scheduled cores ... did not deliver results").
  #
  # AND NEVER FORK UNDER RStudio. Its session process is multithreaded, and on
  # macOS a fork from a process that has initialised the Objective-C runtime
  # aborts the child immediately. The children are killed rather than raising
  # an R condition, so mclapply reports missing results with no error to show:
  # measured here, eight forked children doing nothing but crossprod() on a
  # random matrix returned 0 of 8, while a single fork succeeded. R's own
  # documentation warns against mclapply in GUI front ends for this reason.
  # The parallel path is available by running the pipeline with Rscript from a
  # terminal. Escape hatch, for a session where forking has been fixed at the
  # OS level (e.g. OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES in ~/.Renviron):
  # set ALLOW_FORK_IN_RSTUDIO=1.
  no_fork_gui <- identical(Sys.getenv("RSTUDIO"), "1") &&
    !nzchar(Sys.getenv("ALLOW_FORK_IN_RSTUDIO"))
  if (no_fork_gui && isTRUE(cfg$parallel) && nw > 1L && n >= 4L &&
      !isTRUE(.log$in_worker) && !isTRUE(.log$said_no_fork)) {
    say("  Running sequentially: forking is unreliable under RStudio and the\n",
        "  workers are killed without an R error. For the parallel path run\n",
        "  the pipeline from a terminal:  Rscript 38_run_all.R\n",
        "  (override with ALLOW_FORK_IN_RSTUDIO=1 if forking works for you.)\n")
    .log$said_no_fork <- TRUE
  }
  use_par <- isTRUE(cfg$parallel) && nw > 1L && n >= 4L &&
    .Platform$OS.type == "unix" && !isTRUE(.log$in_worker) && !no_fork_gui
  # The sequential path used to be completely silent, which on a long block
  # (103 influence refits, or a 200-replicate bootstrap) is indistinguishable
  # from a hang. Report a heartbeat roughly every 10% so it is obvious the run
  # is progressing rather than stuck. Semantics are unchanged.
  if (!use_par) {
    if (n < 4L) return(lapply(seq_len(n), run_one))
    sayf("  %s: %d jobs, sequential (no forking)\n", label, n)
    .every <- max(1L, n %/% 10L)
    .t0 <- Sys.time()
    out <- vector("list", n)
    for (i in seq_len(n)) {
      out[[i]] <- run_one(i)
      if (i %% .every == 0L || i == n) {
        .el <- as.numeric(difftime(Sys.time(), .t0, units = "mins"))
        sayf("    %s: %d/%d (%.1f min elapsed, ~%.1f min left)\n",
             label, i, n, .el, .el / i * (n - i))
      }
    }
    return(out)
  }
  say(sprintf("  %s: %d jobs on %d workers (%s)\n", label, n, nw,
              if (preschedule) "prescheduled" else "one fork per job"))
  # Envelope every return value: a legitimate NULL from f() must stay
  # distinguishable from a worker that died (mclapply fills those with NULL or
  # a try-error). Otherwise a dead worker silently becomes a "result".
  wrap <- function(i) { .log$in_worker <- TRUE; list(value = run_one(i)) }
  raw <- parallel::mclapply(seq_len(n), wrap, mc.cores = nw,
                            mc.preschedule = preschedule, mc.set.seed = FALSE)
  ok <- vapply(raw, function(z)
    is.list(z) && length(z) == 1L && identical(names(z), "value"), logical(1))
  if (any(!ok)) {
    msg <- unlist(lapply(raw[!ok], function(z)
      if (inherits(z, "try-error")) conditionMessage(attr(z, "condition"))
      else NULL))
    say(sprintf("  %s: %d of %d parallel job(s) did not return; rerunning them sequentially.%s\n",
                label, sum(!ok), n,
                if (length(msg)) paste0(" First error: ", msg[1]) else ""))
    notes_add(sprintf("par_map[%s]: %d of %d jobs lost in parallel, rerun sequentially",
                      label, sum(!ok), n))
    for (i in which(!ok)) raw[[i]] <- list(value = run_one(i))
  }
  lapply(raw, `[[`, "value")
}
par_map_dfr <- function(x, f, label = "working", cfg = CFG, seed = NULL,
                        preschedule = NULL)
  bind_rows(par_map(x, f, label = label, cfg = cfg, seed = seed,
                    preschedule = preschedule))

# =============================================================================
# 4. CACHE INTEGRITY
# =============================================================================
.md5_of <- function(obj) {
  f <- tempfile(); on.exit(unlink(f), add = TRUE)
  con <- file(f, "wb"); serialize(obj, con, ascii = FALSE); close(con)
  unname(tools::md5sum(f))
}
.file_stamp <- function(paths) {
  paths <- paths[file.exists(paths)]
  if (!length(paths)) return(character(0))
  info <- file.info(paths)
  paste0(basename(paths), ":", unname(tools::md5sum(paths)), ":", info$size)
}
.cfg_for_stamp <- function(cfg) {
  drop <- c("root", "out_dir", "cache_dir", "audit_dir", "log_file",
            "source_files", "cache_check", "parallel", "workers")
  cfg[setdiff(names(cfg), drop)]
}
.canon_list <- function(x) {
  if (!is.list(x) || !length(x) || is.null(names(x))) return(x)
  x[order(names(x))]
}
cache_stamp <- function(cfg = CFG, overrides = list(), extra = NULL) {
  list(cfg_md5   = .md5_of(.cfg_for_stamp(cfg)),
       ov_md5    = .md5_of(.canon_list(overrides)),
       code      = .file_stamp(cfg$source_files),
       data      = .file_stamp(cfg$wells_file),
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
  if (is.null(obj)) return(NULL)
  if (!isTRUE(cfg$cache_check)) {
    warning("cache_check is FALSE: reusing ", basename(path), " unvalidated.")
    return(obj)
  }
  if (!stamp_matches(attr(obj, "stamp"), stamp)) {
    say("  [stale cache discarded] ", basename(path), "\n")
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
# 5. DATA LOADING (memoised)
# =============================================================================
.cache <- new.env(parent = emptyenv())
load_wells <- function(cfg = CFG) {
  if (is.null(.cache$wells)) {
    if (!file.exists(cfg$wells_file))
      stop("Input file not found: ", cfg$wells_file)
    w <- read.csv(cfg$wells_file, stringsAsFactors = FALSE) %>% as_tibble()
    need <- c("GWSlp", "dist_LB_km", "CC", "Aquifer", "TB", "lon", "lat",
              "urbkHaKm2", "CS_max", "LB_river", "prec_mm")
    miss <- setdiff(need, names(w))
    if (length(miss))
      stop("wellsData.csv is missing required column(s): ",
           paste(miss, collapse = ", "))
    .cache$wells <- w %>% mutate(unit_id = paste0(Aquifer, "_", CC))
  }
  .cache$wells
}
cache_clear <- function() rm(list = ls(.cache), envir = .cache)

# =============================================================================
# 6. FIRST STAGE FROM THE WELLS
# =============================================================================
# Per aquifer-country segment:
#   H1  beta_0  mean depletion rate (mm/yr), HC3 SE
#   H2  z       Fisher z of the within-segment depletion-distance correlation
# plus the physical slope beta_1 (table-only), a near-vs-interior binned
# contrast (table-only), segment covariates (well means) and distance geometry.

.hc3_fit <- function(X, y) {
  XtX_inv <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XtX_inv)) return(NULL)
  b <- XtX_inv %*% crossprod(X, y)
  e <- as.vector(y - X %*% b)
  h <- rowSums((X %*% XtX_inv) * X)
  A <- X * (e / (1 - pmin(h, 0.999)))
  se <- sqrt(pmax(diag(XtX_inv %*% crossprod(A) %*% XtX_inv), 0))
  list(X = X, b = as.vector(b), e = e, XtX_inv = XtX_inv, hc_se = se)
}

.beta0_of <- function(y, cfg) {
  switch(cfg$first_stage_h1,
    mean = {
      n <- length(y); e <- y - mean(y)
      list(b = mean(y), se = sqrt(sum((e / (1 - 1 / n))^2)) / n)
    },
    trimmed_mean = {
      b <- mean(y, trim = cfg$trim_frac)
      n <- length(y); e <- y - b
      list(b = b, se = sqrt(sum((e / (1 - 1 / n))^2)) / n)
    },
    robust_mean = {
      if (!have_mass) return(list(b = NA_real_, se = NA_real_))
      f <- tryCatch(MASS::rlm(y ~ 1, maxit = 50), error = function(e) NULL)
      if (is.null(f)) return(list(b = NA_real_, se = NA_real_))
      cf <- summary(f)$coefficients
      list(b = unname(cf[1, 1]), se = unname(cf[1, 2]))
    },
    stop("unknown first_stage_h1: ", cfg$first_stage_h1))
}

.first_stage_unit <- function(y, dd, cfg) {
  n <- length(y)
  out <- list(n_w = n, beta_0 = NA_real_, se_0 = NA_real_,
              beta_1 = NA_real_, se_1 = NA_real_,
              r = NA_real_, r_spear = NA_real_,
              bin_diff = NA_real_, bin_se = NA_real_,
              min_dist = NA_real_, sd_dist = NA_real_, rng_dist = NA_real_)
  if (n < 2) return(out)
  out$min_dist <- min(dd); out$sd_dist <- stats::sd(dd)
  out$rng_dist <- diff(range(dd))
  b0 <- .beta0_of(y, cfg)
  out$beta_0 <- b0$b; out$se_0 <- b0$se
  if (is.finite(out$sd_dist) && out$sd_dist > 0) {
    m1 <- .hc3_fit(cbind(1, dd), y)
    if (!is.null(m1)) { out$beta_1 <- m1$b[2]; out$se_1 <- m1$hc_se[2] }
  }
  if (n >= 10 && is.finite(stats::sd(y)) && stats::sd(y) > 0 &&
      is.finite(out$sd_dist) && out$sd_dist > 0) {
    out$r       <- stats::cor(y, dd)
    out$r_spear <- suppressWarnings(stats::cor(y, dd, method = "spearman"))
  }
  b <- cfg$distance_bins
  near <- dd <= b[1]; far <- dd > b[1] & dd <= b[2]
  if (sum(near) >= cfg$min_wells_per_bin && sum(far) >= cfg$min_wells_per_bin) {
    out$bin_diff <- mean(y[near]) - mean(y[far])
    out$bin_se <- sqrt(stats::var(y[near]) / sum(near) +
                         stats::var(y[far]) / sum(far))
  }
  out
}

# Wells entering the first stage for a given configuration.
fs_wells <- function(cfg = CFG) {
  w <- load_wells(cfg) %>%
    filter(!CC %in% cfg$drop_cc, is.finite(GWSlp), is.finite(dist_LB_km))
  if (!is.null(cfg$distance_window_km))
    w <- w %>% filter(dist_LB_km >= cfg$distance_window_km[1],
                      dist_LB_km <= cfg$distance_window_km[2])
  w
}

.fs_key <- function(cfg) paste(
  paste(cfg$drop_cc, collapse = ","),
  if (is.null(cfg$distance_window_km)) "full"
  else paste(cfg$distance_window_km, collapse = "-"),
  cfg$first_stage_h1, cfg$trim_frac,
  paste(cfg$distance_bins, collapse = "-"), cfg$min_wells_per_bin, sep = "|")

# Memoised: one HC3 pass per segment per (sample, H1-estimator) key. Also
# attaches segment covariates (well means; the aquifer-segment centroid for
# lat/lon), the TB flag, and within-distance well coverage.
first_stage_from_wells <- function(cfg = CFG) {
  key <- .fs_key(cfg)
  hit <- .cache[[paste0("fs_", key)]]
  if (!is.null(hit)) base <- hit else {
    w <- fs_wells(cfg)
    sp <- split(seq_len(nrow(w)), w$unit_id)
    est <- map_dfr(names(sp), function(u) {
      i <- sp[[u]]
      as_tibble(.first_stage_unit(w$GWSlp[i], w$dist_LB_km[i], cfg)) %>%
        mutate(unit_id = u, .before = 1)
    })
    covs <- w %>%
      group_by(unit_id) %>%
      summarise(Aquifer = first(Aquifer), CC = first(CC),
                TB = any(TB %in% c(TRUE, "TRUE", 1, "1")),
                lat_c = mean(lat, na.rm = TRUE),
                lon_c = mean(lon, na.rm = TRUE),
                urbkHaKm2 = mean(urbkHaKm2, na.rm = TRUE),
                CS_max = mean(CS_max, na.rm = TRUE),
                LB_river = mean(LB_river, na.rm = TRUE),
                prec_mm = mean(prec_mm, na.rm = TRUE),
                !!!setNames(
                  lapply(CFG$coverage_km, function(km)
                    quo(mean(dist_LB_km <= !!km, na.rm = TRUE))),
                  paste0("share_within_", CFG$coverage_km, "km")),
                .groups = "drop")
    base <- est %>% left_join(covs, by = "unit_id") %>%
      mutate(aq_id = as.character(Aquifer))
    .cache[[paste0("fs_", key)]] <- base
  }
  base %>% mutate(
    r_use = if (cfg$first_stage_h2 == "spearman") r_spear else r,
    r_use = pmin(pmax(r_use, -0.999), 0.999),
    z     = atanh(r_use),
    z_pref = atanh(pmin(pmax(r, -0.999), 0.999)))
}

# =============================================================================
# 7. FIRST-STAGE SPATIAL (CONLEY) SE INFLATION
# =============================================================================
# HC3 assumes independent wells; monitoring wells are spatially correlated, so
# HC3 understates se. Recompute with a Bartlett spatial-HAC kernel and carry
# the inflation ratio forward. Cached by (sample key, cutoff, seed); a
# restricted-sample specification gets its OWN cache rows.

project_local <- function(lon, lat) {
  R <- 6371000
  lat0 <- mean(lat, na.rm = TRUE) * pi / 180
  cbind(x = R * (lon * pi / 180) * cos(lat0),
        y = R * (lat * pi / 180))
}

# Dense Bartlett kernel (n <= n_max_wells, so at most ~2500 x 2500).
.kernel_dense <- function(xy, cutoff_m) {
  dx <- outer(xy[, 1], xy[, 1], "-")
  dy <- outer(xy[, 2], xy[, 2], "-")
  dd <- sqrt(dx * dx + dy * dy)
  k <- 1 - dd / cutoff_m
  k[k < 0] <- 0
  k
}

conley_se <- function(X, e, W, XtX_inv, hc_se) {
  U <- X * e
  meat <- crossprod(U, W %*% U)
  V <- XtX_inv %*% meat %*% XtX_inv
  s <- sqrt(pmax(diag(V), 0))
  ifelse(is.finite(s) & s > 0, pmax(s, hc_se), hc_se)
}

unit_ratios <- function(dat, cutoff_km, cfg) {
  dat <- dat[is.finite(dat$GWSlp) & is.finite(dat$dist_LB_km) &
               is.finite(dat$lon) & is.finite(dat$lat), , drop = FALSE]
  n <- nrow(dat)
  if (n < 10) return(NULL)
  if (n > cfg$n_max_wells) {
    set.seed(cfg$conley_seed * 1e6 + n)
    dat <- dat[sample.int(n, cfg$n_max_wells), , drop = FALSE]
  }
  xy <- project_local(dat$lon, dat$lat)
  y  <- dat$GWSlp; nn <- nrow(dat)
  m0 <- .hc3_fit(matrix(1, nn, 1), y)
  m1 <- .hc3_fit(cbind(1, dat$dist_LB_km), y)
  if (is.null(m0) || is.null(m1)) return(NULL)
  W <- .kernel_dense(xy, cutoff_km * 1000)
  s0 <- conley_se(m0$X, m0$e, W, m0$XtX_inv, m0$hc_se)[1]
  s1 <- conley_se(m1$X, m1$e, W, m1$XtX_inv, m1$hc_se)[2]
  tibble(cutoff_km = cutoff_km, seed = cfg$conley_seed, n_used = nn,
         hc_se_0 = m0$hc_se[1], conley_se_0 = s0, ratio_0 = s0 / m0$hc_se[1],
         hc_se_1 = m1$hc_se[2], conley_se_1 = s1, ratio_1 = s1 / m1$hc_se[2])
}

# The key must name every configuration field that changes WHICH WELLS enter
# the Conley computation. cfg$min_wells does not: it is a segment-level screen
# applied in build_segments(), long after the ratios are formed, and neither
# fs_wells() nor unit_ratios() consults it. Including it made a min_wells
# override (minw50 / minw100) miss a cache that was in fact valid, forcing a
# full 1562-unit rebuild inside a forked worker.
ratio_sample_key <- function(cfg) {
  paste0("window=",
         if (is.null(cfg$distance_window_km)) "full"
         else paste(cfg$distance_window_km, collapse = "-"),
         "|drop=", paste(sort(cfg$drop_cc), collapse = ","),
         "|nmax=", cfg$n_max_wells)
}
# Rows written before min_wells was removed from the key carry a "|minw=<n>"
# field; strip it so the existing cache stays valid.
.norm_sample_key <- function(k) sub("\\|minw=[^|]*", "", as.character(k))

build_ratios <- function(well_df, cutoff_km, cfg) {
  parts <- well_df %>% group_split(Aquifer, CC, .keep = TRUE)
  say(sprintf("Conley SEs: %d units, %g km cutoff, seed %d (cached afterwards).\n",
              length(parts), cutoff_km, cfg$conley_seed))
  one <- function(p) {
    r <- unit_ratios(p, cutoff_km, cfg)
    if (is.null(r)) return(NULL)
    r$Aquifer <- p$Aquifer[1]; r$CC <- p$CC[1]
    r
  }
  out <- par_map_dfr(parts, one, label = "Conley ratios", cfg = cfg,
                     seed = 1000L + as.integer(cfg$conley_seed))
  out %>% mutate(unit_id = paste0(Aquifer, "_", CC),
                 sample_key = ratio_sample_key(cfg))
}

.ratios_path <- function(cfg) file.path(cfg$cache_dir, "conley_ratios.csv")
.write_ratios <- function(cache, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(path, ".tmp", Sys.getpid())
  write.csv(as.data.frame(cache), tmp, row.names = FALSE)
  ok <- file.rename(tmp, path)
  if (!ok) { file.copy(tmp, path, overwrite = TRUE); unlink(tmp) }
  invisible(path)
}

get_ratios <- function(cutoff_km, cfg, wells = NULL) {
  path  <- .ratios_path(cfg)
  skey  <- ratio_sample_key(cfg)
  cache <- if (file.exists(path)) {
    as_tibble(suppressWarnings(read.csv(path, stringsAsFactors = FALSE)))
  } else tibble()
  pick <- function(k) {
    if (!nrow(cache)) return(tibble())
    cache %>% filter(.data$cutoff_km == !!cutoff_km,
                     .data$seed == !!cfg$conley_seed,
                     .norm_sample_key(.data$sample_key) ==
                       !!.norm_sample_key(k))
  }
  hit <- pick(skey); source_lab <- "exact (this well sample)"
  if (!nrow(hit)) {
    if (is.null(wells))
      stop("No cached Conley ratios for this (cutoff, seed, sample) and no ",
           "well data supplied to build them.")
    new <- build_ratios(wells, cutoff_km, cfg)
    cache <- bind_rows(cache, new)
    .write_ratios(cache, path)
    hit <- new
  }
  out <- hit %>% select(unit_id, ratio_0, ratio_1) %>%
    distinct(unit_id, .keep_all = TRUE)
  attr(out, "ratio_source") <- source_lab
  out
}

# =============================================================================
# 8. THE SEGMENT TABLE
# =============================================================================
build_segments <- function(cfg = CFG) {
  fsw <- first_stage_from_wells(cfg)
  d <- fsw %>%
    filter(!is.na(TB), n_w >= cfg$min_wells,
           is.finite(beta_0), is.finite(se_0), se_0 > 0,
           is.finite(beta_1), is.finite(se_1), se_1 > 0) %>%
    mutate(n = n_w)

  rat <- get_ratios(cfg$conley_km, cfg, wells = fs_wells(cfg))
  ratio_source <- attr(rat, "ratio_source")

  d <- d %>%
    left_join(rat, by = "unit_id", relationship = "many-to-one") %>%
    mutate(
      ratio_0  = if_else(is.finite(ratio_0) & ratio_0 >= 1, ratio_0, 1),
      ratio_1  = if_else(is.finite(ratio_1) & ratio_1 >= 1, ratio_1, 1),
      se_0_con = se_0 * ratio_0,
      se_1_con = se_1 * ratio_1,
      n_eff    = pmax(n_w / ratio_1^cfg$n_eff_power, 5),
      TBn      = as.integer(TB),
      bv_0     = se_0_con^2,
      bv_1     = se_1_con^2,
      bv_z     = 1 / pmax(n_eff - 3, 1),
      bv_bin   = bin_se^2)

  ps_cov <- if (is.null(cfg$ps_covariates)) cfg$cov_match else cfg$ps_covariates

  # Sample defined by BOTH correlation outcomes, so it is fixed across H2
  # variants (a segment with Pearson but no Spearman r, or vice versa, would
  # otherwise move the sample -- and H1 with it -- under an H2-only override).
  n_pre <- nrow(d)
  d <- d %>% filter(is.finite(z), is.finite(z_pref),
                    !is.na(aq_id), nzchar(aq_id), !is.na(CC))
  if (nrow(d) < n_pre)
    notes_add(sprintf("common-outcome filter dropped %d segment(s)", n_pre - nrow(d)))
  if (cfg$min_sd_dist_km > 0)
    d <- d %>% filter(is.finite(sd_dist), sd_dist >= cfg$min_sd_dist_km)
  if (isTRUE(cfg$exclude_river_borders)) {
    d <- d %>% filter(LB_river < 0.5)
    ps_cov <- setdiff(ps_cov, "LB_river")
  }
  if (!is.null(cfg$drop_countries))
    d <- d %>% filter(!CC %in% cfg$drop_countries)
  if (!is.null(cfg$drop_treated_segments))
    d <- d %>% filter(!unit_id %in% cfg$drop_treated_segments)
  if (!is.null(cfg$drop_aquifers))
    d <- d %>% filter(!aq_id %in% cfg$drop_aquifers)

  ps_cov <- intersect(ps_cov, names(d))
  d <- d %>% filter(if_all(all_of(ps_cov), is.finite))
  if (!nrow(d)) return(NULL)
  if (length(unique(d$TBn)) < 2) return(NULL)
  if (anyDuplicated(d$unit_id))
    stop("unit_id is not unique per row; ~1|unit_id assumes it is.")

  p <- cfg$variance_winsor
  d <- d %>% mutate(bv_0 = winsor(bv_0, p), bv_1 = winsor(bv_1, p),
                    bv_z = winsor(bv_z, p),
                    bv_bin = if (all(is.na(bv_bin))) bv_bin else winsor(bv_bin, p))

  # Covariates centred at the TRANSBOUNDARY mean: the intercept is then the
  # fitted outcome of a matched DOMESTIC segment with average shared-aquifer
  # covariates. For H2 that intercept is the placebo-type diagnostic.
  cov_all <- unique(c(cfg$cov_match, cfg$cov_extra, ps_cov))
  cov_all <- intersect(cov_all, names(d))
  for (v in cov_all) {
    mu <- mean(d[[v]][d$TBn == 1], na.rm = TRUE)
    d[[paste0(v, "_ctr")]] <- d[[v]] - mu
  }

  attr(d, "ps_cov")  <- ps_cov
  attr(d, "cov_all") <- cov_all
  attr(d, "ratio_source") <- ratio_source
  attr(d, "h2_source")    <- paste0("wells (", cfg$first_stage_h2, ")")
  d
}

# =============================================================================
# 9. DESIGN WEIGHTS (ATO preferred; ATT via full matching as alternative)
# =============================================================================
safe_glm_ps <- function(dat, ps_cov) {
  form <- reformulate(ps_cov, response = "TBn")
  fit <- tryCatch(suppressWarnings(glm(form, data = dat, family = binomial())),
                  error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$converged)) return(NULL)
  e <- fitted(fit)
  if (length(e) != nrow(dat) || any(!is.finite(e))) return(NULL)
  list(fit = fit, e = pmin(pmax(e, 1e-6), 1 - 1e-6))
}
add_overlap_weights <- function(dat, ps_cov) {
  ps <- safe_glm_ps(dat, ps_cov)
  if (is.null(ps)) return(NULL)
  dat$e     <- ps$e
  dat$w_ato <- ifelse(dat$TBn == 1, 1 - dat$e, dat$e)
  dat
}
# ATT via nearest-neighbour matching (1:2, without replacement). Full optimal
# matching requires 'optmatch', which is unavailable in the canonical run
# environment; nearest-neighbour is the documented alternative and the row is
# marked as a DIFFERENT ESTIMAND either way.
add_att_weights <- function(dat, ps_cov) {
  if (!have_match) stop("Package 'MatchIt' is required for the ATT estimand.")
  m <- tryCatch(MatchIt::matchit(reformulate(ps_cov, response = "TBn"),
                                 data = as.data.frame(dat),
                                 method = "nearest", ratio = 2),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  md <- MatchIt::match.data(m)
  if (!"unit_id" %in% names(md)) md$unit_id <- dat$unit_id[as.integer(rownames(md))]
  dat %>%
    left_join(md %>% transmute(unit_id, w_att = weights,
                               subclass = as.character(subclass)),
              by = "unit_id", relationship = "one-to-one") %>%
    mutate(w_att = if_else(is.finite(w_att), w_att, 0))
}
add_design_weights <- function(d, cfg, with_att = FALSE) {
  ps_cov <- attr(d, "ps_cov")
  d2 <- add_overlap_weights(d, ps_cov)
  if (is.null(d2)) return(NULL)
  if (with_att || cfg$estimand == "ATT") {
    d3 <- add_att_weights(d2, ps_cov)
    if (is.null(d3)) { if (cfg$estimand == "ATT") return(NULL) } else d2 <- d3
  }
  for (a in c("ps_cov", "cov_all", "ratio_source", "h2_source"))
    attr(d2, a) <- attr(d, a)
  d2$w_design <- if (cfg$estimand == "ATT") d2$w_att else d2$w_ato
  d2
}
balance_table <- function(d, cfg) {
  ps_cov  <- attr(d, "ps_cov")
  cov_rep <- unique(c(ps_cov, intersect(cfg$cov_extra, names(d))))
  has_att <- "w_att" %in% names(d)
  tibble(
    covariate  = cov_rep,
    in_ps      = cov_rep %in% ps_cov,
    unweighted = sapply(cov_rep, function(v) smd(d[[v]], d$TBn, rep(1, nrow(d)))),
    matched    = if (has_att)
      sapply(cov_rep, function(v) smd(d[[v]], d$TBn, d$w_att)) else NA_real_,
    ATO        = sapply(cov_rep, function(v) smd(d[[v]], d$TBn, d$w_ato)))
}
balance_summary <- function(d, bal) {
  tibble(
    max_abs_smd_unweighted = max(abs(bal$unweighted[bal$in_ps]), na.rm = TRUE),
    max_abs_smd_ato        = max(abs(bal$ATO[bal$in_ps]), na.rm = TRUE),
    max_abs_smd_att        = suppressWarnings(
      max(abs(bal$matched[bal$in_ps]), na.rm = TRUE)),
    ess_control_ato        = ess(d$w_ato[d$TBn == 0]),
    ess_treated_ato        = ess(d$w_ato[d$TBn == 1]),
    ess_control_att        = if ("w_att" %in% names(d))
      ess(d$w_att[d$TBn == 0]) else NA_real_,
    n_control              = sum(d$TBn == 0),
    n_treated              = sum(d$TBn == 1))
}

# =============================================================================
# 10. SECOND STAGE
# =============================================================================
# rma.mv(y, V = winsor(bv), W = winsor(w_design/(V + tau2)),
#        mods = ~ 1 + TBn + centred covariates,
#        random = list(~1|CC, ~1|unit_id), REML), iterated in tau2 to a fixed
# point. V carries the true sampling variance; W carries the design weight, so
# the realised weight is exactly w/(bv + tau2).

mods_of <- function(d, cfg) {
  ps_cov <- attr(d, "ps_cov")
  tm <- "TBn"
  if (isTRUE(cfg$outcome_adjustment) && length(ps_cov))
    tm <- c(tm, paste0(ps_cov, "_ctr"))
  if (cfg$country_effect == "fixed") tm <- c(tm, "factor(CC)")
  reformulate(tm)
}
rand_of <- function(cfg) {
  if (cfg$country_effect == "fixed") {
    if (isTRUE(cfg$include_segment_re)) list(~ 1 | unit_id) else NULL
  } else {
    if (isTRUE(cfg$include_segment_re)) list(~ 1 | CC, ~ 1 | unit_id)
    else list(~ 1 | CC)
  }
}

fit_stage2 <- function(d, ycol, bvcol, cfg, tag = "") {
  x <- d
  x$.yi <- x[[ycol]]; x$.vi <- x[[bvcol]]; x$.w0 <- x$w_design
  x <- x %>% filter(is.finite(.yi), is.finite(.vi), .vi > 0,
                    is.finite(.w0), .w0 > 0)
  if (nrow(x) < 10 || length(unique(x$TBn)) < 2) return(NULL)
  mods <- mods_of(d, cfg); rand <- rand_of(cfg)
  wp <- cfg$w_winsor
  fit_once <- function(dat) {
    # nlminb first (fast); Nelder-Mead fallback for the rare non-convergence.
    for (ctrl in list(list(optimizer = "nlminb"),
                      list(optimizer = "optim", rel.tol = 1e-8))) {
      f <- quiet(tryCatch(
        rma.mv(.yi, V = .vi, W = .W, mods = mods, random = rand, data = dat,
               method = "REML", sparse = TRUE, control = ctrl),
        error = function(e) NULL), tag)
      if (!is.null(f)) return(f)
    }
    NULL
  }
  if (!isTRUE(cfg$precision_weight)) {
    x$.W <- winsor(x$.w0, wp)
    f <- fit_once(x)
    if (is.null(f)) return(NULL)
    return(list(fit = f, data = x, tau2 = sum(f$sigma2), iters = 1L,
                converged_tau = TRUE, tau2_gap = 0))
  }
  tau2 <- 0; f <- NULL; it <- 0; converged_tau <- FALSE; gap <- NA_real_
  repeat {
    it <- it + 1
    x$.W <- winsor(x$.w0 / (x$.vi + tau2), wp)
    f <- fit_once(x)
    if (is.null(f)) return(NULL)
    tau2_new <- sum(f$sigma2)
    if (!is.finite(tau2_new)) return(NULL)
    gap <- abs(tau2_new - tau2)
    if (gap <= cfg$tol * max(1, abs(tau2_new))) {
      converged_tau <- TRUE; tau2 <- tau2_new; break
    }
    tau2 <- tau2_new
    if (it >= cfg$max_iter) break
  }
  x$.W <- winsor(x$.w0 / (x$.vi + tau2), wp)
  f <- fit_once(x)
  if (is.null(f)) return(NULL)
  if (!converged_tau)
    warning("fit_stage2: tau2 iteration hit cap for ", tag)
  list(fit = f, data = x, tau2 = tau2, iters = it,
       converged_tau = converged_tau, tau2_gap = gap)
}

# ---- inference route 2: CR2 on country --------------------------------------
robust_rma <- function(f, cluster, label = "") {
  n_cl <- length(unique(cluster)); n_ob <- length(cluster)
  cs_ok <- have_cs && (n_cl < 0.8 * n_ob)
  tries <- list(list(cs = TRUE,  tag = "CR2 clubSandwich"),
                list(cs = FALSE, tag = "CR1 metafor"))
  for (a in tries) {
    if (a$cs && !cs_ok) next
    r <- quiet(tryCatch(metafor::robust(f, cluster = cluster, clubSandwich = a$cs),
                        error = function(e) conditionMessage(e)), label)
    if (!is.character(r)) return(list(obj = r, method = a$tag))
  }
  NULL
}
rob_row <- function(rob, j) {
  if (is.null(rob) || is.na(j))
    return(list(se = NA_real_, p = NA_real_, df = NA_real_, method = "FAILED"))
  o <- rob$obj
  df <- if (!is.null(o$dfs)) as.numeric(o$dfs)[j] else
    if (!is.null(o$ddf)) as.numeric(o$ddf)[j] else NA_real_
  list(se = as.numeric(o$se)[j], p = as.numeric(o$pval)[j], df = df,
       method = rob$method)
}

# ---- inference route 3 (PRIMARY p-value): null-imposed wild cluster bootstrap
crve <- function(X, w, res, cl) {
  bread <- tryCatch(solve(crossprod(X * sqrt(w))), error = function(e) NULL)
  if (is.null(bread)) return(NULL)
  meat <- matrix(0, ncol(X), ncol(X))
  for (g in split(seq_along(res), cl)) {
    sg <- crossprod(X[g, , drop = FALSE], w[g] * res[g])
    meat <- meat + tcrossprod(sg)
  }
  G <- nlevels(factor(cl)); n <- length(res); k <- ncol(X)
  ((G / (G - 1)) * ((n - 1) / (n - k))) * bread %*% meat %*% bread
}
wcb_p <- function(y, X, w, cl, jname = "TBn", cfg = CFG) {
  wdraw <- if (cfg$wcb_weights == "webb")
    function(n) sample(c(-sqrt(1.5), -1, -sqrt(0.5),
                         sqrt(0.5), 1, sqrt(1.5)), n, replace = TRUE)
  else function(n) sample(c(-1, 1), n, replace = TRUE)
  j <- which(colnames(X) == jname)
  if (!length(j)) return(NULL)
  cl <- factor(cl)
  if (nlevels(cl) < 3) return(NULL)
  full <- tryCatch(lm.wfit(X, y, w), error = function(e) NULL)
  if (is.null(full)) return(NULL)
  V <- crve(X, w, as.vector(full$residuals), cl)
  if (is.null(V) || !is.finite(V[j, j]) || V[j, j] <= 0) return(NULL)
  t_obs <- unname(full$coefficients[j] / sqrt(V[j, j]))
  Xr <- X[, -j, drop = FALSE]
  fr <- tryCatch(lm.wfit(Xr, y, w), error = function(e) NULL)
  if (is.null(fr)) return(NULL)
  res_r <- as.vector(fr$residuals); yhat_r <- y - res_r
  gs <- split(seq_along(y), cl)
  st_wcb <- .rng_save(); on.exit(.rng_restore(st_wcb), add = TRUE)
  set.seed(cfg$wcb_seed, kind = "Mersenne-Twister")
  tb <- replicate(cfg$n_wcb, {
    v <- wdraw(length(gs))
    ystar <- yhat_r
    for (i in seq_along(gs)) ystar[gs[[i]]] <- yhat_r[gs[[i]]] + v[i] * res_r[gs[[i]]]
    fb <- tryCatch(lm.wfit(X, ystar, w), error = function(e) NULL)
    if (is.null(fb)) return(NA_real_)
    Vb <- crve(X, w, as.vector(fb$residuals), cl)
    if (is.null(Vb) || !is.finite(Vb[j, j]) || Vb[j, j] <= 0) return(NA_real_)
    unname(fb$coefficients[j] / sqrt(Vb[j, j]))
  })
  tb <- tb[is.finite(tb)]
  if (!length(tb)) return(NULL)
  n <- length(tb)
  one <- function(alt) switch(alt,
    greater = (1 + sum(tb >= t_obs)) / (n + 1),
    less    = (1 + sum(tb <= t_obs)) / (n + 1),
    NA_real_)
  list(p = (1 + sum(abs(tb) >= abs(t_obs))) / (n + 1),
       p_greater = one("greater"), p_less = one("less"),
       t_obs = t_obs, B_ok = n, G = nlevels(cl))
}

# Wild-cluster-bootstrap CI by test inversion: the set of null values the
# null-imposed bootstrap does not reject. Consistent with p_wcb by
# construction. Expensive (about 24 bootstraps per coefficient), so gated
# behind cfg$wcb_ci and used only for the reported main specification.
wcb_ci <- function(y, X, w, cl, jname = "TBn", cfg = CFG, level = cfg$ci_level,
                   max_steps = 12) {
  a <- 1 - level
  j <- which(colnames(X) == jname)
  if (!length(j)) return(NULL)
  full <- tryCatch(lm.wfit(X, y, w), error = function(e) NULL)
  if (is.null(full)) return(NULL)
  b_hat <- unname(full$coefficients[j])
  V <- crve(X, w, as.vector(full$residuals), factor(cl))
  if (is.null(V) || !is.finite(V[j, j]) || V[j, j] <= 0) return(NULL)
  se <- sqrt(V[j, j])
  p_at <- function(b) {
    r <- wcb_p(y - b * X[, j], X, w, cl, jname, cfg)
    if (is.null(r)) NA_real_ else r$p
  }
  if (!is.finite(p_at(b_hat))) return(NULL)
  side <- function(dir) {
    lo <- b_hat; hi <- b_hat + dir * 2 * se; k <- 0
    while (k < 8 && !is.na(p_at(hi)) && p_at(hi) > a) {
      lo <- hi; hi <- b_hat + dir * 2^(k + 2) * se; k <- k + 1
    }
    if (is.na(p_at(hi)) || p_at(hi) > a) return(NA_real_)
    for (i in seq_len(max_steps)) {
      mid <- (lo + hi) / 2
      pm <- p_at(mid)
      if (is.na(pm)) return(NA_real_)
      if (pm > a) lo <- mid else hi <- mid
    }
    (lo + hi) / 2
  }
  list(estimate = b_hat, lo = side(-1), hi = side(+1), level = level)
}

# One-sided wild-cluster-bootstrap bound by test inversion: mirrors wcb_ci()
# above but inverts the DIRECTIONAL bootstrap p (p_greater / p_less), and
# searches only the side away from zero in the pre-stated direction. For
# alt = "greater" this returns the lower limit of the one-sided (1 - alpha)
# confidence set [bound, +Inf); for alt = "less" the upper limit of
# (-Inf, bound]. Because it inverts the SAME directional p that one_sided_p()
# reports as the primary p-value, the bound and the p-value cannot disagree
# about significance by construction -- unlike a two-sided CR2/Wald interval,
# which answers a different inference question. Gated behind cfg$wcb_ci, like
# the two-sided version.
wcb_ci_one_sided <- function(y, X, w, cl, jname = "TBn", cfg = CFG, alt,
                             level = 1 - cfg$alpha_one, max_steps = 12) {
  a <- 1 - level
  if (!identical(alt, "greater") && !identical(alt, "less")) return(NULL)
  j <- which(colnames(X) == jname)
  if (!length(j)) return(NULL)
  full <- tryCatch(lm.wfit(X, y, w), error = function(e) NULL)
  if (is.null(full)) return(NULL)
  b_hat <- unname(full$coefficients[j])
  V <- crve(X, w, as.vector(full$residuals), factor(cl))
  if (is.null(V) || !is.finite(V[j, j]) || V[j, j] <= 0) return(NULL)
  se <- sqrt(V[j, j])
  # dir = -1 (search below b_hat) gives the LOWER bound for alt="greater";
  # dir = +1 (search above b_hat) gives the UPPER bound for alt="less".
  dir <- if (identical(alt, "greater")) -1 else 1
  p_at <- function(b) {
    r <- wcb_p(y - b * X[, j], X, w, cl, jname, cfg)
    if (is.null(r)) return(NA_real_)
    if (identical(alt, "greater")) r$p_greater else r$p_less
  }
  if (!is.finite(p_at(b_hat))) return(NULL)
  lo <- b_hat; hi <- b_hat + dir * 2 * se; k <- 0
  while (k < 8 && !is.na(p_at(hi)) && p_at(hi) > a) {
    lo <- hi; hi <- b_hat + dir * 2^(k + 2) * se; k <- k + 1
  }
  if (is.na(p_at(hi)) || p_at(hi) > a) return(NULL)
  for (i in seq_len(max_steps)) {
    mid <- (lo + hi) / 2
    pm <- p_at(mid)
    if (is.na(pm)) return(NULL)
    if (pm > a) lo <- mid else hi <- mid
  }
  list(bound = (lo + hi) / 2, level = level, alt = alt)
}

# ---- pull the TBn row and intercept, all routes attached --------------------
summarise_fit <- function(o, cfg, lab = "", hyp = NA_character_) {
  if (is.null(o)) return(NULL)
  f <- o$fit; x <- o$data
  jt <- grep("^TBn$", rownames(f$beta))[1]
  ji <- grep("^intrcpt$", rownames(f$beta))[1]
  if (is.na(jt)) return(NULL)
  b  <- as.numeric(f$beta); se <- sqrt(diag(vcov(f)))
  aligned <- !is.null(f$X) && nrow(f$X) == nrow(x)
  if (!aligned)
    warning("summarise_fit: model/data row mismatch; clustered inference skipped for ", lab)
  rob <- if (aligned) robust_rma(f, x$CC, lab) else NULL
  rcc <- rob_row(rob, jt)
  icc <- rob_row(rob, ji)
  wb  <- if (aligned) wcb_p(x$.yi, f$X, x$.W, x$CC, cfg = cfg) else NULL
  wci <- if (aligned && isTRUE(cfg$wcb_ci))
    wcb_ci(x$.yi, f$X, x$.W, x$CC, cfg = cfg) else NULL
  alt1 <- if (!is.na(hyp)) cfg$directional_alt[[hyp]] else NULL
  wci1 <- if (aligned && isTRUE(cfg$wcb_ci) && !is.null(alt1))
    wcb_ci_one_sided(x$.yi, f$X, x$.W, x$CC, cfg = cfg, alt = alt1) else NULL
  cv  <- crit_val(NA_real_, cfg$ci_level)
  list(
    lab = lab, fit = f, data = x, tau2 = o$tau2, iters = o$iters,
    converged_tau = isTRUE(o$converged_tau), aligned = aligned,
    est = b[jt], se_model = se[jt],
    p_model = 2 * pnorm(abs(b[jt] / se[jt]), lower.tail = FALSE),
    ci_lo = b[jt] - cv * se[jt], ci_hi = b[jt] + cv * se[jt],
    se_cc = rcc$se, p_cc = rcc$p, df_cc = rcc$df, method_cc = rcc$method,
    cc_lo = if (is.finite(rcc$se)) b[jt] - crit_val(rcc$df, cfg$ci_level) * rcc$se else NA_real_,
    cc_hi = if (is.finite(rcc$se)) b[jt] + crit_val(rcc$df, cfg$ci_level) * rcc$se else NA_real_,
    p_wcb = if (is.null(wb)) NA_real_ else wb$p,
    p_wcb_greater = if (is.null(wb)) NA_real_ else wb$p_greater,
    p_wcb_less    = if (is.null(wb)) NA_real_ else wb$p_less,
    wcb_lo = if (is.null(wci)) NA_real_ else wci$lo,
    wcb_hi = if (is.null(wci)) NA_real_ else wci$hi,
    wcb_bound_one = if (is.null(wci1)) NA_real_ else wci1$bound,
    wcb_bound_one_alt = if (is.null(wci1)) NA_character_ else wci1$alt,
    G_wcb = if (is.null(wb)) NA_integer_ else wb$G,
    b0 = if (!is.na(ji)) b[ji] else NA_real_,
    b0_se_model = if (!is.na(ji)) se[ji] else NA_real_,
    b0_p_model = if (!is.na(ji))
      2 * pnorm(abs(b[ji] / se[ji]), lower.tail = FALSE) else NA_real_,
    b0_se_cc = icc$se, b0_p_cc = icc$p, b0_df_cc = icc$df,
    b0_lo = if (!is.na(ji)) b[ji] - (if (is.finite(icc$se))
      crit_val(icc$df, cfg$ci_level) * icc$se else cv * se[ji]) else NA_real_,
    b0_hi = if (!is.na(ji)) b[ji] + (if (is.finite(icc$se))
      crit_val(icc$df, cfg$ci_level) * icc$se else cv * se[ji]) else NA_real_,
    QE = if (!is.null(f$QE)) as.numeric(f$QE) else NA_real_,
    QE_df = if (!is.null(f$QEdf)) as.numeric(f$QEdf)[1] else NA_real_,
    QE_p = if (!is.null(f$QEp)) as.numeric(f$QEp) else NA_real_,
    sigma2 = f$sigma2)
}

# One-sided p on the pre-stated side. Primary: wild cluster bootstrap
# (few-cluster-valid); fallback: CR2 t. The source used is recorded.
one_sided_p <- function(s, hyp, cfg = CFG) {
  alt <- cfg$directional_alt[[hyp]]
  if (is.null(alt)) return(list(p = NA_real_, source = "none (no direction stated)"))
  pw <- if (identical(alt, "greater")) s$p_wcb_greater else s$p_wcb_less
  if (is.finite(pw))
    return(list(p = pw, source = "wild cluster bootstrap (country, null-imposed)"))
  if (is.finite(s$se_cc) && is.finite(s$df_cc)) {
    t <- s$est / s$se_cc
    p <- if (identical(alt, "greater")) pt(t, s$df_cc, lower.tail = FALSE)
         else pt(t, s$df_cc)
    return(list(p = p, source = "CR2 country-clustered t"))
  }
  t <- s$est / s$se_model
  p <- if (identical(alt, "greater")) pnorm(t, lower.tail = FALSE) else pnorm(t)
  list(p = p, source = "model-based Wald (fallback)")
}

# Preferred SE for reporting/MDE and its provenance.
se_of <- function(s) {
  if (is.finite(s$se_cc)) list(se = s$se_cc, source = "CR2 country-clustered SE")
  else list(se = s$se_model, source = "model-based SE (CR2 unavailable)")
}

# =============================================================================
# 11. MDE -- MINIMUM DETECTABLE EFFECT (H1 detection limit)
# =============================================================================
# For the pre-stated ONE-SIDED test at alpha_one, at power `power`:
#   MDE = (qnorm(1 - alpha_one) + qnorm(power)) * SE.
# This is the smallest true effect in the predicted direction that the design
# would detect with the stated probability. It is NOT the confidence interval,
# and the two are exported as distinct quantities.
mde_of <- function(se, power, cfg = CFG) {
  stopifnot(is.finite(se), se > 0, power > 0.5, power < 1)
  (qnorm(1 - cfg$alpha_one) + qnorm(power)) * se
}
# Assertion required by the brief: verify alpha and power enter as intended.
assert_mde <- function(cfg = CFG) {
  se <- 1
  expect80 <- qnorm(1 - cfg$alpha_one) + qnorm(0.80)
  stopifnot(abs(mde_of(se, 0.80, cfg) - expect80) < 1e-12)
  stopifnot(abs(qnorm(1 - cfg$alpha_one) - qnorm(0.95)) < 1e-12) # alpha_one = .05
  invisible(TRUE)
}

# =============================================================================
# 12. run_specification() -- THE ONE ENTRY POINT
# =============================================================================
run_specification <- function(spec_id,
                              specification = spec_id,
                              family        = "Design",
                              cfg           = CFG,
                              run_h3        = FALSE,
                              with_att      = FALSE,
                              ...) {
  t0 <- Sys.time()
  notes0 <- notes_get()
  cfg <- cfg_with(cfg, ...)

  blank <- list(estimate = NA_real_, ci90_low = NA_real_, ci90_high = NA_real_,
                p_model = NA_real_, p_cr2 = NA_real_, p_wild = NA_real_,
                wcb_lo = NA_real_, wcb_hi = NA_real_,
                wcb_bound_one = NA_real_, wcb_bound_one_alt = NA_character_,
                p_one = NA_real_, p_one_source = NA_character_,
                se_model = NA_real_, se_cr2 = NA_real_, se_source = NA_character_,
                mde80 = NA_real_, mde90 = NA_real_,
                baseline = NA_real_, baseline_p = NA_real_, tau2 = NA_real_)

  fail <- function(msg) list(
    spec_id = spec_id, specification = specification, family = family,
    estimand = cfg$estimand, ok = FALSE, h1 = blank, h2 = blank,
    sample = list(n_wells = NA_integer_, n_segments = NA_integer_,
                  n_treated = NA_integer_, n_countries = NA_integer_,
                  n_aquifers = NA_integer_),
    diagnostics = list(converged = FALSE, warnings = msg,
                       runtime_seconds = as.numeric(difftime(Sys.time(), t0,
                                                             units = "secs"))),
    h3 = NULL, cfg = cfg)

  res <- tryCatch({
    d <- build_segments(cfg)
    if (is.null(d)) return(fail("empty analysis sample"))
    d <- add_design_weights(d, cfg, with_att = with_att)
    if (is.null(d)) return(fail("propensity / matching model failed"))

    h2o <- h2_outcome(cfg)
    o1 <- fit_stage2(d, "beta_0", "bv_0", cfg, tag = paste(spec_id, "H1"))
    o2 <- fit_stage2(d, h2o$y,   h2o$bv,  cfg, tag = paste(spec_id, "H2"))
    s1 <- summarise_fit(o1, cfg, paste(spec_id, "H1"), hyp = "H1")
    s2 <- summarise_fit(o2, cfg, paste(spec_id, "H2"), hyp = "H2")
    if (is.null(s1) && is.null(s2)) return(fail("both second-stage fits failed"))

    pack <- function(s, hyp) {
      if (is.null(s)) return(blank)
      po <- one_sided_p(s, hyp, cfg)
      sp <- se_of(s)
      list(estimate = s$est,
           ci90_low  = if (is.finite(s$cc_lo)) s$cc_lo else s$ci_lo,
           ci90_high = if (is.finite(s$cc_hi)) s$cc_hi else s$ci_hi,
           p_model = s$p_model, p_cr2 = s$p_cc, p_wild = s$p_wcb,
           wcb_lo = s$wcb_lo %||% NA_real_, wcb_hi = s$wcb_hi %||% NA_real_,
           wcb_bound_one = s$wcb_bound_one %||% NA_real_,
           wcb_bound_one_alt = s$wcb_bound_one_alt %||% NA_character_,
           p_one = po$p, p_one_source = po$source,
           se_model = s$se_model, se_cr2 = s$se_cc, df_cr2 = s$df_cc,
           se_source = sp$source,
           mde80 = if (hyp == "H1") mde_of(sp$se, 0.80, cfg) else NA_real_,
           mde90 = if (hyp == "H1") mde_of(sp$se, 0.90, cfg) else NA_real_,
           baseline = s$b0, baseline_p = s$b0_p_cc,
           baseline_p_model = s$b0_p_model,
           baseline_lo = s$b0_lo, baseline_hi = s$b0_hi,
           tau2 = s$tau2, QE = s$QE, QE_df = s$QE_df, QE_p = s$QE_p,
           iters = s$iters, converged_tau = isTRUE(s$converged_tau))
    }

    h3 <- if (isTRUE(run_h3)) h3_analysis(d, cfg) else NULL

    list(spec_id = spec_id, specification = specification, family = family,
         estimand = cfg$estimand, ok = TRUE,
         h1 = pack(s1, "H1"), h2 = pack(s2, "H2"),
         h2_unit = h2o$unit, h2_scale = h2o$scale,
         sample = list(n_wells = sum(d$n, na.rm = TRUE),
                       n_segments = nrow(d),
                       n_treated = sum(d$TBn == 1),
                       n_control = sum(d$TBn == 0),
                       n_countries = n_distinct(d$CC),
                       n_aquifers = n_distinct(d$aq_id)),
         diagnostics = list(
           converged = !is.null(s1) && !is.null(s2) &&
             isTRUE(s1$converged_tau) && isTRUE(s2$converged_tau),
           warnings = paste(setdiff(notes_get(), notes0), collapse = " | "),
           ratio_source = attr(d, "ratio_source"),
           h2_source = attr(d, "h2_source"),
           uncertainty_estimator_matched = cfg$first_stage_h1 == "mean",
           ps_covariates = paste(attr(d, "ps_cov"), collapse = "+"),
           runtime_seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))),
         h3 = h3, data = d, s1 = s1, s2 = s2, cfg = cfg)
  }, error = function(e) fail(paste("ERROR:", conditionMessage(e))))
  res
}

h2_outcome <- function(cfg) {
  switch(cfg$first_stage_h2,
    pearson        = list(y = "z", bv = "bv_z", unit = "Fisher z", scale = "z"),
    spearman       = list(y = "z", bv = "bv_z", unit = "Fisher z (Spearman)", scale = "z"),
    physical_slope = list(y = "beta_1", bv = "bv_1", unit = "mm/yr per km", scale = "slope"),
    distance_bins  = list(y = "bin_diff", bv = "bv_bin",
                          unit = "mm/yr, near minus interior", scale = "bin"),
    stop("unknown first_stage_h2: ", cfg$first_stage_h2))
}

# Flatten one result into a table row.
spec_row <- function(r) {
  g <- function(x, f, default = NA_real_) {
    v <- x[[f]]; if (is.null(v) || !length(v)) default else v[1]
  }
  gc <- function(x, f) { v <- x[[f]]; if (is.null(v)) NA_character_ else v[1] }
  # NULL identifiers must become NA, not vanish: tibble() DROPS a NULL column,
  # so a malformed result would otherwise yield a row with no spec_id and the
  # failure would surface far downstream as "Column `spec_id` doesn't exist".
  if (!is.list(r)) r <- list()
  tibble(
    spec_id = gc(r, "spec_id"), family = gc(r, "family"),
    specification = gc(r, "specification"),
    estimand = gc(r, "estimand"), ok = isTRUE(r$ok),
    h1_estimate = g(r$h1, "estimate"),
    h1_ci90_low = g(r$h1, "ci90_low"), h1_ci90_high = g(r$h1, "ci90_high"),
    h1_p_model = g(r$h1, "p_model"), h1_p_cr2 = g(r$h1, "p_cr2"),
    h1_p_wild = g(r$h1, "p_wild"), h1_p_one = g(r$h1, "p_one"),
    h1_p_one_source = gc(r$h1, "p_one_source"),
    h1_se_source = gc(r$h1, "se_source"),
    h1_mde80 = g(r$h1, "mde80"), h1_mde90 = g(r$h1, "mde90"),
    h1_tau2 = g(r$h1, "tau2"),
    h2_estimate = g(r$h2, "estimate"),
    h2_ci90_low = g(r$h2, "ci90_low"), h2_ci90_high = g(r$h2, "ci90_high"),
    h2_p_model = g(r$h2, "p_model"), h2_p_cr2 = g(r$h2, "p_cr2"),
    h2_p_wild = g(r$h2, "p_wild"), h2_p_one = g(r$h2, "p_one"),
    h2_p_one_source = gc(r$h2, "p_one_source"),
    h2_tau2 = g(r$h2, "tau2"),
    h2_baseline = g(r$h2, "baseline"), h2_baseline_p = g(r$h2, "baseline_p"),
    h2_unit = r$h2_unit %||% NA_character_,
    h2_scale = r$h2_scale %||% NA_character_,
    n_wells = g(r$sample, "n_wells"), n_segments = g(r$sample, "n_segments"),
    n_treated = g(r$sample, "n_treated"),
    n_control = g(r$sample, "n_control", NA_integer_),
    n_countries = g(r$sample, "n_countries"),
    n_aquifers = g(r$sample, "n_aquifers"),
    converged = isTRUE(r$diagnostics$converged),
    ratio_source = r$diagnostics$ratio_source %||% NA_character_,
    h2_source = r$diagnostics$h2_source %||% NA_character_,
    uncertainty_estimator_matched =
      r$diagnostics$uncertainty_estimator_matched %||% NA,
    runtime_seconds = g(r$diagnostics, "runtime_seconds"),
    warnings = r$diagnostics$warnings %||% "")
}

# =============================================================================
# 13. H3 -- HETEROGENEITY ACROSS SHARED AQUIFERS
# =============================================================================
fit_h3_model <- function(d, ycol, bvcol, cfg, tag = "") {
  mods <- if (cfg$h3_deviations == "adjusted") mods_of(d, cfg) else ~ 1 + TBn
  rand <- list(~ 1 | CC, ~ 1 | unit_id)
  if (cfg$h3_variance %in% c("legacy", "precision")) {
    x <- d
    x$.yi <- x[[ycol]]
    x$.vi <- if (cfg$h3_variance == "legacy")
      winsor(x[[bvcol]] / pmax(x$w_ato, 1e-12), cfg$variance_winsor)
    else winsor(x[[bvcol]], cfg$variance_winsor)
    x <- x %>% filter(is.finite(.yi), is.finite(.vi), .vi > 0)
    f <- quiet(tryCatch(rma.mv(.yi, V = .vi, mods = mods, random = rand,
                               data = x, method = "REML", sparse = TRUE),
                        error = function(e) NULL), tag)
    if (is.null(f)) return(NULL)
    o <- list(fit = f, data = x, tau2 = sum(f$sigma2), iters = 1L,
              converged_tau = TRUE, tau2_gap = 0)
  } else {
    cfg2 <- cfg_with(cfg, outcome_adjustment = (cfg$h3_deviations == "adjusted"),
                     country_effect = "random", include_segment_re = TRUE)
    o <- fit_stage2(d, ycol, bvcol, cfg2, tag = tag)
    if (is.null(o)) return(NULL)
  }
  f <- o$fit; x <- o$data
  j <- grep("^TBn$", rownames(f$beta))[1]
  if (is.na(j)) return(NULL)
  rr <- rob_row(robust_rma(f, x$aq_id, tag), j)
  list(fit = f, data = x, est = as.numeric(f$beta)[j],
       gamma_se = if (is.finite(rr$se)) rr$se else sqrt(diag(vcov(f)))[j],
       tau2_cc = f$sigma2[1],
       tau2_unit = if (length(f$sigma2) >= 2) f$sigma2[2] else 0)
}

valid_fdr <- function(f, x, cfg) {
  s    <- f$sigma2
  t_cc <- if (length(s) >= 1) s[1] else 0
  t_un <- if (length(s) >= 2) s[2] else 0
  fixed <- as.vector(f$X %*% as.numeric(f$beta))
  zi   <- (x$.yi - fixed) / sqrt(x$.vi + t_un + t_cc)
  gid  <- x$unit_id
  tb   <- x$TBn == 1
  p    <- 2 * pnorm(abs(zi[tb]), lower.tail = FALSE)
  q    <- p.adjust(p, "BH")
  R    <- sum(q < cfg$q_cut, na.rm = TRUE); m <- length(q)
  lev  <- if (R > 0) 1 - cfg$q_cut * R / m else 0.95
  tibble(.gid = gid[tb], z_std = zi[tb], p_std = p, q_std = q,
         fcr_z = qnorm(1 - (1 - lev) / 2), fcr_level = lev,
         n_selected = R, n_tested = m)
}

blup_direct <- function(f, x, res) {
  n <- nrow(x)
  if (is.null(f$X) || nrow(f$X) != n) return(NULL)
  ccf <- factor(x$CC)
  Zc  <- stats::model.matrix(~ 0 + ccf)
  s2cc <- res$tau2_cc; s2u <- res$tau2_unit
  M <- s2cc * tcrossprod(Zc) + diag(s2u + x$.vi, n, n)
  r <- as.numeric(x$.yi - as.vector(f$X %*% as.numeric(f$beta)))
  Mi_r <- tryCatch(solve(M, r), error = function(e) NULL)
  if (is.null(Mi_r)) return(NULL)
  list(u_cc = as.numeric(Zc %*% (s2cc * crossprod(Zc, Mi_r))),
       u_k  = as.numeric(s2u * Mi_r))
}

validate_blups <- function(f, x, u_cc, u_k, cfg, res = NULL) {
  out <- list(ok = NA, max_abs_diff_unit = NA_real_, max_abs_diff_cc = NA_real_,
              method = "direct GLS", msg = "not attempted")
  if (!isTRUE(cfg$h3_blup_validate)) { out$msg <- "disabled"; return(out) }
  dd <- if (is.null(res)) NULL else blup_direct(f, x, res)
  if (is.null(dd)) { out$msg <- "direct GLS BLUP could not be computed"; return(out) }
  out$max_abs_diff_cc   <- max(abs(dd$u_cc - u_cc), na.rm = TRUE)
  out$max_abs_diff_unit <- max(abs(dd$u_k  - u_k),  na.rm = TRUE)
  worst <- max(c(out$max_abs_diff_cc, out$max_abs_diff_unit), na.rm = TRUE)
  scale <- max(1, stats::sd(x$.yi, na.rm = TRUE))
  out$ok  <- worst <= cfg$h3_blup_tol * scale
  out$msg <- sprintf("max |closed form - GLS| = %.3g (tol %.3g)",
                     worst, cfg$h3_blup_tol * scale)
  if (!isTRUE(out$ok))
    warning("H3 BLUP validation FAILED: ", out$msg)
  out
}

h3_blups <- function(res, label, cfg) {
  if (is.null(res)) return(NULL)
  x <- res$data; f <- res$fit
  b <- as.numeric(f$beta)
  fixed <- as.vector(f$X %*% b)
  r0 <- x$.yi - fixed
  cc <- as.integer(factor(x$CC)); wj <- 1 / (x$.vi + res$tau2_unit)
  u_cc <- ((res$tau2_cc * as.vector(rowsum(r0 * wj, cc))) /
             (1 + res$tau2_cc * as.vector(rowsum(wj, cc))))[cc]
  shrink <- res$tau2_unit / (res$tau2_unit + x$.vi)
  u_k    <- shrink * (r0 - u_cc)
  se_u   <- sqrt(pmax(res$tau2_unit * x$.vi / (res$tau2_unit + x$.vi), 0))
  vb <- validate_blups(f, x, u_cc, u_k, cfg, res = res)
  fdr <- valid_fdr(f, x, cfg)
  out <- x %>%
    transmute(.gid = unit_id,
              unit_id, aq_id, Aquifer, countries = CC,
              label_txt = paste0(Aquifer, " (", CC, ")"),
              n_wells = .data$n, shrinkage = shrink,
              deviation = u_k,
              excess = if (cfg$h3_centre == "domestic") res$est + u_k else u_k,
              post_se = sqrt(se_u^2 + res$gamma_se^2), TBn) %>%
    filter(TBn == 1) %>%
    left_join(fdr, by = ".gid") %>%
    mutate(screen = !is.na(q_std) & q_std < cfg$q_cut,
           exc_lo = excess - coalesce(fcr_z, 1.96) * post_se,
           exc_hi = excess + coalesce(fcr_z, 1.96) * post_se,
           outcome = label, tau2 = res$tau2_unit, centre = cfg$h3_centre) %>%
    arrange(desc(abs(excess)))
  attr(out, "blup_validation") <- vb
  out
}

h3_analysis <- function(d, cfg = CFG) {
  h2o <- h2_outcome(cfg)
  m0 <- fit_h3_model(d, "beta_0", "bv_0", cfg, "H3 level")
  mz <- fit_h3_model(d, h2o$y,   h2o$bv,  cfg, "H3 gradient")
  eb0 <- h3_blups(m0, "level", cfg)
  ebz <- h3_blups(mz, "gradient", cfg)
  if (is.null(eb0) || is.null(ebz)) return(NULL)
  .vrow <- function(x, lab) {
    v <- attr(x, "blup_validation")
    if (is.null(v)) v <- list(ok = NA, max_abs_diff_unit = NA_real_,
                              max_abs_diff_cc = NA_real_,
                              method = "none", msg = "not recorded")
    as_tibble(c(list(outcome = lab), v))
  }
  blup_check <- bind_rows(.vrow(eb0, "level"), .vrow(ebz, "gradient"))

  sc <- eb0 %>%
    select(.gid, unit_id, aq_id, Aquifer, countries, label_txt, n_wells,
           dev_level = deviation, exc_level = excess,
           lo_level = exc_lo, hi_level = exc_hi, flag_level = screen,
           z_level = z_std, q_level = q_std) %>%
    inner_join(ebz %>% select(.gid, dev_gradient = deviation,
                              exc_gradient = excess,
                              lo_gradient = exc_lo, hi_gradient = exc_hi,
                              flag_gradient = screen,
                              z_gradient = z_std, q_gradient = q_std),
               by = ".gid") %>%
    mutate(
      pattern = case_when(
        exc_level >  0 & exc_gradient <  0 ~ "fast + borderward",
        exc_level >  0 & exc_gradient >= 0 ~ "fast, not borderward",
        exc_level <= 0 & exc_gradient <  0 ~ "borderward, not fast",
        TRUE                               ~ "neither"),
      flagged = flag_level | flag_gradient,
      rank_level    = rank(-abs(exc_level),    ties.method = "min"),
      rank_gradient = rank(-abs(exc_gradient), ties.method = "min"),
      shortlist     = rank_level <= cfg$top_k | rank_gradient <= cfg$top_k)

  het <- bind_rows(
    tibble(outcome = "level (mm/yr)",
           tau_unit = sqrt(max(m0$tau2_unit, 0)), tau_cc = sqrt(max(m0$tau2_cc, 0)),
           contrast = m0$est,
           QE = as.numeric(m0$fit$QE %||% NA), QE_df = as.numeric(m0$fit$QEdf %||% NA)[1],
           QE_p = as.numeric(m0$fit$QEp %||% NA)),
    tibble(outcome = paste0("gradient (", h2o$unit, ")"),
           tau_unit = sqrt(max(mz$tau2_unit, 0)), tau_cc = sqrt(max(mz$tau2_cc, 0)),
           contrast = mz$est,
           QE = as.numeric(mz$fit$QE %||% NA), QE_df = as.numeric(mz$fit$QEdf %||% NA)[1],
           QE_p = as.numeric(mz$fit$QEp %||% NA))) %>%
    mutate(ratio_tau_over_contrast = tau_unit / abs(contrast))

  rp <- suppressWarnings(cor(sc$exc_level, sc$exc_gradient, use = "complete.obs"))
  rs <- suppressWarnings(cor(sc$exc_level, sc$exc_gradient, method = "spearman",
                             use = "complete.obs"))
  bs <- with_seed(11, replicate(2000, {
    i <- sample(nrow(sc), replace = TRUE)
    suppressWarnings(cor(sc$exc_level[i], sc$exc_gradient[i], use = "complete.obs"))
  }))
  bs <- bs[is.finite(bs)]
  tab <- table(fast = sc$exc_level > 0, borderward = sc$exc_gradient < 0)
  ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
  cooc <- tibble(
    pearson = rp,
    pearson_lo = if (length(bs)) unname(quantile(bs, 0.025)) else NA_real_,
    pearson_hi = if (length(bs)) unname(quantile(bs, 0.975)) else NA_real_,
    spearman = rs,
    fisher_p = if (is.null(ft)) NA_real_ else ft$p.value,
    odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate))

  list(scatter = sc, eb_level = eb0, eb_gradient = ebz, het = het,
       blup_validation = blup_check, cooccurrence = cooc,
       model_level = m0, model_gradient = mz,
       contrast_level = m0$est, contrast_gradient = mz$est,
       h2_unit = h2o$unit, variance_path = cfg$h3_variance)
}

# ---- the coded H3 selection rule --------------------------------------------
raw_descriptives <- function(d, cfg = CFG) {
  dom <- d %>% filter(TBn == 0)
  c0  <- mean(dom$beta_0, na.rm = TRUE)
  cz  <- mean(dom$z,      na.rm = TRUE)
  out <- d %>%
    filter(TBn == 1) %>%
    transmute(.gid = unit_id, unit_id, aq_id, Aquifer, CC,
              label_txt = paste0(Aquifer, " (", CC, ")"),
              n_wells = n,
              depletion = beta_0, gradient = z,
              dep_vs_domestic  = beta_0 - c0,
              grad_vs_domestic = z - cz) %>%
    mutate(
      zl = as.numeric(scale(dep_vs_domestic)),
      zg = as.numeric(scale(grad_vs_domestic)),
      dev_raw   = sqrt(zl^2 + zg^2),
      rank_raw  = rank(-dev_raw, ties.method = "min"),
      quadrant  = dep_vs_domestic > 0 & grad_vs_domestic < 0) %>%
    arrange(rank_raw)
  attr(out, "domestic_mean") <- c(level = c0, gradient = cz)
  out
}

path_scatters <- function(d, cfg = CFG, paths) {
  map_dfr(paths, function(v) {
    h <- tryCatch(h3_analysis(d, cfg_with(cfg, h3_variance = v)),
                  error = function(e) NULL)
    if (is.null(h)) { say("  H3 variance path '", v, "' failed to fit.\n")
                      return(NULL) }
    h$scatter %>%
      transmute(.gid, label_txt, Aquifer, countries, n_wells, path = v,
                exc_level, exc_gradient, rank_level, rank_gradient,
                q_level, q_gradient,
                zl = as.numeric(scale(exc_level)),
                zg = as.numeric(scale(exc_gradient))) %>%
      mutate(dev = sqrt(zl^2 + zg^2),
             rank_dev = rank(-dev, ties.method = "min"),
             quadrant = exc_level > 0 & exc_gradient < 0)
  })
}

consistency_selection <- function(long, voting, k) {
  nv <- length(voting)
  long %>%
    filter(path %in% voting) %>%
    group_by(.gid, label_txt, Aquifer, countries, n_wells) %>%
    summarise(
      n_paths          = n(),
      always_top_level = sum(rank_level    <= k) == nv,
      always_top_grad  = sum(rank_gradient <= k) == nv,
      always_top_dev   = sum(rank_dev      <= k) == nv,
      always_quadrant  = sum(quadrant) == nv,
      worst_rank_dev   = max(rank_dev),
      median_exc_level = median(exc_level),
      median_exc_grad  = median(exc_gradient),
      min_q_level      = suppressWarnings(min(q_level, na.rm = TRUE)),
      min_q_gradient   = suppressWarnings(min(q_gradient, na.rm = TRUE)),
      .groups = "drop") %>%
    mutate(consistent = always_top_level | always_top_grad |
             always_top_dev | always_quadrant,
           n_criteria = always_top_level + always_top_grad +
             always_top_dev + always_quadrant) %>%
    arrange(worst_rank_dev)
}

identify_segments <- function(d, cfg = CFG, k = cfg$top_k,
                              paths_voting = cfg$h3_paths_voting,
                              paths_shown  = cfg$h3_paths_shown) {
  raw  <- raw_descriptives(d, cfg)
  long <- path_scatters(d, cfg, paths_shown)
  if (!nrow(long)) return(NULL)
  voting <- intersect(paths_voting, unique(long$path))
  if (length(voting) < 2)
    stop("fewer than two voting variance paths fitted; no consistency to assess.")
  sel <- consistency_selection(long, voting, k)
  final <- raw %>%
    select(.gid, label_txt, n_wells, dep_vs_domestic, grad_vs_domestic,
           rank_raw, quadrant_raw = quadrant) %>%
    left_join(sel %>% select(.gid, consistent, n_criteria, always_top_level,
                             always_top_grad, always_top_dev, always_quadrant,
                             worst_rank_dev, median_exc_level, median_exc_grad,
                             min_q_level, min_q_gradient),
              by = ".gid") %>%
    mutate(raw_top      = rank_raw <= k,
           consistent   = coalesce(consistent, FALSE),
           selected     = raw_top & consistent) %>%
    arrange(desc(selected), rank_raw)
  list(raw = raw, long = long, consistency = sel, table = final,
       selected = final %>% filter(selected),
       k = k, voting = voting, shown = unique(long$path),
       domestic_mean = attr(raw, "domestic_mean"))
}

# =============================================================================
# 14. H3 RANK-STABILITY BOOTSTRAP (country-level, full-pipeline)
# =============================================================================
.resample_first_stage <- function(wells_by_unit, unit, cfg) {
  w <- wells_by_unit[[unit]]
  if (is.null(w)) return(NULL)
  i <- sample.int(nrow(w), nrow(w), replace = TRUE)
  o <- .first_stage_unit(w$GWSlp[i], w$dist_LB_km[i], cfg)
  o$r_use  <- if (cfg$first_stage_h2 == "spearman") o$r_spear else o$r
  o$r_pref <- o$r
  o
}

h3_bootstrap <- function(cfg = CFG) {
  wells <- fs_wells(cfg)
  base_d <- add_design_weights(build_segments(cfg), cfg)
  if (is.null(base_d)) stop("base segment table could not be built")
  base <- h3_analysis(base_d, cfg)
  if (is.null(base)) stop("base H3 analysis failed")

  wells_by_unit <- split(wells[, c("GWSlp", "dist_LB_km")], wells$unit_id)
  ps_cov <- attr(base_d, "ps_cov")
  units_by_cc <- split(seq_len(nrow(base_d)), base_d$CC)
  ccs <- names(units_by_cc)

  one_rep <- function(b) {
    if (cfg$h3_boot_level == "country") {
      draw <- sample(ccs, length(ccs), replace = TRUE)
      pieces <- lapply(seq_along(draw), function(k) {
        z <- base_d[units_by_cc[[draw[k]]], , drop = FALSE]
        z$CC_orig      <- z$CC
        z$unit_id_orig <- z$unit_id
        z$CC      <- paste0(z$CC, "#b", k)
        z$unit_id <- paste0(z$unit_id, "#b", k)
        z$aq_id   <- paste0(z$aq_id, "#b", k)
        z
      })
      x <- bind_rows(pieces)
    } else {
      idx <- sample(seq_len(nrow(base_d)), nrow(base_d), replace = TRUE)
      x <- base_d[idx, , drop = FALSE]
      x$CC_orig <- x$CC; x$unit_id_orig <- x$unit_id
      x$unit_id <- paste0(x$unit_id, "#", seq_len(nrow(x)))
    }
    if (!nrow(x)) return(NULL)
    if (isTRUE(cfg$h3_boot_resample_wells)) {
      fs <- map(x$unit_id_orig, ~ .resample_first_stage(wells_by_unit, .x, cfg))
      ok <- !vapply(fs, is.null, logical(1))
      x <- x[ok, , drop = FALSE]; fs <- fs[ok]
      if (!nrow(x)) return(NULL)
      x$beta_0 <- vapply(fs, function(o) o$beta_0, numeric(1))
      x$se_0   <- vapply(fs, function(o) o$se_0,   numeric(1))
      rr <- vapply(fs, function(o) o$r_use,  numeric(1))
      rp <- vapply(fs, function(o) o$r_pref, numeric(1))
      x$z      <- atanh(pmin(pmax(rr, -0.999), 0.999))
      x$z_pref <- atanh(pmin(pmax(rp, -0.999), 0.999))
      x <- x %>%
        mutate(se_0_con = se_0 * ratio_0,
               bv_0 = winsor(se_0_con^2, cfg$variance_winsor),
               bv_z = winsor(1 / pmax(n_eff - 3, 1), cfg$variance_winsor)) %>%
        filter(is.finite(z), is.finite(z_pref), is.finite(bv_0), bv_0 > 0)
      if (nrow(x) < 20 || length(unique(x$TBn)) < 2) return(NULL)
    }
    xw <- tryCatch(add_overlap_weights(x, ps_cov), error = function(e) NULL)
    if (is.null(xw)) return(NULL)
    xw$w_design <- xw$w_ato
    attr(xw, "ps_cov") <- ps_cov
    cfgb <- cfg_with(cfg, h3_variance = cfg$h3_boot_variance %||% cfg$h3_variance,
                     h3_blup_validate = FALSE)
    h <- tryCatch(h3_analysis(xw, cfgb), error = function(e) NULL)
    if (is.null(h)) return(NULL)
    h$scatter %>%
      transmute(gid = sub("#b?[0-9]+$", "", unit_id),
                rep = b, exc_level, exc_gradient,
                rank_level, rank_gradient, q_level, q_gradient)
  }

  reps <- par_map(seq_len(cfg$h3_boot_reps), one_rep,
                  label = "H3 bootstrap", cfg = cfg, seed = cfg$h3_boot_seed)
  n_ok <- sum(!vapply(reps, is.null, logical(1)))
  draws <- bind_rows(reps)
  if (!nrow(draws)) stop("no valid H3 bootstrap replicates")

  dd <- draws %>%
    group_by(gid, rep) %>%
    summarise(exc_level = mean(exc_level), exc_gradient = mean(exc_gradient),
              rank_level = min(rank_level), rank_gradient = min(rank_gradient),
              .groups = "drop")
  stab <- dd %>%
    group_by(gid) %>%
    summarise(
      n_present            = n(),
      median_rank_level    = median(rank_level),
      median_rank_gradient = median(rank_gradient),
      median_exc_level     = median(exc_level),
      lo_exc_level         = unname(quantile(exc_level, 0.05)),
      hi_exc_level         = unname(quantile(exc_level, 0.95)),
      median_exc_gradient  = median(exc_gradient),
      lo_exc_gradient      = unname(quantile(exc_gradient, 0.05)),
      hi_exc_gradient      = unname(quantile(exc_gradient, 0.95)),
      .groups = "drop") %>%
    mutate(
      top10_level_cond = vapply(gid, function(g)
        mean(dd$rank_level[dd$gid == g] <= 10), numeric(1)),
      top10_gradient_cond = vapply(gid, function(g)
        mean(dd$rank_gradient[dd$gid == g] <= 10), numeric(1)),
      top10_level_uncond = top10_level_cond * n_present / n_ok,
      top10_gradient_uncond = top10_gradient_cond * n_present / n_ok)

  list(base = base, draws = dd, stability = stab,
       n_valid = n_ok, n_requested = cfg$h3_boot_reps,
       level = cfg$h3_boot_level,
       scope_note = paste("Recomputes well-level coefficients, weighting,",
                          "second stage, shrinkage and ranking; holds the",
                          "pre-estimated Conley inflation ratios fixed."))
}

# =============================================================================
# 15. DESCRIPTIVE LOCALIZATION (Panel B and Block C)
# =============================================================================
# Within-TBA depletion profile by distance to the nearest international border.
# DESCRIPTIVE: transboundary segments only; no treated-control bin test is
# valid because domestic controls lack near-border support. Weighting: each
# well carries its segment's overlap weight divided by the segment's well
# count, so a segment's total influence equals its design weight. Uncertainty:
# country-cluster bootstrap.
profile_by_distance <- function(d, cfg = CFG, bins = cfg$profile_bins,
                                weighted = TRUE, centred = FALSE,
                                min_wells_bin = 1L, drop_cc_extra = NULL,
                                boot = cfg$profile_boot,
                                seed = cfg$wcb_seed) {
  wells <- fs_wells(cfg)
  seg <- d %>% select(unit_id, seg_CC = CC, TBn, w_ato, n_seg = n)
  w <- wells %>%
    inner_join(seg, by = "unit_id") %>%
    filter(TBn == 1)
  if (!is.null(drop_cc_extra)) w <- w %>% filter(!seg_CC %in% drop_cc_extra)
  w <- w %>%
    mutate(bin = cut(dist_LB_km, breaks = bins, right = FALSE,
                     include.lowest = TRUE),
           wt  = if (weighted) w_ato / pmax(n_seg, 1) else 1) %>%
    filter(!is.na(bin), is.finite(wt), wt > 0)
  if (!nrow(w)) return(NULL)
  if (centred)
    w <- w %>% group_by(unit_id) %>%
      mutate(GWSlp = GWSlp - wmean(GWSlp, wt)) %>% ungroup()

  point <- w %>%
    group_by(bin) %>%
    summarise(mean_depletion = wmean(GWSlp, wt),
              n_wells = n(), n_segments = n_distinct(unit_id),
              weight_contrib = sum(wt),
              .groups = "drop") %>%
    filter(n_wells >= min_wells_bin) %>%
    mutate(share_wells = n_wells / sum(n_wells),
           share_weight = weight_contrib / sum(weight_contrib))

  idx <- split(seq_len(nrow(w)), w$seg_CC)
  bootres <- with_seed(seed, {
    map_dfr(seq_len(boot), function(b) {
      draw <- sample(names(idx), length(idx), replace = TRUE)
      z <- w[unlist(idx[draw], use.names = FALSE), , drop = FALSE]
      z %>% group_by(bin) %>%
        summarise(m = wmean(GWSlp, wt), .groups = "drop") %>%
        mutate(rep = b)
    })
  })
  al <- (1 - cfg$ci_level) / 2
  ci <- bootres %>% group_by(bin) %>%
    summarise(lo = quantile(m, al, na.rm = TRUE),
              hi = quantile(m, 1 - al, na.rm = TRUE), .groups = "drop")
  point %>% left_join(ci, by = "bin")
}

# Compact localization summary: within vs beyond the influence zone, plus a
# continuous distance trend (weighted well-level regression of depletion on
# distance, within-segment centred), country-bootstrap uncertainty.
localization_summary <- function(d, cfg = CFG, weighted = TRUE,
                                 centred_trend = TRUE,
                                 drop_cc_extra = NULL, drop_aq_extra = NULL,
                                 boot = cfg$profile_boot,
                                 seed = cfg$wcb_seed + 3L,
                                 label = "preferred") {
  wells <- fs_wells(cfg)
  seg <- d %>% select(unit_id, seg_CC = CC, aq_id, TBn, w_ato, n_seg = n)
  w <- wells %>% inner_join(seg, by = "unit_id") %>% filter(TBn == 1)
  if (!is.null(drop_cc_extra)) w <- w %>% filter(!seg_CC %in% drop_cc_extra)
  if (!is.null(drop_aq_extra)) w <- w %>% filter(!aq_id %in% drop_aq_extra)
  w <- w %>% mutate(wt = if (weighted) w_ato / pmax(n_seg, 1) else 1) %>%
    filter(is.finite(wt), wt > 0, is.finite(GWSlp), is.finite(dist_LB_km))
  if (!nrow(w)) return(NULL)
  zone <- cfg$influence_zone_km

  stat_of <- function(z) {
    inn <- z$dist_LB_km <= zone
    ctr <- z %>% group_by(unit_id) %>%
      mutate(gw_c = GWSlp - wmean(GWSlp, wt),
             d_c  = dist_LB_km - wmean(dist_LB_km, wt)) %>% ungroup()
    tr <- if (centred_trend) {
      ok <- is.finite(ctr$gw_c) & is.finite(ctr$d_c)
      if (sum(ok) > 10 && sd(ctr$d_c[ok]) > 0)
        sum(ctr$wt[ok] * ctr$d_c[ok] * ctr$gw_c[ok]) /
          sum(ctr$wt[ok] * ctr$d_c[ok]^2)
      else NA_real_
    } else {
      ok <- is.finite(z$GWSlp) & is.finite(z$dist_LB_km)
      xz <- z$dist_LB_km[ok] - wmean(z$dist_LB_km[ok], z$wt[ok])
      yz <- z$GWSlp[ok]
      if (sum(ok) > 10 && sd(xz) > 0)
        sum(z$wt[ok] * xz * yz) / sum(z$wt[ok] * xz^2)
      else NA_real_
    }
    c(within = wmean(z$GWSlp[inn], z$wt[inn]),
      beyond = wmean(z$GWSlp[!inn], z$wt[!inn]),
      trend  = tr)
  }
  obs <- stat_of(w)
  idx <- split(seq_len(nrow(w)), w$seg_CC)
  bs <- with_seed(seed, {
    t(replicate(boot, {
      draw <- sample(names(idx), length(idx), replace = TRUE)
      stat_of(w[unlist(idx[draw], use.names = FALSE), , drop = FALSE])
    }))
  })
  al <- (1 - cfg$ci_level) / 2
  qs <- apply(bs, 2, quantile, probs = c(al, 1 - al), na.rm = TRUE)
  tibble(
    specification = label,
    mean_within_100km  = obs["within"],
    within_lo = qs[1, "within"], within_hi = qs[2, "within"],
    mean_beyond_100km  = obs["beyond"],
    beyond_lo = qs[1, "beyond"], beyond_hi = qs[2, "beyond"],
    contrast_within_minus_beyond = obs["within"] - obs["beyond"],
    contrast_lo = quantile(bs[, "within"] - bs[, "beyond"], al, na.rm = TRUE),
    contrast_hi = quantile(bs[, "within"] - bs[, "beyond"], 1 - al, na.rm = TRUE),
    share_wells_within_100km = sum(w$dist_LB_km <= zone) / nrow(w),
    weighted_contrib_within_100km =
      sum(w$wt[w$dist_LB_km <= zone]) / sum(w$wt),
    n_segments_within = n_distinct(w$unit_id[w$dist_LB_km <= zone]),
    n_segments_beyond = n_distinct(w$unit_id[w$dist_LB_km > zone]),
    trend_mm_yr_per_km = obs["trend"],
    trend_lo = qs[1, "trend"], trend_hi = qs[2, "trend"],
    n_wells = nrow(w), n_countries = n_distinct(w$seg_CC),
    inference = "descriptive; country-cluster bootstrap interval")
}

# Cumulative distance-window gradient summary (within-TBA; SI block C).
window_summary <- function(d, cfg = CFG) {
  wells <- fs_wells(cfg)
  seg <- d %>% select(unit_id, seg_CC = CC, TBn, w_ato)
  w <- wells %>% inner_join(seg, by = "unit_id") %>% filter(TBn == 1)
  map_dfr(cfg$gradient_windows, function(W) {
    z <- w %>% filter(dist_LB_km <= W) %>%
      group_by(unit_id, seg_CC, w_ato) %>%
      summarise(n_in = n(), sd_d = sd(dist_LB_km),
                level = mean(GWSlp),
                r = suppressWarnings(cor(GWSlp, dist_LB_km)),
                .groups = "drop") %>%
      filter(n_in >= cfg$gradient_min_wells)
    if (!nrow(z)) return(NULL)
    z <- z %>% mutate(zz = if_else(is.finite(sd_d) & sd_d > 0 & is.finite(r),
                                   atanh(pmin(pmax(r, -0.999), 0.999)), NA_real_))
    tibble(window_km = W,
           n_segments = nrow(z),
           mean_level = wmean(z$level, z$w_ato),
           mean_gradient_z = wmean(z$zz, z$w_ato),
           n_segments_gradient = sum(is.finite(z$zz)))
  })
}
