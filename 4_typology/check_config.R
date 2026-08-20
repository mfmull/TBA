# Verification harness. 52_core.R needs dplyr for the plumbing, but the RULE
# itself -- classify_dyads() -- and 51_config.R are pure base R. This sources
# the real config, extracts the real classify_dyads() and zone_col() from
# 52_core.R without loading the packages, reimplements only the dyad-zone join
# in base R, and checks every frozen number in CFG.
#
#   Rscript check_config.R
#
# It exercises the actual AUDIT_FINAL table and the actual manual_class guard,
# so a stale `from` or a mistyped class stops it exactly as the pipeline would.

here <- dirname(normalizePath(sub("^--file=", "",
  commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
source(file.path(here, "51_config.R"))
cfg <- CFG

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a
ex <- parse(file.path(here, "52_core.R"))
for (e in ex) {
  if (is.call(e) && identical(as.character(e[[1]]), "<-") &&
      is.name(e[[2]]) &&
      any(as.character(e[[2]]) == c("classify_dyads", "zone_col",
                                    "require_cols", "cfg_with")))
    eval(e, envir = globalenv())
}

d <- read.csv(cfg$dyad_file, stringsAsFactors = FALSE)

# ---- the dyad-zone overlay, in base R ---------------------------------------
read_zone <- function(src) {
  parts <- sort(Sys.glob(file.path(cfg$dyad_zone_dir,
                                   cfg$dyad_zone_files[[src]])))
  stopifnot(length(parts) > 0)
  need <- c("code", "CC", "CC_other")
  z <- do.call(rbind, lapply(parts, function(p) {
    x <- read.csv(p, stringsAsFactors = FALSE)
    x[, c(need, intersect(cfg$dyad_zone_vars, names(x))), drop = FALSE]
  }))
  k <- paste(z$code, z$CC, z$CC_other, sep = "\r")
  v <- intersect(cfg$dyad_zone_vars, names(z))
  for (nm in v) {
    lo <- tapply(z[[nm]], k, min); hi <- tapply(z[[nm]], k, max)
    den <- pmax(abs(hi), .Machine$double.eps)
    stopifnot(all((hi - lo) / den <= cfg$dyad_zone_dup_tol))
  }
  z[!duplicated(k), , drop = FALSE]
}

overlay <- function(d) {
  covered <- rep(FALSE, nrow(d))
  for (src in names(cfg$dyad_zone_files)) {
    z <- read_zone(src)
    v <- intersect(cfg$dyad_zone_vars, names(z))
    for (side in 1:2) {
      other <- if (side == 1L) 2L else 1L
      k  <- paste(d$code, d[[paste0("CC_", side)]], d[[paste0("CC_", other)]],
                  sep = "\r")
      zk <- paste(z$code, z$CC, z$CC_other, sep = "\r")
      i  <- match(k, zk)
      h  <- !is.na(i)
      for (nm in v) {
        oldc <- zone_col(nm, src, side)
        if (!oldc %in% names(d)) next
        d[[oldc]][h] <- z[[nm]][i[h]]
        covered <- covered | h
      }
    }
  }
  attr(d, "covered") <- covered
  d
}

dc <- overlay(d)
cnt <- function(x) as.integer(table(factor(x, cfg$class_levels)))

build <- function(cfgx) {
  dd <- if (isTRUE(cfgx$dyad_zones)) dc else d
  list(cur = cnt(classify_dyads(dd, cfgx, FALSE)),
       fut = cnt(classify_dyads(dd, cfgx, TRUE)))
}

fail <- 0L
chk <- function(what, got, want) {
  ok <- identical(as.integer(got), as.integer(want))
  cat(sprintf("%-22s %-24s %s\n", what, paste(got, collapse = ", "),
              if (ok) "OK" else paste("MISMATCH, expected",
                                      paste(want, collapse = ", "))))
  if (!ok) fail <<- fail + 1L
}

cat("dyads", nrow(d), "| measured against own border",
    sum(attr(dc, "covered")), "\n\n")

cat("---- preferred (finalised, stage 3) ----\n")
p <- build(cfg)
chk("class",  p$cur, cfg$reference_class[cfg$class_levels])
chk("classF", p$fut, cfg$reference_classF[cfg$class_levels])

cat("\n---- published (stage 1) ----\n")
q <- build(cfg_with(cfg, dyad_zones = FALSE, manual_class = NULL))
chk("class",  q$cur, cfg$reference_class_published[cfg$class_levels])
chk("classF", q$fut, cfg$reference_classF_published[cfg$class_levels])

cat("\n---- registry ----\n")
for (nm in names(cfg$variants)) {
  v  <- cfg$variants[[nm]]
  ov <- v[setdiff(names(v), c("label", "si"))]
  r  <- build(do.call(cfg_with, c(list(cfg), ov)))
  ref <- cfg$reference_variants[[nm]]
  if (is.null(ref)) { cat(nm, ": no frozen counts\n"); next }
  chk(paste0(nm, " class"),  r$cur, ref$class)
  chk(paste0(nm, " classF"), r$fut, ref$classF)
}

cat("\n---- audit tables ----\n")
cat("AUDIT_FINAL   ", nrow(AUDIT_FINAL), "rows;",
    sum(AUDIT_FINAL$tier == 2), "tier 2\n")
cat("AUDIT_FINAL_T1", nrow(AUDIT_FINAL_T1), "rows\n")
stopifnot(!anyDuplicated(paste(AUDIT_FINAL$code, AUDIT_FINAL$CC_1,
                               AUDIT_FINAL$CC_2)))
stopifnot(all(unlist(AUDIT_FINAL[c("from", "to", "fromF", "toF")])
              %in% cfg$class_levels))
cat("keys unique, all classes valid.\n")

cat(if (fail) sprintf("\n*** %d CHECK(S) FAILED\n", fail)
    else "\nAll checks passed.\n")
quit(status = if (fail) 1L else 0L)
