# =============================================================================
# missing_dyad_codes.R -- WHICH DYAD-ZONE ROWS ARE STILL MISSING
#
#   Rscript 4_typology/gee/missing_dyad_codes.R
#
# Builds the full set of (code, CC, CC_other) keys the correction export
# should contain, from outTBA_Typol.csv -- the same segment table GEE reads --
# and diffs it against whatever dyadB*_*.csv parts are currently in
# 1_data/geeOut/.
#
# It checks KEYS, not filenames, so it catches a chunk that ran but was never
# downloaded just as well as one that never finished, and it does not depend
# on how the export happened to be chunked.
#
# Prints a ready-to-paste ONLY_CODES line for
# gee/dyad_border_corrections.js. Re-running by aquifer code is stable;
# re-running "chunk 05" is not, because chunk numbering changed when the
# partition moved from list position to aquifer.
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.here <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
root <- normalizePath(file.path(.here, ".."), mustWork = FALSE)
if (!file.exists(file.path(root, "51_config.R"))) root <- getwd()

gee_dir <- file.path(root, "1_data", "geeOut")
seg_file <- file.path(gee_dir, "outTBA_Typol.csv")
if (!file.exists(seg_file))
  stop("Cannot find ", seg_file, ". Run this from 4_typology/gee/.")

MIN_RIPARIANS <- 3L
WIDTHS <- list(
  list(tag = "10k",  pattern = "^dyadB10k(_[A-Za-z0-9]+)*[.]csv$"),
  list(tag = "100k", pattern = "^dyadB100k(_[A-Za-z0-9]+)*[.]csv$"))

# ---- expected keys -----------------------------------------------------------
seg <- read.csv(seg_file, stringsAsFactors = FALSE)
seg <- unique(seg[, c("code", "CC")])
per <- split(seg$CC, seg$code)
per <- per[vapply(per, length, 1L) >= MIN_RIPARIANS]

expected <- do.call(rbind, lapply(names(per), function(cd) {
  cc <- sort(unique(per[[cd]]))
  g  <- expand.grid(CC = cc, CC_other = cc, stringsAsFactors = FALSE)
  g  <- g[g$CC != g$CC_other, , drop = FALSE]
  data.frame(code = cd, g, stringsAsFactors = FALSE)
}))
exp_key <- paste(expected$code, expected$CC, expected$CC_other, sep = "|")

cat(sprintf("Expected: %d dyad-side rows across %d multi-riparian aquifers.\n",
            length(exp_key), length(per)))

# ---- what is on disk ---------------------------------------------------------
missing_codes <- character(0)
for (w in WIDTHS) {
  parts <- list.files(gee_dir, pattern = w$pattern, full.names = TRUE)
  if (!length(parts)) {
    cat(sprintf("\n%-5s no parts on disk -- every aquifer is missing.\n", w$tag))
    missing_codes <- union(missing_codes, names(per))
    next
  }
  got <- do.call(rbind, lapply(parts, function(p)
    read.csv(p, stringsAsFactors = FALSE)[, c("code", "CC", "CC_other")]))
  got_key <- paste(got$code, got$CC, got$CC_other, sep = "|")

  dup <- got_key[duplicated(got_key)]
  miss <- setdiff(exp_key, got_key)
  extra <- setdiff(got_key, exp_key)

  cat(sprintf("\n%-5s %d part(s), %d row(s): %d of %d keys present.\n",
              w$tag, length(parts), nrow(got),
              length(unique(got_key)), length(exp_key)))
  if (length(dup))
    cat(sprintf("  %d DUPLICATE key(s) -- a part is on disk twice. ",
                length(dup)),
        "apply_dyad_zones() will stop on this; delete the extra file.\n")
  if (length(extra))
    cat(sprintf("  %d row(s) not in the expected set -- check MIN_RIPARIANS.\n",
                length(extra)))
  if (length(miss)) {
    mc <- sort(unique(sub("[|].*$", "", miss)))
    cat(sprintf("  %d missing row(s) across %d aquifer(s).\n",
                length(miss), length(mc)))
    missing_codes <- union(missing_codes, mc)
  } else {
    cat("  complete.\n")
  }
}

# ---- the answer --------------------------------------------------------------
cat("\n")
if (!length(missing_codes)) {
  cat("Nothing missing. Both widths are complete; run 59_dyad_correction.R.\n")
} else {
  missing_codes <- sort(missing_codes)
  # A partial re-run must not overwrite the parts already downloaded, so the
  # RUN_TAG reminder travels with the codes.
  cat(sprintf(paste0("Paste into gee/dyad_border_corrections.js, set RUN_TAG ",
                     "to something new (e.g. 'r1'), and drop NGROUPS to 1-2:\n\n")))
  cat("var ONLY_CODES = ['", paste(missing_codes, collapse = "','"), "'];\n",
      sep = "")
  cat(sprintf("\n(%d aquifer(s). Both widths are re-run for these codes; if ",
              length(missing_codes)),
      "only one width was\nincomplete, comment the other out of WIDTHS.)\n",
      sep = "")
}
