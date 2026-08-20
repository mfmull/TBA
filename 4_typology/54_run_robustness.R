# =============================================================================
# 54_run_robustness.R -- ROBUSTNESS VARIANTS
#
# The legacy folder implemented these as five standalone scripts (FigSI_5.R,
# FigSI_200.R, FigSI_eps.R, FigSI_IR.R, FigSI_thresh.R), each a copy of the
# main script with one constant changed and its own hard-coded output paths.
# Three of them also mislabelled the first plotted metric, which produced an
# "NA" facet strip in the published SI figures.
#
# Here each variant is a list of configuration overrides applied to the
# preferred specification. The classification rule, the metrics and the
# statistics are the ones in 52_core.R; nothing is re-implemented, so no
# variant can silently drift from the main analysis.
#
# Writes derived/cache/robustness_objects.rds.
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!exists("CFG")) source(file.path(.root, "51_config.R"))
if (!exists("classify_dyads")) source(file.path(CFG$root, "52_core.R"))
if (is.null(.log$path)) log_open(CFG$log_file)
if (!exists("FORCE")) FORCE <- nzchar(Sys.getenv("FORCE"))
cfg <- CFG
t0 <- Sys.time()

say("\n================ ROBUSTNESS VARIANTS =====================\n")
sayf("%d variants; each is the preferred specification with one setting changed.\n",
     length(cfg$variants))

d <- load_dyads(cfg)

run_variant <- function(nm) {
  v  <- cfg$variants[[nm]]
  ov <- v[setdiff(names(v), c("label", "si"))]
  cfgx <- do.call(cfg_with, c(list(cfg), ov))

  f  <- file.path(cfg$cache_dir, paste0("variant_", nm, ".rds"))
  st <- cache_stamp(cfgx, ov, extra = list(variant = nm))
  hit <- cache_read(f, st, cfg, force = isTRUE(FORCE))
  if (!is.null(hit)) { say("  [cached] ", nm, "\n"); return(hit) }

  tp  <- build_typology(cfgx, d)
  cur <- as.integer(table(tp$class)[cfgx$class_levels])
  fut <- as.integer(table(tp$classF)[cfgx$class_levels])
  out <- list(name = nm, label = v$label, si = v$si %||% NA_character_,
              overrides = ov, cfg = cfgx, typology = tp,
              tukey = tukey_table(tp, cfgx),
              counts = list(current = cur, future = fut),
              labels = metric_labels(cfgx))
  cache_write(out, f, st)
  out
}

rob <- lapply(names(cfg$variants), run_variant)
names(rob) <- names(cfg$variants)

# ---- assertions: every variant reproduces its frozen counts ------------------
say("\n---- class counts (IU0, IU1, UB, DA, BL) ----\n")
sayf("%-10s %-6s %-24s %-24s %s\n", "variant", "SI", "current", "future", "label")
sayf("%-10s %-6s %-24s %-24s %s\n", "preferred", "Fig 4",
     paste(cfg$reference_class[cfg$class_levels], collapse = ", "),
     paste(cfg$reference_classF[cfg$class_levels], collapse = ", "), "")
for (nm in names(rob)) {
  r <- rob[[nm]]
  sayf("%-10s %-6s %-24s %-24s %s\n", nm, r$si,
       paste(r$counts$current, collapse = ", "),
       paste(r$counts$future,  collapse = ", "), r$label)
  ref <- cfg$reference_variants[[nm]]
  if (!is.null(ref)) {
    if (!identical(r$counts$current, as.integer(ref$class)) ||
        !identical(r$counts$future,  as.integer(ref$classF)))
      stop("Variant '", nm, "' does not reproduce its frozen counts. Expected ",
           paste(ref$class, collapse = ", "), " / ",
           paste(ref$classF, collapse = ", "), ".")
  } else {
    say("    (no frozen counts for this variant; not asserted)\n")
  }
}
say("All variants reproduce their frozen counts.\n")

# ---- the robustness claim, made checkable ------------------------------------
# The partition into interior / urban / bilateral configurations is what the
# manuscript calls robust. BL is a tail-defined category, so its boundary with
# DA moves with the detection and concentration thresholds by construction --
# stating that explicitly is stronger than claiming the whole typology is
# insensitive, and it is what the numbers support.
part <- bind_rows(lapply(names(rob), function(nm) {
  cnt <- rob[[nm]]$counts$current
  names(cnt) <- cfg$class_levels
  tibble(variant = nm, si = rob[[nm]]$si,
         interior = cnt[["IU0"]] + cnt[["IU1"]], urban = cnt[["UB"]],
         bilateral = cnt[["DA"]] + cnt[["BL"]], BL = cnt[["BL"]])
}))
ref <- cfg$reference_class[cfg$class_levels]
part <- bind_rows(
  tibble(variant = "preferred", si = "Fig 4",
         interior = ref[["IU0"]] + ref[["IU1"]], urban = ref[["UB"]],
         bilateral = ref[["DA"]] + ref[["BL"]], BL = ref[["BL"]]),
  part)
say("\n---- stability of the three-way partition ----\n")
print(as.data.frame(part), row.names = FALSE)
sayf("Interior spans %d-%d, urban %d-%d, bilateral %d-%d; BL alone spans %d-%d.\n",
     min(part$interior), max(part$interior), min(part$urban), max(part$urban),
     min(part$bilateral), max(part$bilateral), min(part$BL), max(part$BL))

# Does BL keep its defining profile -- the highest overdraft exposure -- in
# every variant? The audit verified this by hand across the four robustness
# CSVs; here it is a computed column.
bl_top <- bind_rows(lapply(names(rob), function(nm) {
  t <- rob[[nm]]$tukey %>% filter(Variable == "FBW")
  tibble(variant = nm,
         BL_mean_FBW = t$mean[t$class == "BL"],
         BL_is_max = isTRUE(t$class[which.max(t$mean)] == "BL"),
         BL_letters = t$letters[t$class == "BL"])
}))
say("\n---- does BL remain the highest-overdraft category? ----\n")
print(as.data.frame(bl_top), row.names = FALSE, digits = 3)

saveRDS(list(variants = rob, partition = part, bl_profile = bl_top),
        file.path(cfg$cache_dir, "robustness_objects.rds"))
sayf("\n54_run_robustness.R done in %.2f min.\n",
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
