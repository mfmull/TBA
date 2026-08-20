# =============================================================================
# 53_run_main.R -- PREFERRED CLASSIFICATION AND FROZEN BENCHMARK
#
# Classifies every dyad under the published settings, computes the Fig. 4B/C
# statistics and the two mosaic tables, verifies the frozen benchmark, and
# writes derived/cache/main_objects.rds for 55-57.
#
# Canonical outputs are written by 55-57; this script writes none.
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

say("\n================ PREFERRED CLASSIFICATION =====================\n")
sayf("Rule: %s, near %s, far %s, eps %g (= %g ha), threshold %g\n",
     cfg$land_var, cfg$near_zone, cfg$far_zone, cfg$eps, cfg$eps * 1e4,
     cfg$thresh)
sayf("Stage: dyad-border zones %s, audit corrections %s.\n",
     if (isTRUE(cfg$dyad_zones)) "ON" else "off",
     if (is.null(cfg$manual_class)) "off"
     else sprintf("ON (%d rows)", nrow(cfg$manual_class)))

d  <- load_dyads(cfg)
tp <- build_typology(cfg, d)

# ---- assertions: the sample --------------------------------------------------
n_countries <- length(unique(c(tp$CC_1, tp$CC_2)))
n_aquifers  <- length(unique(tp$code))
sayf("%d dyads, %d countries, %d aquifers.\n",
     nrow(tp), n_countries, n_aquifers)
stopifnot(nrow(tp) == cfg$reference_n_dyads,
          n_countries == cfg$reference_n_countries,
          n_aquifers  == cfg$reference_n_aquifers)
# The rule must be exhaustive: an unclassified dyad is a bug, not a category.
stopifnot(!any(is.na(tp$class)), !any(is.na(tp$classF)))

# ---- assertion: the frozen benchmark ----------------------------------------
cur <- as.integer(table(tp$class)[cfg$class_levels])
fut <- as.integer(table(tp$classF)[cfg$class_levels])
names(cur) <- names(fut) <- cfg$class_levels
say("Current: ", paste(sprintf("%s=%d", names(cur), cur), collapse = "  "), "\n")
say("Future : ", paste(sprintf("%s=%d", names(fut), fut), collapse = "  "), "\n")
# Compare UNNAMED: `cur` and `fut` carry names for the printed line above,
# while as.integer() strips them from the reference, and identical() compares
# attributes as well as values -- so a named-vs-unnamed pair of otherwise equal
# vectors compares FALSE. The message also reports what was observed, not only
# what was expected: a failure that prints only the expectation is unreadable
# when the observed values are sitting right above it and look correct.
ref_cur <- as.integer(cfg$reference_class[cfg$class_levels])
ref_fut <- as.integer(cfg$reference_classF[cfg$class_levels])
if (!identical(unname(cur), ref_cur) || !identical(unname(fut), ref_fut))
  stop("The classification does not reproduce the frozen benchmark.\n",
       "  current  observed ", paste(unname(cur), collapse = "/"),
       "   expected ", paste(ref_cur, collapse = "/"), "\n",
       "  future   observed ", paste(unname(fut), collapse = "/"),
       "   expected ", paste(ref_fut, collapse = "/"))
say("Finalised classification reproduces its frozen benchmark.\n")

# ---- assertion: the PUBLISHED classification is still reproducible ----------
# The corrections are adopted, not baked in. Rebuilding stage 1 from the same
# code path and asserting the published counts is what keeps that true: if a
# later edit made the segment-based classification unreachable, or quietly
# changed it, this stops the run. It costs one extra pass over 568 rows.
cfg_pub <- cfg_with(cfg, dyad_zones = FALSE, manual_class = NULL)
tp_pub  <- build_typology(cfg_pub, d)
cur_pub <- as.integer(table(tp_pub$class)[cfg$class_levels])
fut_pub <- as.integer(table(tp_pub$classF)[cfg$class_levels])
ref_cur_pub <- as.integer(cfg$reference_class_published[cfg$class_levels])
ref_fut_pub <- as.integer(cfg$reference_classF_published[cfg$class_levels])
if (!identical(cur_pub, ref_cur_pub) || !identical(fut_pub, ref_fut_pub))
  stop("The published (stage 1) classification is no longer reproducible.\n",
       "  current  observed ", paste(cur_pub, collapse = "/"),
       "   expected ", paste(ref_cur_pub, collapse = "/"), "\n",
       "  future   observed ", paste(fut_pub, collapse = "/"),
       "   expected ", paste(ref_fut_pub, collapse = "/"))
say("Published stage-1 classification still reproduces 85/118/207/113/45.\n")
sayf("Correction moved %d of %d dyads on class and %d on classF.\n",
     sum(as.character(tp$class)  != as.character(tp_pub$class)),
     nrow(tp),
     sum(as.character(tp$classF) != as.character(tp_pub$classF)))

# ---- Fig. 4B/C statistics ----------------------------------------------------
say("\n---- Fig. 4B/C: group means, Tukey groupings ----\n")
sayf("Tukey at %.0f%% (letters at alpha = %.2f).\n",
     100 * cfg$tukey_conf, cfg$tukey_alpha)
tuk <- tukey_table(tp, cfg)
boot <- cluster_bootstrap(tp, cfg)
if (!is.null(boot)) {
  tuk <- tuk %>% left_join(boot, by = c("Variable", "class"))
  sayf("Country-clustered bootstrap: %d replicates over %d countries.\n",
       cfg$cluster_boot_reps, n_countries)
}
if (isTRUE(cfg$report_kruskal)) {
  kw <- tuk %>% distinct(Variable, kruskal_p)
  for (i in seq_len(nrow(kw)))
    sayf("  Kruskal-Wallis %s: p = %.3g (distribution-free check)\n",
         kw$Variable[i], kw$kruskal_p[i])
}
print(as.data.frame(tuk %>% select(Variable, class, n, mean, se, letters)),
      row.names = FALSE, digits = 3)

# ---- mosaic tables -----------------------------------------------------------
say("\n---- contingency tables ----\n")
tp <- tp %>% mutate(CoopType = factor(ArrangementScore),
                    dist_class = distance_class(d, cfg))
mosaics <- list(
  cooperation = list(var = "CoopType",
                     title = "Cooperation arrangement level"),
  distance    = list(var = "dist_class",
                     title = "Distance class of aquifer extent"),
  river       = list(var = "has_river_border",
                     title = "River-defined border"))
for (nm in names(mosaics)) {
  # Name the dimnames. table(tp$class, tp[[var]]) leaves them EMPTY, because
  # the second argument is a computed expression rather than a symbol -- and
  # vcd's labelling draws the variable names from dimnames(), so an unnamed
  # table makes mosaic() fail at draw time with nothing useful to show for it.
  tab <- table(tp$class, tp[[mosaics[[nm]]$var]])
  names(dimnames(tab)) <- c("class", mosaics[[nm]]$var)
  mosaics[[nm]]$table <- tab
  mosaics[[nm]]$test  <- contingency_test(tab, cfg)
  sayf("  %-12s p = %.4g  (%s; %s)\n", nm, mosaics[[nm]]$test$p_value,
       mosaics[[nm]]$test$method, mosaics[[nm]]$test$note)
}

main <- list(typology = tp, tukey = tuk, boot = boot, mosaics = mosaics,
             counts = list(current = cur, future = fut,
                           n_dyads = nrow(tp), n_countries = n_countries,
                           n_aquifers = n_aquifers),
             cfg = cfg)
saveRDS(main, file.path(cfg$cache_dir, "main_objects.rds"))
sayf("\n53_run_main.R done in %.2f min.\n",
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
