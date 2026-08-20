# =============================================================================
# 33_run_main.R -- PREFERRED SPECIFICATION AND FROZEN BENCHMARK
#
# Fits the preferred specification (ATO x precision, country + segment random
# intercepts), the country-fixed-effects identification check, the H3 panel
# and the coded aquifer selection; validates the design (balance, ESS); and
# freezes the result so every later script has something to be compared with.
#
# Writes (derived/cache): main_objects.rds, main_reference.rds
# Canonical outputs are written by 35-37; this script writes none.
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!exists("CFG")) source(file.path(.root, "31_config.R"))
if (!exists("run_specification")) source(file.path(CFG$root, "32_core.R"))
log_open(CFG$log_file)
assert_mde(CFG)

say("\n================ PREFERRED SPECIFICATION =====================\n")
t0 <- Sys.time()

pref <- run_specification("pref", "Preferred specification (ATO)", "Design",
                          run_h3 = TRUE, with_att = TRUE, wcb_ci = TRUE)
if (!isTRUE(pref$ok))
  stop("The preferred specification failed: ", pref$diagnostics$warnings)

d   <- pref$data
s1  <- pref$s1
s2  <- pref$s2
cfg <- pref$cfg

# ---- assertions: samples, directions, placebo -------------------------------
stopifnot(identical(sort(unique(d$TBn)), c(0L, 1L)))
stopifnot(all(d$n >= cfg$min_wells))
stopifnot(!any(d$CC %in% cfg$drop_cc))
# H1 and H2 must use the SAME segments (common-outcome filter):
stopifnot(nrow(s1$data) == nrow(s2$data))
# Directional alternatives as pre-stated:
stopifnot(identical(cfg$directional_alt[["H1"]], "greater"),
          identical(cfg$directional_alt[["H2"]], "less"))
# Within-distance shares in [0, 1]:
for (km in cfg$coverage_km) {
  v <- d[[paste0("share_within_", km, "km")]]
  stopifnot(all(v >= 0 & v <= 1, na.rm = TRUE))
}

# ---- frozen benchmark -------------------------------------------------------
if (!is.finite(cfg$reference_h1) || !is.finite(cfg$reference_h2)) {
  if (!isTRUE(cfg$allow_unlocked_reference))
    stop("CFG$reference_h1/h2 are NA and allow_unlocked_reference is FALSE. ",
         "Paste the frozen benchmark into 31_config.R (printed below) or set ",
         "allow_unlocked_reference = TRUE for an exploratory pass.")
  say("*** UNLOCKED RUN: no frozen benchmark in 31_config.R yet. ***\n")
}
check_against <- function(label, target, value, tol = cfg$reproduce_tol) {
  if (is.null(target) || !all(is.finite(target))) return(invisible(NULL))
  if (any(abs(value - target) > tol * pmax(1, abs(target))))
    stop(sprintf("PREFERRED SPECIFICATION NOT REPRODUCED. %s: got %s, expected %s.",
                 label, paste(signif(value, 10), collapse = ", "),
                 paste(signif(target, 10), collapse = ", ")))
  say(sprintf("  %s reproduces the frozen value.\n", label))
}
check_against("H1", cfg$reference_h1, pref$h1$estimate)
check_against("H2", cfg$reference_h2, pref$h2$estimate)
check_against("H1 90% CI", cfg$reference_h1_ci,
              c(pref$h1$ci90_low, pref$h1$ci90_high), cfg$reproduce_tol_ci)
check_against("H2 90% CI", cfg$reference_h2_ci,
              c(pref$h2$ci90_low, pref$h2$ci90_high), cfg$reproduce_tol_ci)
check_against("segments",  cfg$reference_n_segments,  pref$sample$n_segments, 0)
check_against("treated",   cfg$reference_n_treated,   pref$sample$n_treated, 0)
check_against("countries", cfg$reference_n_countries, pref$sample$n_countries, 0)

say("\nPaste into 31_config.R to lock (if not already locked):\n")
sayf("  reference_h1 = %.10g,  reference_h2 = %.10g,\n",
     pref$h1$estimate, pref$h2$estimate)
sayf("  reference_h1_ci = c(%.10g, %.10g),\n", pref$h1$ci90_low, pref$h1$ci90_high)
sayf("  reference_h2_ci = c(%.10g, %.10g),\n", pref$h2$ci90_low, pref$h2$ci90_high)
sayf("  reference_n_segments = %dL, reference_n_treated = %dL, reference_n_countries = %dL\n",
     pref$sample$n_segments, pref$sample$n_treated, pref$sample$n_countries)

if (!isTRUE(pref$diagnostics$converged))
  stop("tau2 fixed point not reached; raise CFG$max_iter.")

# ---- design validation ------------------------------------------------------
bal  <- balance_table(d, cfg)
bsum <- balance_summary(d, bal)
say("\nBalance (SMD, treated SD denominator):\n")
print(as.data.frame(bal), digits = 3, row.names = FALSE)
if (bsum$max_abs_smd_ato > cfg$balance_tol_substantive)
  stop("Substantive balance failure under ATO: max |SMD| = ",
       signif(bsum$max_abs_smd_ato, 3))

# ---- headline ----------------------------------------------------------------
mde80 <- pref$h1$mde80; mde90 <- pref$h1$mde90
say("\n---------------- HEADLINE ----------------\n")
sayf("H1  %+.2f mm/yr, 90%% CI [%+.2f, %+.2f]; one-sided p = %s (%s)\n",
     pref$h1$estimate, pref$h1$ci90_low, pref$h1$ci90_high,
     pstar(pref$h1$p_one), pref$h1$p_one_source)
sayf("    80%% MDE = %.1f mm/yr | 90%% MDE = %.1f mm/yr (SE source: %s)\n",
     mde80, mde90, pref$h1$se_source)
sayf("H2  %+.4f Fisher z, 90%% CI [%+.4f, %+.4f]; one-sided p = %s (%s)\n",
     pref$h2$estimate, pref$h2$ci90_low, pref$h2$ci90_high,
     pstar(pref$h2$p_one), pref$h2$p_one_source)
sayf("    placebo (intercept): %+.4f [%+.4f, %+.4f], two-sided p = %s\n",
     pref$h2$baseline, pref$h2$baseline_lo, pref$h2$baseline_hi,
     pstar(pref$h2$baseline_p))

# ---- identification check: country fixed effects ----------------------------
cfe <- run_specification("cfe", "Country fixed effects", "Design",
                         country_effect = "fixed")

# ---- H3 panel and coded selection -------------------------------------------
h3 <- pref$h3
if (is.null(h3)) stop("H3 could not be estimated under the preferred path.")
if (any(h3$blup_validation$ok %in% FALSE))
  stop("Closed-form H3 BLUPs disagree with the direct GLS BLUP.")
ident <- identify_segments(d, cfg)
if (is.null(ident)) stop("H3 selection could not be computed.")
sayf("\nH3: %d transboundary segments; %d selected by the coded rule (%s).\n",
     nrow(h3$scatter), nrow(ident$selected),
     paste(ident$selected$label_txt, collapse = "; "))

# ---- Panel B profile and localization summary -------------------------------
prof <- profile_by_distance(d, cfg)
loc  <- localization_summary(d, cfg)

# ---- descriptive distance-window curve (SI block C) -------------------------
wsum <- window_summary(d, cfg)

# ---- persist -----------------------------------------------------------------
ref <- list(h1 = pref$h1$estimate, h2 = pref$h2$estimate,
            h1_ci = c(pref$h1$ci90_low, pref$h1$ci90_high),
            h2_ci = c(pref$h2$ci90_low, pref$h2$ci90_high),
            n_segments = pref$sample$n_segments,
            n_treated = pref$sample$n_treated,
            n_countries = pref$sample$n_countries,
            locked = is.finite(cfg$reference_h1) && is.finite(cfg$reference_h2),
            stamp = format(Sys.time()))
cache_write(ref, file.path(cfg$cache_dir, "main_reference.rds"),
            cache_stamp(cfg, list(spec = "pref")))

main <- list(pref = pref, cfe = cfe, h3 = h3, ident = ident,
             balance = bal, balance_summary = bsum,
             mde80 = mde80, mde90 = mde90,
             profile = prof, localization = loc, windows = wsum,
             reference = ref)
cache_write(main, file.path(cfg$cache_dir, "main_objects.rds"),
            cache_stamp(cfg, list(spec = "pref"), extra = list(object = "main")))
sayf("\n33_run_main.R done in %.1f min.\n",
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
