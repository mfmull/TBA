# =============================================================================
# 49_asymmetry.R -- SUPPLEMENTAL: ECONOMIC-ASYMMETRY HETEROGENEITY
#
#   Rscript 49_asymmetry.R          (or run via 48_run_all.R)
#
# Theory (SI Section S1 / Mullen et al. 2022): transboundary overdraft is
# enhanced when cross-border asymmetries in the energy (pumping) price and in
# the marginal value of water are BOTH large AND aligned -- i.e., the side
# facing more expensive energy also faces a higher water value. Misaligned
# asymmetries dampen the distortion.
#
# Empirical translation ("is there a signal?"):
#   - Energy price: national gasoline pump price (USD/L; GlobalPetrolPrices
#     via TheGlobalEconomy.com, fetched 2026-08). Pump fuel is the relevant
#     marginal energy price for diesel pumping and tracks national energy
#     subsidy regimes that also shape electricity tariffs.
#   - Water value proxies: (main) agricultural value added per hectare of
#     agricultural land (WB NV.AGR.TOTL.CD / AG.LND.AGRI.K2, mrnev);
#     (alt) GDP per capita (WB NY.GDP.PCAP.CD, mrnev).
#   - For each treated TBA-country segment with GW > 0, partner countries are
#     the other riparians of the same aquifer (from Data S2 metadata). For
#     each dyad: e = log(E_A/E_B), w = log(W_A/W_B). The PRIMARY dyad of a
#     segment is the partner with the largest |e| (sensitivity: mean over
#     partners). Dyad-level (orientation-free) quantities:
#       energy_mag = |e|;  aligned = (e * w > 0)
#   - 2x2 cells: energy_mag {Low, High; median split among treated segments}
#     x {Aligned, Misaligned}. Theory predicts the transboundary contrast is
#     strongest in the High x Aligned cell.
#   - Within each cell, the H1/H2 GW-specific contrasts are re-estimated with
#     the SAME machinery as the main pipeline (fit_one: full matching + lmer,
#     country RE, river-border covariate), over ENS control realizations.
#       H1: GW irrigation intensity (outcome GW), matched on cov_int
#       H2: borderward GW Gini   (outcome giniGw), matched on cov_gini
#     Controls: full pool each realization (controls carry no dyad).
#
# Outputs (the pipeline's canonical directories):
#   output/country_econ.csv                 country energy-price / water-value
#   output/asymmetry_dyads_<proxy>.csv      segment-level dyad metrics + cells
#   output/asymmetry_cells_<proxy>.csv      per-cell ensemble results
#   output/asymmetry_quadrant_tests_<proxy>.csv
#                                           FORMAL tests of differences ACROSS
#                                           quadrants (section 5)
#   figure/figure_asymmetry_<proxy>.pdf     2-panel cell-matrix figure
#   derived/cache/asym_*.rds                per-cell ensemble caches
#   derived/cache/quadtest_*.rds            per-design cross-quadrant test caches
# Inputs beyond the main data: 1_data/{gasoline_prices,wb_gdp_pc,wb_agva,
#   wb_agland}.csv
# =============================================================================

.args <- commandArgs(FALSE)
.f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.root <- if (length(.f)) dirname(normalizePath(.f[1])) else getwd()
if (!exists("CFG")) source(file.path(.root, "41_config.R"))
if (!exists("run_spec")) source(file.path(CFG$root, "42_core.R"))
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

# This analysis used to live in asymmetry_analysis/ with its own 1_data/,
# cache/, figure/, output/ and run log -- a second copy of the pipeline's
# directory layout. It now writes to the canonical directories like every
# other stage, so there is one cache, one output folder and one run log.
if (is.null(.log$path)) log_open(CFG$log_file)
for (d in c(CFG$out_dir, CFG$fig_dir, CFG$cache_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
cfg <- CFG
if (!exists("FORCE")) FORCE <- nzchar(Sys.getenv("FORCE"))
if (isTRUE(FORCE)) say("*** FORCE: rebuilding cached asymmetry cells. ***\n")

# cairo_pdf writes NOTHING (with only a warning) on a build without a working
# cairo -- on macOS, without XQuartz -- while the calling code reports success.
# Probe it by writing a real file, fall back to base pdf() with an encoding
# that can represent en dashes, and verify every output. Defined here rather
# than reused from 46_figures.R so this script still runs standalone.
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
if (.cairo_ok) {
  say("  PDF device: cairo_pdf (fonts embedded).\n")
} else {
  say("  PDF device: base pdf() -- cairo unavailable; fonts NOT embedded.\n")
}
save_pdf_asym <- function(p, name, w, h) {
  path <- file.path(cfg$fig_dir, paste0(name, ".pdf"))
  if (file.exists(path)) unlink(path)
  if (.cairo_ok) ggsave(path, p, width = w, height = h,
                        device = grDevices::cairo_pdf)
  else           ggsave(path, p, width = w, height = h,
                        device = grDevices::pdf, encoding = "WinAnsi.enc")
  if (!file.exists(path) || file.size(path) < 1000)
    stop("save_pdf_asym: ", basename(path), " was not written (device failure).")
  sayf("  wrote %s.pdf\n", name)
}
ENS <- 100L          # control realizations per cell (signal-hunting scale)
MIN_TREAT <- 8L      # minimum treated segments to estimate a cell

# =============================================================================
# 1. Country-level economic data + crosswalk
# =============================================================================
gas  <- read.csv(file.path(cfg$root, "1_data", "gasoline_prices.csv"))
gdp  <- read.csv(file.path(cfg$root, "1_data", "wb_gdp_pc.csv"))
agva <- read.csv(file.path(cfg$root, "1_data", "wb_agva.csv"))
agld <- read.csv(file.path(cfg$root, "1_data", "wb_agland.csv"))

# gasoline table country name -> iso3 (only names that differ from WB English)
gas_iso <- c(
  "UK"="GBR","USA"="USA","Domin. Rep."="DOM","C.A. Republic"="CAF",
  "R. of Congo"="COG","DR Congo"="COD","Tr.&Tobago"="TTO","UA Emirates"="ARE",
  "Burma"="MMR","Bosnia & Herz."="BIH","Ivory Coast"="CIV","Swaziland"="SWZ",
  "Cape Verde"="CPV","Czechia"="CZE","South Korea"="KOR","North Korea"="PRK",
  "Palestine"="PSE","North Macedonia"="MKD","Russia"="RUS","Iran"="IRN",
  "Syria"="SYR","Laos"="LAO","Vietnam"="VNM","Venezuela"="VEN","Egypt"="EGY",
  "Turkey"="TUR","Bolivia"="BOL","Moldova"="MDA","Tanzania"="TZA",
  "Brunei"="BRN","Hong Kong"="HKG","Slovakia"="SVK","Kyrgyzstan"="KGZ")
iso_of_gas <- function(nm) {
  hit <- gas_iso[nm]
  if (!is.na(hit)) return(unname(hit))
  # default: WB-style English name match against a small builtin table below
  unname(name_iso[nm])
}
# meta CountryName (IGRAC/GADM naming) -> iso3
name_iso <- c(
  "Afghanistan"="AFG","Albania"="ALB","Algeria"="DZA","Angola"="AGO",
  "Argentina"="ARG","Armenia"="ARM","Austria"="AUT","Azerbaijan"="AZE",
  "Bangladesh"="BGD","Belarus"="BLR","Belgium"="BEL","Belize"="BLZ",
  "Benin"="BEN","Bhutan"="BTN","Bolivia"="BOL","Bosnia and Herzegovina"="BIH",
  "Botswana"="BWA","Brazil"="BRA","Bulgaria"="BGR","Burkina Faso"="BFA",
  "Burundi"="BDI","Cambodia"="KHM","Cameroon"="CMR","Canada"="CAN",
  "Central African Republic"="CAF","Chad"="TCD","Chile"="CHL","China"="CHN",
  "Colombia"="COL","Costa Rica"="CRI","Croatia"="HRV","Cuba"="CUB",
  "Czech Republic"="CZE","Czechia"="CZE",
  "Democratic Republic of the Congo"="COD",
  "Denmark"="DNK","Djibouti"="DJI","Dominican Republic"="DOM",
  "Ecuador"="ECU","Egypt"="EGY","El Salvador"="SLV","Eritrea"="ERI",
  "Estonia"="EST","Ethiopia"="ETH","Finland"="FIN","France"="FRA",
  "Gabon"="GAB","Gambia"="GMB","Gaza Strip"="PSE","Georgia"="GEO",
  "Germany"="DEU","Ghana"="GHA","Greece"="GRC","Guatemala"="GTM",
  "Guinea"="GIN","Guinea-Bissau"="GNB","Guyana"="GUY","Haiti"="HTI",
  "Honduras"="HND","Hungary"="HUN","India"="IND","Indonesia"="IDN",
  "Iran"="IRN","Iraq"="IRQ","Ireland"="IRL","Israel"="ISR","Italy"="ITA",
  "Ivory Coast"="CIV","Cote d'Ivoire"="CIV","Jordan"="JOR",
  "Kazakhstan"="KAZ","Kenya"="KEN","Kuwait"="KWT","Kyrgyzstan"="KGZ",
  "Laos"="LAO","Latvia"="LVA","Lebanon"="LBN","Lesotho"="LSO",
  "Liberia"="LBR","Libya"="LBY","Lithuania"="LTU","Luxembourg"="LUX",
  "Macedonia"="MKD","North Macedonia"="MKD","Malawi"="MWI","Malaysia"="MYS",
  "Mali"="MLI","Mauritania"="MRT","Mexico"="MEX","Moldova"="MDA",
  "Mongolia"="MNG","Montenegro"="MNE","Morocco"="MAR","Mozambique"="MOZ",
  "Myanmar"="MMR","Burma"="MMR","Namibia"="NAM","Nepal"="NPL",
  "Netherlands"="NLD","Nicaragua"="NIC","Niger"="NER","Nigeria"="NGA",
  "North Korea"="PRK","Norway"="NOR","Oman"="OMN","Pakistan"="PAK",
  "Panama"="PAN","Papua New Guinea"="PNG","Paraguay"="PRY","Peru"="PER",
  "Poland"="POL","Portugal"="PRT","Qatar"="QAT","Republic of the Congo"="COG",
  "Republic of Congo"="COG","Congo"="COG","Romania"="ROU","Russia"="RUS",
  "Rwanda"="RWA","Saudi Arabia"="SAU","Senegal"="SEN","Serbia"="SRB",
  "Sierra Leone"="SLE","Slovakia"="SVK","Slovenia"="SVN","Somalia"="SOM",
  "South Africa"="ZAF","South Korea"="KOR","South Sudan"="SSD","Spain"="ESP",
  "Sudan"="SDN","Suriname"="SUR","Swaziland"="SWZ","Eswatini"="SWZ",
  "Sweden"="SWE","Switzerland"="CHE","Syria"="SYR","Tajikistan"="TJK",
  "Tanzania"="TZA","Thailand"="THA","Togo"="TGO","Tunisia"="TUN",
  "Turkey"="TUR","Turkmenistan"="TKM","Uganda"="UGA","Ukraine"="UKR",
  "United Arab Emirates"="ARE","United Kingdom"="GBR","United States"="USA",
  "Uruguay"="URY","Uzbekistan"="UZB","Venezuela"="VEN","Vietnam"="VNM",
  "West Bank"="PSE","Yemen"="YEM","Zambia"="ZMB","Zimbabwe"="ZWE",
  "Ecuador "="ECU","Timor-Leste"="TLS","Eritrea "="ERI")

gas$iso3 <- vapply(gas$country, iso_of_gas, character(1))
gas <- gas[!is.na(gas$iso3), ]
# Venezuela's administered price (~0) makes log-ratios degenerate; drop it and
# let segments fall back to their next partner.
gas$gasoline_usd_l[gas$iso3 == "VEN"] <- NA_real_

econ <- Reduce(function(a, b) merge(a, b, by = "iso3", all = TRUE), list(
  data.frame(iso3 = gas$iso3, energy = gas$gasoline_usd_l),
  data.frame(iso3 = gdp$iso3, gdp_pc = gdp$value),
  data.frame(iso3 = agva$iso3, agva = agva$value),
  data.frame(iso3 = agld$iso3, agland_km2 = agld$value)))
econ$agva_ha <- with(econ, agva / (agland_km2 * 100))   # USD per hectare
write.csv(econ, file.path(cfg$out_dir, "country_econ.csv"), row.names = FALSE)

# =============================================================================
# 2. Segment-level dyad asymmetries
# =============================================================================
meta <- load_meta(cfg)
meta$iso3 <- unname(name_iso[meta$CountryName])
unmapped <- sort(unique(meta$CountryName[is.na(meta$iso3)]))
if (length(unmapped)) say("  [crosswalk] unmapped countries (segments in these
  countries keep only mapped partners): ", paste(unmapped, collapse = "; "), "\n")

d  <- load_main_data(cfg)
tr <- d %>% filter(type == "treat", !is.na(GW), GW > 0)

# riparian sets per aquifer
rip <- meta %>% distinct(aq_id, CountryName, iso3)

dyad_metrics <- function(proxy = c("agva_ha", "gdp_pc")) {
  proxy <- match.arg(proxy)
  ec <- econ[, c("iso3", "energy", proxy)]; names(ec)[3] <- "wval"
  rows <- list()
  for (i in seq_len(nrow(tr))) {
    a_id <- tr$aq_id[i]; cA <- tr$CntrName[i]
    isoA <- name_iso[cA]; if (is.na(isoA)) next
    eA <- ec$energy[ec$iso3 == isoA]; wA <- ec$wval[ec$iso3 == isoA]
    if (!length(eA) || is.na(eA) || !length(wA) || is.na(wA)) next
    partners <- rip %>% filter(aq_id == a_id, CountryName != cA, !is.na(iso3))
    if (!nrow(partners)) next
    pm <- lapply(unique(partners$iso3), function(isoB) {
      eB <- ec$energy[ec$iso3 == isoB]; wB <- ec$wval[ec$iso3 == isoB]
      if (!length(eB) || is.na(eB) || !length(wB) || is.na(wB)) return(NULL)
      data.frame(isoB = isoB, e = log(eA / eB), w = log(wA / wB))
    })
    pm <- do.call(rbind, pm); if (is.null(pm) || !nrow(pm)) next
    j <- which.max(abs(pm$e))                    # primary dyad: largest |e|
    rows[[length(rows) + 1]] <- data.frame(
      aq_id = a_id, CntrName = cA, isoA = unname(isoA),
      partner = pm$isoB[j], n_partners = nrow(pm),
      e = pm$e[j], w = pm$w[j],
      e_mean = mean(pm$e), w_mean = mean(pm$w))
  }
  do.call(rbind, rows)
}

classify <- function(dy) {
  dy$energy_mag <- abs(dy$e)
  dy$aligned    <- dy$e * dy$w > 0
  med <- median(dy$energy_mag, na.rm = TRUE)
  dy$E_bin <- ifelse(dy$energy_mag > med, "High energy asym.", "Low energy asym.")
  dy$A_bin <- ifelse(dy$aligned, "Aligned", "Misaligned")
  dy$cell  <- paste(dy$E_bin, dy$A_bin, sep = " x ")
  attr(dy, "med_split") <- med
  dy
}

# =============================================================================
# 3. Per-cell estimation (H1: GW intensity; H2: GW borderward Gini)
# =============================================================================
HYPS <- list(
  H1 = list(outcome = "GW",     cov_match = "cov_int",
            title = "H1: GW irrigation intensity | Rivers"),
  H2 = list(outcome = "giniGw", cov_match = "cov_gini",
            title = "H2: Borderward GW Gini | Rivers"))

run_cell <- function(dy, cell, hyp, proxy, cfgl) {
  keep <- dy[dy$cell == cell, c("aq_id", "CntrName")]
  spec <- list(outcome = HYPS[[hyp]]$outcome, filter = "GW > 0",
               cov_match = HYPS[[hyp]]$cov_match, cov_reg = "river_cov",
               group = "asym", title = paste(hyp, cell))
  dd <- load_main_data(cfgl) %>% filter(!!rlang::parse_expr(spec$filter))
  dd <- dd %>% filter(type != "treat" |
                        paste(aq_id, CntrName) %in% paste(keep$aq_id, keep$CntrName))
  vars_keep <- unique(c(spec$outcome, "type", "aq_id", "CntrName",
                        resolve_cov(spec$cov_match, cfgl),
                        resolve_cov(spec$cov_reg, cfgl)))
  dd <- dd %>% select(all_of(vars_keep[vars_keep %in% names(dd)])) %>%
    filter(complete.cases(across(everything())))
  trdat <- dd %>% filter(type == "treat"); cndat <- dd %>% filter(type != "treat")
  n_tr <- nrow(trdat)
  if (n_tr < MIN_TREAT)
    return(list(cell = cell, hyp = hyp, n_treat = n_tr, estimable = FALSE))

  key  <- paste0("asym_", proxy, "_", hyp, "_", gsub("[^A-Za-z]", "", cell))
  path <- file.path(cfg$cache_dir, paste0(key, ".rds"))
  stamp <- cache_stamp(cfgl, spec = c(key, spec,
    list(ids = sort(paste(keep$aq_id, keep$CntrName)), ens = ENS)))
  hit <- if (isTRUE(FORCE)) NULL else cache_read(path, stamp, cfgl)
  if (!is.null(hit)) { say("  [cache] ", key, "\n"); return(hit) }

  sets <- load_control_sets(cfgl)[seq_len(ENS)]
  oplan <- future::plan(future::multisession, workers = par_workers(cfgl))
  on.exit(future::plan(oplan), add = TRUE)
  results <- future_map(sets, function(ids)
    fit_one(trdat, cndat, ids, spec, cfgl),
    .options = furrr_options(seed = TRUE))
  ens <- summarise_ensemble(results, cfgl)
  res <- list(cell = cell, hyp = hyp, n_treat = n_tr, estimable = TRUE,
              summary_df = ens$summary_df, counts = ens$counts,
              mean_se = median(ens$summary_df$treat_se, na.rm = TRUE))
  cache_write(res, path, stamp)
  sayf("  %-42s n_treat=%3d  eff=%+.4g  sig(neg/pos)=%d/%d%%\n",
       paste(hyp, cell), n_tr, res$counts$mean_eff,
       round(100 * res$counts$share_neg_sig), round(100 * res$counts$share_pos_sig))
  res
}

run_proxy <- function(proxy) {
  say("\n==== proxy: ", proxy, " ====\n")
  dy <- classify(dyad_metrics(proxy))
  write.csv(dy, file.path(cfg$out_dir,
            paste0("asymmetry_dyads_", proxy, ".csv")), row.names = FALSE)
  say("  segments classified: ", nrow(dy), "  (median |log E-ratio| split = ",
      round(attr(dy, "med_split"), 3), ")\n")
  print(table(dy$cell))
  cells <- sort(unique(dy$cell))
  out <- list()
  for (hyp in names(HYPS)) for (cl in cells)
    out[[paste(hyp, cl)]] <- run_cell(dy, cl, hyp, proxy, cfg)
  rows <- lapply(out, function(r) {
    base <- data.frame(hyp = r$hyp, cell = r$cell, n_treat = r$n_treat,
                       estimable = r$estimable)
    if (isTRUE(r$estimable)) cbind(base, r$counts) else base
  })
  tab <- dplyr::bind_rows(rows)
  write.csv(tab, file.path(cfg$out_dir,
            paste0("asymmetry_cells_", proxy, ".csv")), row.names = FALSE)
  list(dy = dy, tab = tab, res = out)
}

# =============================================================================
# 4. Figure: 2x2 cell matrices for H1 and H2
# =============================================================================
theme_paper <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 8, hjust = 0),
        plot.subtitle = element_text(size = 6.5, colour = "grey30"),
        plot.tag = element_text(face = "bold", size = 10),
        legend.position = "bottom", legend.key.height = unit(0.28, "cm"),
        legend.title = element_text(size = 6.5), legend.text = element_text(size = 6),
        axis.text = element_text(size = 7))

cell_panel <- function(tab, hyp, pp = FALSE, tag = NULL) {
  tt <- tab[tab$hyp == hyp, ]
  tt$E <- ifelse(grepl("^High", tt$cell), "High", "Low")
  tt$A <- ifelse(grepl("Aligned$", tt$cell) & !grepl("Misaligned$", tt$cell),
                 "Aligned", "Misaligned")
  sc <- if (pp) 100 else 1
  tt$lab <- ifelse(!tt$estimable | is.na(tt$mean_eff),
    sprintf("n = %d\n(not estimated)", tt$n_treat),
    sprintf("%+.3g\n[%+.3g, %+.3g]\n%d%% sig (n=%d)",
            tt$mean_eff * sc, tt$q05_eff * sc, tt$q95_eff * sc,
            round(100 * pmax(tt$share_neg_sig, tt$share_pos_sig)), tt$n_treat))
  tt$fill <- ifelse(tt$estimable, tt$mean_eff * sc, NA_real_)
  lim <- max(abs(tt$fill), na.rm = TRUE)
  ggplot(tt, aes(A, E)) +
    geom_tile(aes(fill = fill), colour = "white", linewidth = 1.2) +
    geom_tile(data = tt[tt$E == "High" & tt$A == "Aligned", ],
              fill = NA, colour = "grey15", linewidth = 0.9, linetype = 1) +
    geom_text(aes(label = lab), size = 2.3, lineheight = 1.0, colour = "grey10") +
    scale_fill_gradient2(if (pp) "Mean TBA difference (pp)" else "Mean TBA difference",
      low = "#0072B2", mid = "#F7F7F7", high = "#D55E00",
      midpoint = 0, limits = c(-lim, lim), na.value = "grey92") +
    scale_y_discrete(limits = c("Low", "High"),
                     labels = c("Low", "High")) +
    scale_x_discrete(limits = c("Misaligned", "Aligned")) +
    labs(x = "Energy-price vs water-value asymmetry",
         y = "Energy-price asymmetry  |log E-ratio|", title = HYPS[[hyp]]$title,
         subtitle = "outlined cell: theory-predicted maximum", tag = tag) +
    coord_fixed(0.8) + theme_paper
}

# -----------------------------------------------------------------------------
# Cross-quadrant test read-out, drawn UNDER each cell matrix (section 5)
# -----------------------------------------------------------------------------
# The cell matrix alone invites the reading it cannot support: four cells, four
# shades, and a reader concludes that the dark one differs from the pale one.
# It does not say so -- each cell is a separate contrast against its own matched
# controls, and "significant here, not there" is not a test of a difference.
# This strip puts the actual cross-quadrant tests on the same page as the cells,
# including the cases where they are NOT estimable, so the absence of a test is
# visible rather than inferred.
#
# Values are the CR2 tier, Rubin-pooled over the ensemble, from
# test_quadrant_differences(). H1 is rescaled to percentage points to match the
# cell labels; H2 (a Gini) is not.
.fmt_p <- function(p) if (!is.finite(p)) "NA" else
  if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)

test_strip <- function(qt, hyp, pp = FALSE) {
  sc  <- if (pp) 100 else 1
  rows <- list()
  hdr  <- "not available"
  if (!is.null(qt) && nrow(qt)) {
    q <- qt[qt$tier == "CR2" & qt$hyp == hyp, ]
    if (nrow(q)) {
      hdr <- sprintf(
        "Cross-quadrant tests -- CR2, country-clustered; Rubin-pooled over %d realizations",
        max(q$m_used, na.rm = TRUE))
      pull <- function(design, quantity) {
        r <- q[q$design == design & q$quantity == quantity, ]
        if (!nrow(r)) return(NA)
        r[1, ]
      }
      cellstr <- function(r, is_omni) {
        if (!is.data.frame(r)) return("--")
        flag <- if (nzchar(r$note)) " *" else ""
        if (is_omni)
          sprintf("F(%.0f, %.1f) = %.2f, %s%s", r$df1, r$df2, r$stat,
                  .fmt_p(r$p), flag)
        else
          sprintf("%+.3g (SE %.2g), %s%s", r$est * sc, r$se * sc,
                  .fmt_p(r$p), flag)
      }
      # "all four equal" is wrong when a quadrant was dropped -- the omnibus is
      # then over the surviving quadrants only, and saying so is the point.
      nq_used <- length(strsplit(q$quadrants[1], "\\+")[[1]])
      lab <- c(omnibus_all_equal = sprintf("all %d quadrants equal", nq_used),
               align = "alignment", mag = "energy magnitude",
               inter = "alignment x magnitude")
      for (qn in names(lab)) {
        is_o <- qn == "omnibus_all_equal"
        rs <- pull("stacked", qn); rp <- pull("pooled", qn)
        if (!is.data.frame(rs) && !is.data.frame(rp)) {
          rows[[qn]] <- c(lab[[qn]], "not estimable", "not estimable")
        } else {
          rows[[qn]] <- c(lab[[qn]], cellstr(rs, is_o), cellstr(rp, is_o))
        }
      }
      drop <- setdiff(QUADS, strsplit(q$quadrants[1], "\\+")[[1]])
      if (length(drop))
        hdr <- paste0(hdr, sprintf("\nquadrant%s %s excluded (n < %d): the 2x2 contrasts are not estimable",
                                   if (length(drop) > 1) "s" else "",
                                   paste(QUAD_LAB[drop], collapse = ", "), MIN_TREAT))
    }
  }
  if (!length(rows))
    rows <- list(c("all four equal", "not available", "not available"))
  n <- length(rows)
  dd <- do.call(rbind, lapply(seq_len(n), function(i)
    data.frame(row = i, col = 1:3, txt = rows[[i]])))
  dd$x <- c(0, 0.30, 0.655)[dd$col]
  dd$y <- -dd$row
  hd <- data.frame(x = c(0.30, 0.655), y = 0,
                   txt = c("stacked (primary)", "pooled (robustness)"))
  ggplot() +
    geom_text(data = hd, aes(x, y, label = txt), hjust = 0, size = 2.15,
              fontface = "bold", colour = "grey20") +
    geom_text(data = dd[dd$col == 1, ], aes(x, y, label = txt), hjust = 0,
              size = 2.15, colour = "grey20") +
    geom_text(data = dd[dd$col != 1, ], aes(x, y, label = txt), hjust = 0,
              size = 2.15, colour = "grey10") +
    scale_x_continuous(limits = c(0, 1.02)) +
    scale_y_continuous(limits = c(-n - 0.6, 1.3)) +
    labs(title = hdr) +
    theme_void(base_size = 8) +
    theme(plot.title = element_text(size = 6.2, colour = "grey30", hjust = 0,
                                    lineheight = 1.15,
                                    margin = margin(b = 2)),
          plot.margin = margin(2, 2, 2, 2))
}

make_fig <- function(rr, proxy, wlab, qt = NULL) {
  p1 <- cell_panel(rr$tab, "H1", pp = TRUE,  tag = "A")
  p2 <- cell_panel(rr$tab, "H2", pp = FALSE, tag = "B")
  t1 <- test_strip(qt, "H1", pp = TRUE)
  t2 <- test_strip(qt, "H2", pp = FALSE)
  fig <- ((p1 | p2) / (t1 | t2)) +
    plot_layout(heights = c(1, 0.42)) +
    plot_annotation(
    caption = paste0(
      "Dyad asymmetries: e = log energy-price ratio (pump gasoline, USD/L), ",
      "w = log ", wlab, " ratio; primary dyad = largest |e| partner.\n",
      "Cells: median split of |e| x sign-alignment of (e, w). Each cell: ",
      "transboundary-control contrast re-estimated over ", ENS,
      " control realizations\n(full matching + weighted mixed model, country ",
      "random intercept, river-border covariate), GW > 0 sample; ",
      "[5th, 95th] percentile of ensemble estimates; sig at p < 0.1.\n",
      "Lower strip: tests of DIFFERENCES between cells (section 5). ",
      "The three contrasts are mutually orthogonal and jointly equivalent to ",
      "the omnibus.\n\"*\" marks a denominator df below 5, where the ",
      "p-value carries little information. Cell shading is not a test of a ",
      "difference between cells."),
    theme = theme(plot.caption = element_text(size = 6, hjust = 0,
                                              colour = "grey30")))
  save_pdf_asym(fig, paste0("figure_asymmetry_", proxy), 8.2, 6.3)
}

# =============================================================================
# 5. TESTING DIFFERENCES *ACROSS* ASYMMETRY QUADRANTS
# =============================================================================
# WHY THIS SECTION EXISTS. Sections 1-4 estimate a transboundary contrast
# SEPARATELY in each quadrant and report, per quadrant, the ensemble mean and
# the share of realizations with p < sig_level. That answers "is the effect
# non-zero in THIS quadrant?". It does not answer "are the effects DIFFERENT
# ACROSS quadrants?", which is the claim the Case-3 theory actually makes.
# Four separate significance statements cannot be combined into a statement
# about their differences: the per-cell ensembles carry no joint covariance,
# and "significant here, not significant there" is not a test of a difference
# (Gelman & Stern 2006). Nor is the ensemble spread of a per-cell coefficient a
# sampling distribution -- it is control-selection variability with the treated
# sample held fixed. This section therefore estimates the quadrant effects
# JOINTLY, so that every cross-quadrant comparison is a linear contrast of
# coefficients from ONE fit with ONE covariance matrix.
#
# TWO DESIGNS, deliberately different in what they assume.
#
#   "stacked" (primary). Keep the per-quadrant matching exactly as sections 1-4
#   do it -- one matchit() per quadrant per realization, so each quadrant gets
#   its OWN matched control group and its own balance. Then STACK the four
#   matched datasets and fit a single model in which EVERY regressor is
#   quadrant-specific:
#
#       y ~ 0 + b_Q(4 quadrant control baselines)
#             + t_Q(4 quadrant transboundary effects)
#             + river_cov_Q(4 quadrant nuisance slopes) + (1|CntrName)
#
#   The design matrix is then block diagonal, and the weighted OLS fit
#   REPRODUCES THE FOUR PER-QUADRANT REGRESSIONS EXACTLY -- verified at machine
#   precision on every run by .verify_block_equivalence(). So the joint fit
#   changes no point estimate; it adds exactly one thing, a covariance BETWEEN
#   the quadrants, which is what a cross-quadrant contrast needs and what four
#   separate fits cannot supply. (Constraining the nuisance slope to be common
#   across quadrants was tried and rejected: it moved the quadrant effects by
#   more than the effects themselves.)
#
#   The coefficient t_Q is quadrant Q's transboundary-vs-matched-control
#   contrast; t_Q - t_Q' is therefore a DIFFERENCE IN DIFFERENCES. This is the
#   estimand we want: quadrants differ systematically in geography and aridity,
#   and each quadrant's own matched controls absorb those differences before the
#   comparison is taken. Cost: a control basin can serve several quadrants, so
#   rows are reused across strata. That reuse is WITHIN country (a basin has one
#   country), so CR2 standard errors clustered on country subsume it -- which is
#   also the clustering the rest of this pipeline uses (42_core.R::refit_best).
#   The country random intercept of the "mixed" tier IS shared across quadrants,
#   so its t_Q are partially pooled and differ from the section 1-3 per-cell
#   values; the CR2 tier does not have that property and is the one to quote.
#
#   "pooled" (robustness). One matchit() on ALL treated GW > 0 segments per
#   realization -- i.e. the matched dataset of the existing GWRivsInt (H1) and
#   GWRivsGini (H2) specifications -- then
#
#       y ~ 1 + t_Q(4) + river_cov + (1|CntrName)
#
#   with a single shared control baseline (the intercept). The pooled treatment
#   effect implied by this fit is the reported main-text GW contrast, so the
#   quadrant effects here decompose a number already in the paper. Because the
#   control baseline is common, t_Q - t_Q' reduces to a comparison AMONG treated
#   segments, adjusted only by the pooled matching weights and river_cov: fewer
#   moving parts, but no quadrant-specific covariate adjustment. Agreement
#   between the two designs is the thing worth reporting.
#
# WHAT IS TESTED. With quadrants ordered (LM, LA, HM, HA) = (Low x Misaligned,
# Low x Aligned, High x Misaligned, High x Aligned), the 3-df null "all four
# quadrant effects are equal" is decomposed into three MUTUALLY ORTHOGONAL 1-df
# contrasts that exactly span it, each with a theory reading:
#
#   align  = (LA + HA)/2 - (LM + HM)/2   alignment main effect
#   mag    = (HM + HA)/2 - (LM + LA)/2   energy-magnitude main effect
#   inter  = (HA - HM) - (LA - LM)       alignment x magnitude interaction
#
# `inter` is the sharp Case-3 prediction: amplification requires the cost
# asymmetry to dominate the value asymmetry (B > A), which only happens among
# ALIGNED dyads with a LARGE energy asymmetry -- so the theory predicts an
# interaction, not a main effect (README_asymmetry.md, points 2-3). Because the
# three are orthogonal and jointly equivalent to the omnibus, this is three
# tests, not six pairwise ones: no multiplicity correction beyond reporting the
# omnibus alongside them. The omnibus F is reported first and should be read as
# the gatekeeper.
#
# COMBINING THE ENSEMBLE. Each realization r gives a contrast estimate and its
# within-realization variance. These are pooled by RUBIN'S RULES -- the standard
# rule for combining an estimate across repeated matched/imputed samples:
#
#   T = Ubar + (1 + 1/m) B,  Ubar = mean within-variance, B = between-variance
#
# so the total standard error counts BOTH the uncertainty of a single matched
# comparison and the control-selection spread that the per-cell ensembles have
# been displaying all along. Degrees of freedom are Barnard-Rubin. The omnibus
# uses the multivariate analogue (the D1 statistic of Li, Raghunathan & Rubin
# 1991), which is invariant to the choice of basis for the contrast space. The
# realizations are alternative non-overlapping control assignments rather than
# independent posterior draws, so the pooled p is an approximation and should be
# described as such; it is the conservative direction, since B is added to Ubar.
#
# The ensemble mean / 5th / 95th percentile of the per-realization contrast and
# the share of realizations with p < sig_level are ALSO reported, in the same
# idiom as the rest of the figure -- but they are DESCRIPTIVE columns, not the
# test. The test is the Rubin-pooled estimate, SE, df and p.
#
# TWO INFERENCE TIERS per design, exactly as elsewhere in this pipeline:
#   "mixed"  weighted lmer with a country random intercept, Satterthwaite df
#            for the contrast (lmerTest::contest);
#   "CR2"    the SAME weighted fixed-effect model fitted by OLS with CR2
#            country-clustered SEs (clubSandwich), which is what handles the
#            reuse of control rows across quadrant strata in the stacked design.
# CR2 is the one to quote; the mixed tier is reported so the numbers can be
# lined up against the per-cell estimates in sections 1-3.
#
# LIMITS, to be stated wherever these numbers are used. Per-quadrant treated
# samples are 10-49 segments; quadrant membership is built from CONTEMPORARY
# national aggregates and from a median split; the ensemble is a
# control-selection device, not a resample of the treated population. A
# non-significant omnibus here means "these data cannot resolve the difference",
# not "the quadrants are the same".
#
# Outputs:
#   output/asymmetry_quadrant_tests_<proxy>.csv
#   derived/cache/quadtest_<proxy>_<hyp>.rds
# =============================================================================

QUADS    <- c("LM", "LA", "HM", "HA")
QUAD_LAB <- c(LM = "Low x Misaligned",  LA = "Low x Aligned",
              HM = "High x Misaligned", HA = "High x Aligned")

# cell strings are paste(E_bin, A_bin, sep = " x "); test "Misaligned" FIRST,
# because "Aligned$" also matches "Misaligned".
quad_code <- function(cell) {
  paste0(ifelse(grepl("^High", cell), "H", "L"),
         ifelse(grepl("Misaligned$", cell), "M", "A"))
}

# The orthogonal basis. Rows are contrasts, columns are quadrants.
QUAD_CONTRASTS <- rbind(
  align = c(LM = -0.5, LA =  0.5, HM = -0.5, HA =  0.5),
  mag   = c(LM = -0.5, LA = -0.5, HM =  0.5, HA =  0.5),
  inter = c(LM =  1.0, LA = -1.0, HM = -1.0, HA =  1.0))
stopifnot(all(abs(QUAD_CONTRASTS %*% t(QUAD_CONTRASTS) -
                  diag(diag(QUAD_CONTRASTS %*% t(QUAD_CONTRASTS)))) < 1e-12))

# -----------------------------------------------------------------------------
# 5a. Rubin's rules
# -----------------------------------------------------------------------------
# Scalar pooling with Barnard-Rubin degrees of freedom. df_com is the
# complete-data (single-realization) df; when it is unavailable the classical
# Rubin df is used, which is the larger and therefore less conservative of the
# two -- so df_com is passed in wherever the fitting tier can supply it.
mi_pool_scalar <- function(est, var, df_com = NULL) {
  ok  <- is.finite(est) & is.finite(var) & var > 0
  est <- est[ok]; var <- var[ok]
  m   <- length(est)
  na  <- list(est = NA_real_, se = NA_real_, df = NA_real_, stat = NA_real_,
              p = NA_real_, lambda = NA_real_, m = m,
              se_within = NA_real_, sd_between = NA_real_)
  if (m < 2L) return(na)
  qbar <- mean(est); ubar <- mean(var); b <- stats::var(est)
  tvar <- ubar + (1 + 1/m) * b
  if (!is.finite(tvar) || tvar <= 0) return(na)
  lambda <- min(max((1 + 1/m) * b / tvar, 1e-10), 1 - 1e-10)
  df_old <- (m - 1) / lambda^2
  dfc    <- if (is.null(df_com)) NA_real_ else stats::median(df_com, na.rm = TRUE)
  df <- if (is.finite(dfc) && dfc > 0) {
    df_obs <- ((dfc + 1) / (dfc + 3)) * dfc * (1 - lambda)
    df_old * df_obs / (df_old + df_obs)
  } else df_old
  t <- qbar / sqrt(tvar)
  # `lambda` is the share of the TOTAL variance contributed by control
  # selection -- mice calls this lambda, not fmi (its `fmi` adds the df term).
  list(est = qbar, se = sqrt(tvar), df = df, stat = t,
       p = 2 * stats::pt(-abs(t), df), lambda = lambda, m = m,
       se_within = sqrt(ubar), sd_between = sqrt(b))
}

# Multivariate pooling: the D1 statistic (Li, Raghunathan & Rubin 1991), the
# multiple-imputation analogue of a Wald test. est_list holds the k-vector of
# contrast estimates per realization, U_list the k x k within-realization
# covariance. D1 is invariant to a non-singular change of basis of the contrast
# space, so any spanning set of the "all quadrants equal" null gives the same F.
mi_pool_wald <- function(est_list, U_list, df_com = NULL) {
  keep <- vapply(seq_along(est_list), function(i)
    all(is.finite(est_list[[i]])) && all(is.finite(U_list[[i]])), logical(1))
  est_list <- est_list[keep]; U_list <- U_list[keep]
  m <- length(est_list)
  na <- list(stat = NA_real_, df1 = NA_real_, df2 = NA_real_, p = NA_real_,
             rbar = NA_real_, m = m)
  if (m < 2L) return(na)
  k    <- length(est_list[[1]])
  # With a single constraint the omnibus IS that contrast, so route it through
  # the scalar pooler: same F (D1 reduces to t^2 exactly -- verified) but the
  # Barnard-Rubin df rather than the Li et al. one, which keeps the omnibus row
  # and the single contrast row from carrying two different p-values for one
  # hypothesis.
  if (k == 1L) {
    ps <- mi_pool_scalar(vapply(est_list, function(z) z[1], numeric(1)),
                         vapply(U_list,   function(z) z[1, 1], numeric(1)),
                         df_com)
    return(list(stat = ps$stat^2, df1 = 1, df2 = ps$df, p = ps$p,
                rbar = NA_real_, m = ps$m))
  }
  Q    <- Reduce(`+`, est_list) / m
  Ubar <- Reduce(`+`, U_list)   / m
  B    <- stats::cov(do.call(rbind, est_list))
  Uinv <- tryCatch(solve(Ubar), error = function(e) NULL)
  if (is.null(Uinv)) return(na)
  rbar <- (1 + 1/m) * sum(diag(B %*% Uinv)) / k
  if (!is.finite(rbar) || rbar <= 0) rbar <- 1e-8
  D1 <- as.numeric(t(Q) %*% Uinv %*% Q) / (k * (1 + rbar))
  tt <- k * (m - 1)
  nu <- if (tt > 4) 4 + (tt - 4) * (1 + (1 - 2/tt) / rbar)^2
        else        tt * (1 + 1/k) * (1 + 1/rbar)^2 / 2
  # The Li et al. denominator df diverges as the between-realization variance
  # goes to zero: with rbar ~ 0 it returns 1e8-ish, which is arithmetically
  # right (no control-selection uncertainty left, so the pooled test degenerates
  # to the single-realization test) but reads as nonsense in a table. In that
  # limit the correct reference distribution is the COMPLETE-DATA one, so cap
  # the pooled df at the (median) complete-data denominator df. The cap can only
  # shrink df, i.e. it is never anti-conservative. Same role as the Reiter (2007)
  # small-sample adjustment to D1.
  dfc <- if (is.null(df_com)) NA_real_ else stats::median(df_com, na.rm = TRUE)
  if (is.finite(dfc) && dfc > 0) nu <- min(nu, dfc)
  list(stat = D1, df1 = k, df2 = nu,
       p = stats::pf(D1, k, nu, lower.tail = FALSE), rbar = rbar, m = m)
}

# Regression guard, pure arithmetic and instant. The reference numbers are
# mice::pool.scalar(Q = 1:5/10, U = rep(0.04, 5), n = 51, k = 1), i.e. an
# independent implementation of Rubin's rules with Barnard-Rubin df; and the
# fact that D1 with a single constraint reduces exactly to the squared pooled t.
local({
  p <- mi_pool_scalar(c(0.1, 0.2, 0.3, 0.4, 0.5), rep(0.04, 5), rep(50, 5))
  stopifnot(abs(p$est - 0.3) < 1e-12,
            abs(p$se  - 0.264575131106459) < 1e-10,
            abs(p$df  - 12.1520095309879)  < 1e-8)
  w <- mi_pool_wald(as.list(c(0.1, 0.2, 0.3, 0.4, 0.5)),
                    rep(list(matrix(0.04, 1, 1)), 5), rep(50, 5))
  stopifnot(abs(w$stat - p$stat^2) < 1e-10, abs(w$p - p$p) < 1e-12)
})

# -----------------------------------------------------------------------------
# 5b. One realization: build the joint design, fit both tiers, return contrasts
# -----------------------------------------------------------------------------
# Explicit dummy columns rather than a `quad * type` formula: the coefficient
# NAMES are then known in advance ("t_HA", ...), so the contrast matrix can be
# placed by name and a rank-deficiency drop by lme4/lm shows up as a missing
# name instead of silently shifting which column a contrast weight lands on.
.quad_design <- function(md, quads, spec, cfgl, stacked) {
  md$.tb <- as.integer(md$type == "treat")
  for (q in quads) {
    if (stacked) md[[paste0("b_", q)]] <- as.numeric(md$quad == q)
    md[[paste0("t_", q)]] <- as.numeric(md$quad == q & md$.tb == 1L)
  }
  cov_reg <- resolve_cov(spec$cov_reg, cfgl)
  cov_reg <- cov_reg[cov_reg %in% names(md)]

  if (!stacked) {
    # Pooled: one shared control baseline (the intercept) and one shared
    # nuisance adjustment. Controls carry no quadrant, so a quadrant-specific
    # river coefficient is not identified here -- that parsimony is the point of
    # this design, and is why it is the robustness column rather than the
    # primary one.
    return(list(md = md,
                form = reformulate(c(paste0("t_", quads), cov_reg),
                                   response = spec$outcome),
                tnames = paste0("t_", quads)))
  }

  # Stacked: EVERY regressor is made quadrant-specific -- baseline, treatment
  # AND the nuisance covariates. Two reasons, one statistical and one practical.
  #
  #   Statistical. With quadrant-specific nuisance columns the design matrix is
  #   block diagonal (each row belongs to exactly one quadrant), so the weighted
  #   OLS fit reproduces the four per-quadrant regressions EXACTLY, coefficient
  #   for coefficient. The joint fit therefore adds nothing to the point
  #   estimates and exactly one thing to the inference: a covariance BETWEEN the
  #   quadrants, which is what a cross-quadrant contrast needs and what four
  #   separate fits cannot supply. Sharing the river coefficient across
  #   quadrants instead was tried and rejected: it moved the quadrant effects by
  #   up to 1.4 pp (H1), i.e. more than the effects themselves, so a "difference
  #   across quadrants" would partly have been an artefact of the constraint.
  #
  #   Practical. `0 + b_LM + ... + river_cov` looks equivalent and is a trap:
  #   with no intercept term model.matrix() gives the first FACTOR or LOGICAL
  #   predictor -- river_cov is logical -- FULL dummy coding, and those two
  #   columns sum to 1, which is collinear with the four b_Q that also sum to 1.
  #   lme4 then quietly drops a column. Building every column as an explicit
  #   numeric interaction removes every factor from the formula, so the
  #   no-intercept parameterisation is safe and the coefficient NAMES are known
  #   in advance.
  #
  # A quadrant-specific nuisance column is dropped when it is degenerate inside
  # its own block (no rows, or all rows, near a river border) -- exactly the
  # condition under which the per-quadrant regression would have dropped it.
  rcols <- character(0)
  for (q in quads) for (v in cov_reg) {
    nm <- paste0("x_", v, "_", q)
    md[[nm]] <- as.numeric(md$quad == q) * as.numeric(md[[v]])
    nz <- sum(md[[nm]] != 0)
    if (nz == 0 || nz == sum(md$quad == q)) md[[nm]] <- NULL else rcols <- c(rcols, nm)
  }
  list(md = md,
       form = reformulate(c("0", paste0("b_", quads), paste0("t_", quads), rcols),
                          response = spec$outcome),
       tnames = paste0("t_", quads))
}

# Small-sample (HTZ / Tipton-Pustejovsky) denominator df for each contrast row,
# in ONE Wald_test call with the already-computed CR2 matrix -- the adjustment
# matrices are the expensive part and are then built once per realization.
# Two version hazards, both handled rather than assumed away: the denominator-df
# column has been named both `df` and `df_denom`, and tidy = TRUE dispatches
# through mapply on names(constraints), so the list MUST be named. 42_core.R
# already carries the same kind of version-compat code for coef_test.
.wald_dfs <- function(ols, L, V) {
  cl <- lapply(seq_len(nrow(L)), function(j) L[j, , drop = FALSE])
  names(cl) <- rownames(L)
  w <- tryCatch(as.data.frame(
    clubSandwich::Wald_test(ols, constraints = cl, vcov = V,
                            test = "HTZ", tidy = TRUE)),
    error = function(e) NULL)
  nm <- if (is.null(w)) character(0) else intersect(c("df_denom", "df"), names(w))
  if (!length(nm) || nrow(w) != nrow(L)) return(rep(NA_real_, nrow(L)))
  as.numeric(w[[nm[1]]])
}

# Returns, per inference tier, the contrast estimates, their within-realization
# variances and df, and the full k x k covariance of the omnibus rows.
.fit_joint_one <- function(md, quads, spec, cfgl, stacked, Call, omni_rows) {
  bad <- function(why) list(ok = FALSE, why = why)
  dz  <- .quad_design(md, quads, spec, cfgl, stacked)

  ## --- tier 1: weighted mixed model, Satterthwaite df -----------------------
  # Warnings are suppressed, not caught: a singular country variance or a
  # convergence nudge is routine on a 200-row stacked design and is NOT a reason
  # to discard the realization. Errors are.
  fit <- tryCatch(suppressWarnings(
    lmerTest::lmer(update(dz$form, . ~ . + (1 | CntrName)),
                   data = dz$md, weights = weights_combined)),
    error = function(e) NULL)
  if (is.null(fit)) return(bad("lmer failed"))
  fx <- lme4::fixef(fit)
  if (!all(dz$tnames %in% names(fx))) return(bad("rank-deficient lmer design"))
  Lm <- matrix(0, nrow(Call), length(fx),
               dimnames = list(rownames(Call), names(fx)))
  Lm[, match(dz$tnames, names(fx))] <- Call
  ct <- tryCatch(as.data.frame(
    lmerTest::contest(fit, L = Lm, joint = FALSE, confint = FALSE,
                      check_estimability = TRUE)),
    error = function(e) NULL)
  if (is.null(ct) || nrow(ct) != nrow(Call)) return(bad("contest failed"))
  Vm <- as.matrix(stats::vcov(fit))
  Um <- Lm %*% Vm %*% t(Lm)
  jm <- tryCatch(as.data.frame(
    lmerTest::contest(fit, L = Lm[omni_rows, , drop = FALSE], joint = TRUE)),
    error = function(e) NULL)
  mixed <- list(est = ct[["Estimate"]],
                var = ct[["Std. Error"]]^2,
                df  = ct[["df"]],
                U   = Um[omni_rows, omni_rows, drop = FALSE],
                df_omni = if (is.null(jm)) NA_real_ else
                  suppressWarnings(as.numeric(jm[["DenDF"]][1])))

  ## --- tier 2: same fixed-effect model by weighted OLS + CR2 on country -----
  ols <- tryCatch(stats::lm(dz$form, data = dz$md, weights = weights_combined),
                  error = function(e) NULL)
  cr2 <- NULL
  if (!is.null(ols)) {
    bo <- stats::coef(ols)
    if (all(dz$tnames %in% names(bo)) && !any(is.na(bo))) {
      V <- tryCatch(clubSandwich::vcovCR(ols, cluster = dz$md$CntrName,
                                         type = "CR2"),
                    error = function(e) NULL)
      if (!is.null(V) && nrow(V) == length(bo)) {
        Lo <- matrix(0, nrow(Call), length(bo),
                     dimnames = list(rownames(Call), names(bo)))
        Lo[, match(dz$tnames, names(bo))] <- Call
        Uo  <- Lo %*% V %*% t(Lo)
        dfo <- .wald_dfs(ols, Lo, V)
        if (all(!is.finite(dfo)))
          dfo <- rep(dplyr::n_distinct(dz$md$CntrName) - 1, nrow(Lo))
        wo <- tryCatch(as.data.frame(clubSandwich::Wald_test(
          ols, constraints = Lo[omni_rows, , drop = FALSE], vcov = V,
          test = "HTZ")), error = function(e) NULL)
        nmo <- if (is.null(wo)) character(0) else
          intersect(c("df_denom", "df"), names(wo))
        cr2 <- list(est = as.vector(Lo %*% bo), var = diag(Uo), df = dfo,
                    U = Uo[omni_rows, omni_rows, drop = FALSE],
                    df_omni = if (length(nmo)) as.numeric(wo[[nmo[1]]][1])
                              else NA_real_)
      }
    }
  }
  if (is.null(cr2)) return(bad("CR2 tier failed"))

  list(ok = TRUE, mixed = mixed, CR2 = cr2,
       n_row = nrow(dz$md), n_country = dplyr::n_distinct(dz$md$CntrName))
}

# stacked: one matchit per quadrant, then bind. Each quadrant keeps its own
# matched control group, so `quad` labels the treated segment AND the controls
# matched to it.
.realization_stacked <- function(tr_by_quad, cndat, ids, spec, cfgl, quads,
                                 Call, omni_rows) {
  parts <- vector("list", length(quads)); names(parts) <- quads
  smd   <- stats::setNames(rep(NA_real_, length(quads)), quads)
  for (q in quads) {
    b <- fit_one(tr_by_quad[[q]], cndat, ids, spec, cfgl, keep_fit = TRUE)
    if (!isTRUE(b$success))
      return(list(ok = FALSE, why = paste0("quadrant ", q, ": ", b$error)))
    smd[q]     <- b$max_smd
    parts[[q]] <- b$mdat %>%
      dplyr::select(-dplyr::any_of(c("subclass", "distance", "weights"))) %>%
      dplyr::mutate(quad = q)
  }
  out <- .fit_joint_one(dplyr::bind_rows(parts), quads, spec, cfgl,
                        stacked = TRUE, Call = Call, omni_rows = omni_rows)
  out$max_smd <- smd
  out
}

# pooled: one matchit on ALL treated segments -- the matched dataset of the
# existing GWRivsInt / GWRivsGini specification -- with a single shared control
# baseline. Controls carry quad = "_ctrl" (a level, not NA, so nothing in
# MatchIt or model.frame has to drop rows).
.realization_pooled <- function(trdat_all, cndat, ids, spec, cfgl, quads,
                                Call, omni_rows) {
  b <- fit_one(trdat_all, cndat, ids, spec, cfgl, keep_fit = TRUE)
  if (!isTRUE(b$success)) return(list(ok = FALSE, why = b$error))
  md <- b$mdat %>%
    dplyr::select(-dplyr::any_of(c("subclass", "distance", "weights")))
  out <- .fit_joint_one(md, quads, spec, cfgl, stacked = FALSE,
                        Call = Call, omni_rows = omni_rows)
  out$max_smd <- stats::setNames(rep(b$max_smd, length(quads)), quads)
  out
}

# Assertion, run once per hypothesis rather than per realization. The stacked
# design claims to BE the four per-quadrant regressions with a covariance
# attached; this checks that claim on the first control realization instead of
# asserting it in a comment. If it ever fails, a nuisance covariate has been
# constrained across quadrants and the quadrant effects are no longer the ones
# sections 1-3 report.
.verify_block_equivalence <- function(tr_by_quad, cndat, ids, spec, cfgl, quads) {
  parts <- list(); sep <- stats::setNames(rep(NA_real_, length(quads)), quads)
  cov_reg <- resolve_cov(spec$cov_reg, cfgl)
  for (q in quads) {
    b <- fit_one(tr_by_quad[[q]], cndat, ids, spec, cfgl, keep_fit = TRUE)
    if (!isTRUE(b$success)) return(invisible(NA))
    sep[q] <- unname(stats::coef(stats::lm(
      reformulate(c("type", cov_reg), response = spec$outcome),
      data = b$mdat, weights = weights_combined))[2])
    parts[[q]] <- b$mdat %>%
      dplyr::select(-dplyr::any_of(c("subclass", "distance", "weights"))) %>%
      dplyr::mutate(quad = q)
  }
  dz  <- .quad_design(dplyr::bind_rows(parts), quads, spec, cfgl, stacked = TRUE)
  ols <- stats::lm(dz$form, data = dz$md, weights = weights_combined)
  joint <- stats::coef(ols)[dz$tnames]
  d <- max(abs(joint - sep))
  sayf("    block-diagonality check: max |joint - per-quadrant OLS| = %.2e\n", d)
  if (!(is.finite(d) && d < 1e-8))
    stop("49_asymmetry: the stacked design no longer reproduces the ",
         "per-quadrant regressions (max |diff| = ", format(d), ").")
  invisible(d)
}

# -----------------------------------------------------------------------------
# 5c. Driver: one hypothesis x one design x the whole ensemble
# -----------------------------------------------------------------------------
.pool_design <- function(res, Call, row_kind, quads, alpha) {
  ok  <- vapply(res, function(r) isTRUE(r$ok), logical(1))
  res <- res[ok]
  m   <- length(res)
  if (!m) return(NULL)
  omni_rows <- which(row_kind == "omnibus_basis")
  out <- list()
  for (tier in c("mixed", "CR2")) {
    E <- do.call(rbind, lapply(res, function(r) r[[tier]]$est))
    V <- do.call(rbind, lapply(res, function(r) r[[tier]]$var))
    D <- do.call(rbind, lapply(res, function(r) r[[tier]]$df))
    rows <- lapply(seq_len(nrow(Call)), function(j) {
      pl <- mi_pool_scalar(E[, j], V[, j], D[, j])
      pr <- 2 * stats::pt(-abs(E[, j]) / sqrt(V[, j]), D[, j])
      tibble::tibble(
        tier = tier, quantity = rownames(Call)[j], kind = row_kind[j],
        est = pl$est, se = pl$se, stat = pl$stat, stat_type = "t",
        df1 = NA_real_, df2 = pl$df, p = pl$p,
        lambda = pl$lambda, se_within = pl$se_within, sd_between = pl$sd_between,
        m_used = pl$m,
        ens_mean = mean(E[, j], na.rm = TRUE),
        ens_q05  = unname(stats::quantile(E[, j], 0.05, na.rm = TRUE)),
        ens_q95  = unname(stats::quantile(E[, j], 0.95, na.rm = TRUE)),
        ens_share_sig = mean(pr < alpha, na.rm = TRUE))
    })
    om <- mi_pool_wald(lapply(res, function(r) r[[tier]]$est[omni_rows]),
                       lapply(res, function(r) r[[tier]]$U),
                       vapply(res, function(r) {
                         z <- r[[tier]]$df_omni
                         if (is.null(z)) NA_real_ else as.numeric(z)[1] },
                         numeric(1)))
    rows[[length(rows) + 1]] <- tibble::tibble(
      tier = tier, quantity = "omnibus_all_equal", kind = "omnibus_test",
      est = NA_real_, se = NA_real_, stat = om$stat, stat_type = "F",
      df1 = om$df1, df2 = om$df2, p = om$p,
      lambda = NA_real_, se_within = NA_real_, sd_between = NA_real_,
      m_used = om$m, ens_mean = NA_real_, ens_q05 = NA_real_,
      ens_q95 = NA_real_, ens_share_sig = NA_real_)
    out[[tier]] <- dplyr::bind_rows(rows)
  }
  smd <- do.call(rbind, lapply(res, function(r) r$max_smd))
  list(tab = dplyr::bind_rows(out),
       m = m,
       max_smd_med = apply(smd, 2, stats::median, na.rm = TRUE),
       n_country = stats::median(vapply(res, function(r) r$n_country, numeric(1))))
}

# Main entry point. Returns a tidy table; also cached and written to output/.
#   dy     the classified dyad table from classify(dyad_metrics(proxy))
#   proxy  "agva_ha" | "gdp_pc"  (labels the outputs only)
#   ens    number of control realizations (default ENS, as in sections 1-3)
test_quadrant_differences <- function(dy, proxy, cfgl = cfg, ens = ENS,
                                      designs = c("stacked", "pooled")) {
  dy$quad <- quad_code(dy$cell)
  dy$key  <- paste(dy$aq_id, dy$CntrName)
  say("\n---- quadrant difference tests: ", proxy, " ----\n")

  sets <- load_control_sets(cfgl)
  stopifnot(length(sets) >= ens)
  sets <- sets[seq_len(ens)]
  alpha <- cfgl$sig_level
  all_tabs <- list()

  for (hyp in names(HYPS)) {
    spec <- list(outcome = HYPS[[hyp]]$outcome, filter = "GW > 0",
                 cov_match = HYPS[[hyp]]$cov_match, cov_reg = "river_cov",
                 group = "asym", title = paste(hyp, "quadrant tests"))

    ## ---- spec data, complete cases, quadrant labels ------------------------
    dd <- load_main_data(cfgl) %>% filter(!!rlang::parse_expr(spec$filter))
    vars_keep <- unique(c(spec$outcome, "type", "aq_id", "CntrName",
                          resolve_cov(spec$cov_match, cfgl),
                          resolve_cov(spec$cov_reg, cfgl)))
    dd <- dd %>% select(all_of(vars_keep[vars_keep %in% names(dd)])) %>%
      filter(complete.cases(across(everything())))
    trd <- dd %>% filter(type == "treat") %>%
      mutate(quad = unname(stats::setNames(dy$quad, dy$key)[paste(aq_id, CntrName)]))
    cnd <- dd %>% filter(type != "treat")
    n_unclassified <- sum(is.na(trd$quad))
    trd <- trd %>% filter(!is.na(quad))

    nq  <- table(factor(trd$quad, levels = QUADS))
    used <- QUADS[nq >= MIN_TREAT]
    sayf("  %s: %d treated classified (%d dropped: no economic data), n = %s\n",
         hyp, nrow(trd), n_unclassified,
         paste(sprintf("%s:%d", QUADS, as.integer(nq)), collapse = " "))
    if (length(used) < 2L) {
      say("    fewer than two estimable quadrants -- no test possible.\n")
      next
    }
    if (length(used) < length(QUADS))
      sayf("    quadrants below MIN_TREAT=%d excluded from the joint model: %s\n",
           MIN_TREAT, paste(setdiff(QUADS, used), collapse = ", "))

    ## ---- contrast rows: per-quadrant effects, named contrasts, omnibus -----
    k    <- length(used)
    Cid  <- diag(k); rownames(Cid) <- paste0("beta_", used)
    # a named contrast survives only if every quadrant it weights is estimable
    gone  <- setdiff(QUADS, used)
    keepC <- if (!length(gone)) rep(TRUE, nrow(QUAD_CONTRASTS)) else
      apply(QUAD_CONTRASTS[, gone, drop = FALSE], 1, function(z) all(z == 0))
    Cnm  <- QUAD_CONTRASTS[keepC, used, drop = FALSE]
    # A (k-1)-row basis for the "all quadrants equal" null. These rows are also
    # readable on their own as pairwise differences against the first quadrant;
    # they are the omnibus BASIS, not a multiplicity-controlled pairwise family,
    # and are labelled as such. D1 is invariant to the choice of basis.
    Comni <- do.call(rbind, lapply(2:k, function(j) {
      z <- rep(0, k); z[1] <- 1; z[j] <- -1; z }))
    rownames(Comni) <- paste0("diff_", used[1], "_minus_", used[-1])
    Call <- rbind(Cid, Cnm, Comni)
    colnames(Call) <- used
    row_kind <- c(rep("quadrant_effect", nrow(Cid)),
                  rep("named_contrast",  nrow(Cnm)),
                  rep("omnibus_basis",   nrow(Comni)))
    if (!all(rownames(QUAD_CONTRASTS) %in% rownames(Cnm)))
      sayf("    named contrasts not estimable (need a dropped quadrant): %s\n",
           paste(setdiff(rownames(QUAD_CONTRASTS), rownames(Cnm)),
                 collapse = ", "))

    tr_by_quad <- lapply(used, function(q) trd %>% filter(quad == q))
    names(tr_by_quad) <- used
    trd_all <- trd %>% filter(quad %in% used)
    cnd_p   <- cnd %>% mutate(quad = "_ctrl")

    for (design in designs) {
      key   <- sprintf("quadtest_%s_%s_%s", proxy, hyp, design)
      path  <- file.path(cfgl$cache_dir, paste0(key, ".rds"))
      stamp <- cache_stamp(cfgl, spec = c(key, spec, list(
        ids = sort(paste(trd_all$aq_id, trd_all$CntrName, trd_all$quad)),
        ens = ens, quads = used, contrasts = Call, min_treat = MIN_TREAT)))
      hit <- if (isTRUE(FORCE)) NULL else cache_read(path, stamp, cfgl)
      if (!is.null(hit)) {
        say("  [cache] ", key, "\n"); all_tabs[[key]] <- hit; next
      }

      if (design == "stacked")
        .verify_block_equivalence(tr_by_quad, cnd, sets[[1]], spec, cfgl, used)

      oplan <- future::plan(future::multisession, workers = par_workers(cfgl))
      res <- with_seed(cfgl$ensemble_seed, future_map(sets, function(ids) {
        if (design == "stacked")
          .realization_stacked(tr_by_quad, cnd, ids, spec, cfgl, used,
                               Call, which(row_kind == "omnibus_basis"))
        else
          .realization_pooled(trd_all, cnd_p, ids, spec, cfgl, used,
                              Call, which(row_kind == "omnibus_basis"))
      }, .options = furrr_options(seed = TRUE)))
      future::plan(oplan)

      n_fail <- sum(!vapply(res, function(r) isTRUE(r$ok), logical(1)))
      if (n_fail)
        sayf("    %s/%s: %d of %d realizations failed and are excluded.\n",
             hyp, design, n_fail, ens)
      pl <- .pool_design(res, Call, row_kind, used, alpha)
      if (is.null(pl)) { say("    no usable realizations.\n"); next }

      # A denominator df below ~2 means the CR2 small-sample correction has run
      # out of independent country clusters for that contrast: the p-value is
      # arithmetically defined but carries no information. Say so in the table
      # rather than leaving a reader to notice df2 = 0.9 on their own.
      tab <- pl$tab %>% mutate(
        note = dplyr::case_when(
          !is.finite(df2)               ~ "df not available",
          df2 < 2                       ~ "df < 2: not resolvable at this design",
          df2 < 5                       ~ "df < 5: read with caution",
          TRUE                          ~ ""),
        proxy = proxy, hyp = hyp, design = design,
        outcome = spec$outcome,
        quadrants = paste(used, collapse = "+"),
        n_treat_total = nrow(trd_all),
        n_treat_detail = paste(sprintf("%s:%d", used,
                                       as.integer(nq[used])), collapse = " "),
        n_country_med = pl$n_country,
        max_smd_med = paste(sprintf("%s:%.3f", used,
                                    pl$max_smd_med[used]), collapse = " "),
        n_realizations = ens, n_failed = n_fail,
        .before = 1)
      cache_write(tab, path, stamp)
      all_tabs[[key]] <- tab

      om <- tab %>% filter(tier == "CR2", quantity == "omnibus_all_equal")
      sayf("  %-4s %-8s omnibus (CR2): F(%.0f, %.1f) = %.2f, p = %s%s\n",
           hyp, design, om$df1, om$df2, om$stat, pstar(om$p),
           if (nzchar(om$note)) paste0("   [", om$note, "]") else "")
      nc <- tab %>% filter(tier == "CR2", kind == "named_contrast")
      for (i in seq_len(nrow(nc)))
        sayf("       %-6s %+.4g (SE %.4g, df %.1f)  p = %s\n",
             nc$quantity[i], nc$est[i], nc$se[i], nc$df2[i], pstar(nc$p[i]))
    }
  }

  if (!length(all_tabs)) { say("  nothing estimable.\n"); return(invisible(NULL)) }
  out <- dplyr::bind_rows(all_tabs)
  f <- file.path(cfgl$out_dir, paste0("asymmetry_quadrant_tests_", proxy, ".csv"))
  write.csv(out, f, row.names = FALSE)
  sayf("  wrote %s\n", basename(f))
  invisible(out)
}

# =============================================================================
# RUN
# =============================================================================
# The cross-quadrant tests run BEFORE the figure, because the figure now reports
# them (test_strip()). Sections 1-3 are untouched by this ordering.
rr_main <- run_proxy("agva_ha")
qt_main <- test_quadrant_differences(rr_main$dy, "agva_ha")
make_fig(rr_main, "agva_ha", "agricultural value added per ha", qt_main)

rr_alt <- run_proxy("gdp_pc")
qt_alt <- test_quadrant_differences(rr_alt$dy, "gdp_pc")
make_fig(rr_alt, "gdp_pc", "GDP per capita", qt_alt)

say("\n49_asymmetry done.\n")
