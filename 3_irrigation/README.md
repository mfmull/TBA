# Transboundary irrigation / land-use analysis — final pipeline (41-series)

Reproducible pipeline for the irrigated-cropland analysis (manuscript Fig. 3,
Fig. S8 counterparts, Table S5 replacement, balance diagnostics, and the
river-border exclusion robustness). Mirrors the structure of the wells
pipeline (`2_wells`, 31-series): one config, one core library, staged
run scripts, cached heavy objects, manifest + assertions.

## Required inputs (`1_data/`)

- `_dataMain.csv` — polygon-level master table from the current GEE build
  (929 treated TBA–country segments + 6,304 HydroSHEDS control basins;
  produced by `../1_buildData/2_buildDataset.R`).
- `CtrlNoOverlapHYBAS_B.rds` — 500 non-overlapping control realizations
  (produced by `../1_buildData/3_GenNonOverlappingCtrl.R`).
- `DataS2_meta.csv` — IGRAC metadata (names, codes, countries, region) keyed
  on `aq_id` × `CountryName`; used only to regenerate `Dataset_S3.csv` and to
  label the Lorenz examples.
- `landpct.csv`, `irpct.csv` — distance-to-border percentile curves (treated
  segments), used only for the Fig. 3C Lorenz illustration.

## Package setup

R >= 4.3 with: `dplyr`, `tidyr`, `purrr`, `tibble`, `readr`, `stringr`,
`MatchIt` (+ `optmatch`), `lme4`, `lmerTest`, `broom.mixed`, `clubSandwich`,
`future`, `furrr`, `ggplot2`, `patchwork`, `rlang`. Exact versions for the
canonical run are recorded in `output/session_info.txt`.

## One-command execution

```
Rscript 48_run_all.R
```

| script | role |
|---|---|
| `41_config.R` | every analytical choice (specs, covariate sets, best-member rule, SMD scope, frozen benchmark) |
| `42_core.R`   | core library: ensemble matching + weighted mixed models, balance extraction, CR2 refit, Lorenz curves, SI Dataset S3 build, caching |
| `43_run_main.R` | seven main specifications (Fig. 3), design assertions, frozen-benchmark check |
| `44_run_robustness.R` | SI block (Fig. S8 counterparts, GW-attribution sensitivity) + river-border exclusion block |
| `45_summaries.R` | methods numbers, sample flow, SI `Dataset_S3.csv` (segment level) |
| `46_figures.R` | canonical figures (PDF only, in `figure/`) |
| `47_tables.R`  | canonical tables (CSV + TEX, in `output/`) |
| `48_run_all.R` | orchestration, manifest, `session_info.txt`, final assertions |

Per-specification ensembles (500 matchit + lmer fits each) are cached under
`derived/cache/spec_<label>.rds` with content stamps covering the input data,
the spec definition and `42_core.R`; stale caches are refused and recomputed.
First full run: roughly 3–4 h on 2 cores (scales with cores; ~17 specs ×
500 realizations); cached re-run: minutes.

## Preferred specification

Treated units: TBA–country segments (`type == "treat"`); controls: HydroSHEDS
basins from one non-overlapping realization. Sample: `Ir > 0` polygons for
irrigation outcomes, `GW > 0` for groundwater-specific Gini outcomes. Full
optimal matching (`MatchIt`, `method = "full"`; estimand: ATT) on centroid
lat/lon, distance to major road, and CSI (intensity outcomes) or the CSI
borderward Gini (Gini outcomes). Outcome model: weighted linear mixed model
(full-matching weights; country random intercept); ensemble p-values are
Satterthwaite. Across the 500 realizations we report the coefficient
distribution and the share significant at p < 0.1. The best-balanced member
(lowest **max** |SMD|, propensity-distance row excluded) is re-fit and
reported with model-based AND CR2 country-clustered SEs (weighted-OLS route,
per MatchIt guidance; `clubSandwich` does not support weighted `lmer` fits).

## Changes vs. the legacy pipeline

See `CHANGELOG.md`. Model fits are v3-identical; only balance accounting,
diagnostics, inference add-ons, the regenerated `Dataset_S3.csv` (GW>0 count
207 → 159), and figure design changed.

## Causal vs descriptive outputs

The matched contrasts are associational estimates of the transboundary
*status* contrast (ATT under selection-on-observables within the ensemble);
figure and table titles use associational language throughout (reviewer
R1.3). Ensemble significance shares are descriptive (control-selection
uncertainty), not a multiple-testing family.

## Known limitations / dependency notes

- `clubSandwich` cannot compute CR2 on `lmer` fits with prior weights; the
  CR2 tier therefore uses the same weighted fixed-effect regression fit by
  OLS (standard post-matching practice), clustered on country.
- The river-border exclusion block drops every polygon within 50 km of a
  GSRB river border; identifying *border-perpendicular* rivers specifically
  would require additional GIS work outside this pipeline.
- The 500 realizations are overlapping draws from one candidate pool: the
  violin spread quantifies control-selection uncertainty, not sampling
  uncertainty.
