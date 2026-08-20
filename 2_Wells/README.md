# 2_wells — well-based groundwater-depletion analysis (30-series)

Overlap-weighted two-stage meta-regression comparing groundwater-depletion
trends in transboundary aquifer–country segments with matched domestic
segments, plus the specification registry, influence analyses and the
heterogeneity/selection machinery behind MS Fig. 2.

## Scripts

| script | role |
|---|---|
| `31_config.R` | every path, parameter, frozen reference value and tolerance |
| `32_core.R` | the function library (`run_specification`, `fit_stage2`, caching, Conley/WCB inference, `h3_analysis`) |
| `33_run_main.R` | preferred specification, balance, MDEs, BLUPs, the selection rule |
| `34_run_robustness.R` | the specification registry, leave-one-out and influence sweeps, **and the augmented specifications** (below) |
| `35_summaries.R` | derived summaries consumed by figures and tables |
| `36_figures.R` | MS Fig. 2 A–E, SI Fig. S1 and SI Fig. S2 |
| `37_tables.R` | main and SI tables, plus `methods_numbers` |
| `38_run_all.R` | runs 33→37, then the output manifest and the verification assertions |
| `make_table_s1.R` | standalone entry point for SI Table S1; needs only `main_objects.rds` from 33, so it does not wait on the robustness blocks |

Two standalone data-preparation helpers, both slow and both optional:

- `classify_river_crossings.R` — classifies which shared aquifers have their
  two national parts joined by a border-crossing river, from the GloRiC v1.0
  network. Criteria: reach discharge ≥ 10 m³ s⁻¹, ≥ 5 km of channel inside the
  segment polygon, and ≥ 10 km of penetration beyond the shared border on both
  sides. Writes `1_data/river_crossings_by_unit.csv` plus a Q ∈ {10, 30, 100,
  300} sweep. All five files ship with this repository.
- `get_wtd_fan.R` — downloads the Fan et al. (2013) natural water-table depth
  grids and samples them at well locations, writing `1_data/wtd_by_well.csv`.
  That file is several MB and is **not** shipped; the corresponding robustness
  row self-skips without it.

## The augmented specifications

Four rows of SI Table S4 refit the design on a segment frame that the registry
cannot express, because the added covariate does not live in `wellsData.csv`.
They are estimated in `34_run_robustness.R` (section
"Augmented specifications and hydraulic-connectivity subsets") and appended to
the registry by `37_tables.R`:

| row | needs |
|---|---|
| Well-level regression, random intercepts | `lme4` |
| Well-level regression, country fixed effects | `clubSandwich`, `sandwich` |
| Augmented covariate set (+ governance proxy) | `1_data/governance_proxy.csv` |
| Augmented covariate set (+ natural water-table depth) | `1_data/wtd_by_well.csv` |

The same section also builds the hydraulic-connectivity subsets behind SI
Fig. S2, using transmissivity and storativity from
`1_data/IGRAC_Properties.csv`, and the crossing-river exclusion row and its
complement (SI Table S4, note c).

Each block self-skips with a message when its optional input or package is
absent, so the canonical paper outputs never depend on something a user has
not downloaded.

> These blocks were previously carried in separate reviewer-response scripts
> (`39_responses.R`, `39b_responses_addons.R`) writing to a `ResponsesSupport/`
> folder. Both files and that folder are gone: everything the paper or SI uses
> now lives in the numbered pipeline and writes to `output/`, and everything
> else has been deleted.

## Inputs

Shipped here: `governance_proxy.csv` (ordinal 0–2 national groundwater
governance, hand-coded from national statutes with a per-country source —
this is also the content of SI Table S5 and cannot be regenerated),
`IGRAC_Properties.csv` (transmissivity/storativity behind SI Table S10), and
the `river_crossings_*.csv` classification.

Not shipped, see the data table in the top-level README: `wellsData.csv`
(Jasechko et al. 2024 well trends), `wells_points.csv`, `wtd_by_well.csv`,
and the shapefiles under `map/` needed by the Fig. 2 D–E maps.

## Caching

Every expensive object is cached under the module's cache directory and
stamped with the code and configuration it was built from, so a rerun only
recomputes what actually changed. `FORCE=1` rebuilds everything;
`FORCE_RATIOS=1` rebuilds only the Conley inflation ratios, moving the
existing table aside rather than deleting it.

## Verification

`38_run_all.R` finishes by asserting that the figures and tables agree on the
main results and that the preferred estimates still match the frozen reference
values in `31_config.R` to within `reproduce_tol`. A run that changes a
headline number fails loudly.
