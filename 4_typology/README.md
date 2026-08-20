# 4_typology — transboundary aquifer dyad typology (50-series)

A self-contained rewrite of `4_Typology` in the same architecture as the wells
(30-series) and irrigation (40-series) pipelines. **`4_Typology` is untouched**
and remains the reference for the published outputs.

```
Rscript 58_run_all.R          # or, from RStudio:  source("58_run_all.R")
```

The whole pipeline is a deterministic classification over 568 rows: a full run
takes seconds. Nothing here needs parallelism, so none of the fork problems that
affect the wells pipeline apply.

## Layout

| file | role |
|---|---|
| `51_config.R` | every analytical choice, the robustness registry, frozen benchmarks |
| `52_core.R` | the classification rule, metrics, statistics, caching, PDF device |
| `53_run_main.R` | preferred classification + frozen-benchmark assertions |
| `54_run_robustness.R` | the five variants, as configuration overrides |
| `55_summaries.R` | methods numbers, SI2 table, sample flow |
| `56_figures.R` | Fig. 4 and the SI panels |
| `57_tables.R` | canonical tables (`output/`) |
| `58_run_all.R` | one-command pipeline, manifest, verification |
| `1_data/` | `FinalDiads.csv` and the GEE exports it was built from |

Outputs go to `output/` (tables, SI `Dataset_S4.csv`) and `figure/` (PDFs). Cached
objects live in `derived/cache/`; the run log in `derived/audit/`.

## What changed from `4_Typology`

**The rule exists once.** The legacy folder implemented the same ordered
`case_when` classification seven times — `Fig4_and_mosaic.R`, `Fig4_scatter.R`
and five `FigSI_*.R` scripts differing only in a buffer suffix, a threshold or a
variable name, each with its own hard-coded output paths. Here it is
`52_core.R::classify_dyads()`, and a robustness variant is a list of overrides:

```r
variants = list(
  near5    = list(near_zone = "_B11"),   # 5 km instead of 10 km
  far200   = list(far_zone  = "_B22"),   # 200 km instead of 100 km
  eps_lo   = list(eps       = 1e-5),     # 0.1 ha instead of 10 ha
  irrig    = list(land_var  = "IR"),     # all irrigation, not groundwater
  thresh50 = list(thresh    = 0.50))     # 50% instead of 99%
```

**Four defects from the audit are fixed.**

1. *NA facet strip* — `FigSI_5/IR/200.R` named the first metric `FIR` and then
   applied `factor(Variable, c("FGW","FBW"))`, turning every row of that facet
   into `NA`. S11–S13 carry an "NA" strip as a result. Metric names are now
   fixed and only the axis *label* varies, so it cannot recur.
2. *SE denominator* — `sd(x)/sqrt(n())` counted rows where the value was `NA`.
   Now `sqrt(sum(!is.na(x)))`.
3. *Tukey level* — the code used 90% while the Fig. 4 caption said 95%. Now
   95%, matching the caption. The letter groupings are identical at both
   levels, so the figure does not change.
4. *Sparse-cell chi-square* — 11 of 20 cells in the class × cooperation table
   have expected counts below 5. Monte Carlo p-values are used, and the cell
   sparsity is recorded in `table_si_contingency_tests.csv`.

**Additions the audit recommended**, each one setting away in `51_config.R`:
a Kruskal–Wallis check alongside the ANOVA (`report_kruskal`), and a
country-clustered bootstrap of the class means (`cluster_boot`) — dyads are not
independent, since one country appears in up to 39 of them.

**Fig. 4B/C is now the scatter** that reviewer R1.7 asked for — national
relevance against overdraft exposure, class centroids with confidence intervals
and Tukey letters on both axes. The bar panels move to
`figure_si_typology_bars.pdf` so the letters stay reported.

## Verification

`53` halts unless the classification reproduces the frozen benchmark: 568 dyads,
142 countries, 409 aquifers, and

| | IU0 | IU1 | UB | DA | BL |
|---|---|---|---|---|---|
| current | 85 | 118 | 207 | 113 | 45 |
| future | 70 | 95 | 99 | 204 | 100 |

`54` asserts the same for every variant, and `58` re-checks both plus the
exported tables. The benchmarks were taken from an independent reimplementation
that reproduces `4_Typology/SI2.csv` 568/568 on both `class` and `classF`.

`58` also checks that each canonical output exists **and is non-trivial in
size** — existence alone is not enough, because a failed PDF device writes
nothing while reporting success.

## Two things to know

**The dyad table is frozen input.** `1_data/FinalDiads.csv` is shipped as-is.
Rebuilding it from geometry needs the IGRAC aquifer split, the GSRB river layer
and the Fraser agreement table; the latter two are not in the repository, and
the legacy `1_makeDyad.R` also pointed at a misspelled `0_PreporcessIGRAC`.
Point `rivers_shapefile` and `agreements_file` in `51_config.R` at real files to
enable the rebuild.

**The published palette has a legibility problem.** Checked with a
colour-vision validator, IU0 `#ABABAB` and IU1 `#A6CEE3` separate by only
ΔE 10.2 for *normal* vision (15 is the usual floor) and 9.4 under deuteranopia —
two adjacent categories that are hard to tell apart for everyone. The published
palette is the default so Fig. 4 reproduces as printed; `palette = "cvd"` in
`51_config.R` switches to a validated alternative that passes for normal,
deuteranopic, protanopic and tritanopic vision.
