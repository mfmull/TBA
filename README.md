# Geography of groundwater use and depletion in internationally shared aquifers

Code and processed data for:

> Müller, M.F., Davis, K.F., Faure, J., Chiarelli, D.D., & Hung, F.
> *Geography of groundwater use and depletion in internationally shared aquifers.*
> (in revision, 2026)

Contact: Marc F. Müller (marc.mueller@eawag.ch)

---

## Repository structure

```
0_preprocess_igrac/     IGRAC transboundary-aquifer polygons intersected with
                        national boundaries -> the 990-row segment inventory
1_gee/                  Google Earth Engine scripts (JavaScript) that extract
                        every land-use, irrigation, well and border-distance
                        layer the downstream analyses consume
2_wells/                Well-based depletion analysis (H1/H2): overlap-weighted
                        two-stage meta-regression, specification registry,
                        heterogeneity/BLUP analysis and the selection rule
3_irrigation/           Irrigated-cropland analysis: 500-realization matched
                        control ensemble, borderward Gini indices, river-border
                        adjustment and exclusion, economic-asymmetry analysis
4_typology/             Dyad typology: ordered classification, dyad-border
                        correction, external-audit corrections (Stages 1-3),
                        robustness registry
supplementary_datasets/ Datasets S1-S5 exactly as cited in the paper, with a
                        column-by-column CODEBOOK.md
assemble_SI_datasets.R  The single place a pipeline output becomes a numbered
                        SI dataset. Run last.
```

Each analysis module is self-contained: `NN_config.R` holds every path and
parameter, `NN_core.R` the functions, then numbered driver scripts, and
`NN_run_all.R` regenerates the module end to end. Each `output/` folder holds
the exact tables shipped with the paper plus an `output_manifest.csv` and a
`session_info.txt` recording the run that produced them.

## Reproducing the paper

```bash
# 1. Extract the raster/vector layers (Google Earth Engine account required)
#    Run each script in 1_gee/ and export the CSVs it defines.

Rscript 0_preprocess_igrac/0_preprocess.R   # segment inventory

# 2. Optional data preparation for two robustness rows of SI Table S4.
#    Both are slow, need large third-party downloads, and self-skip if absent.
Rscript 2_wells/classify_river_crossings.R  # GloRiC border-crossing rivers
Rscript 2_wells/get_wtd_fan.R               # Fan et al. natural water-table depth

# 3. The three analyses
Rscript 2_wells/38_run_all.R        # MS Fig. 2; SI Figs. S1-S2; SI Tables S2-S6
Rscript 3_irrigation/48_run_all.R   # MS Fig. 3; SI Figs. S3, S5, S6; Tables S7-S9
Rscript 3_irrigation/49_asymmetry.R # SI Fig. S4
Rscript 4_typology/58_run_all.R     # MS Fig. 4; SI Fig. S7; Tables S10-S13

# 4. Numbered supplementary datasets
Rscript assemble_SI_datasets.R      # -> supplementary_datasets/S1-S5 + MANIFEST
```

`FORCE=1` rebuilds every cache; without it each module reuses cached
intermediates whose stamps still match the code and configuration.

## What produces what

| Paper artefact | Produced by |
|---|---|
| Fig. 1 global map | `1_gee/getMainMap.txt`, `1_gee/getMainMapData.txt` |
| Fig. 2 A-C | `2_wells/36_figures.R` → `figure_main_three_panel.pdf` |
| Fig. 2 D well-density map | `2_wells/36_figures.R` → `figure_main_wells_map.pdf` |
| Fig. 2 E segment insets | `2_wells/36_figures.R` → `figure_main_wells_insets.pdf` |
| Fig. 3 irrigation ensemble | `3_irrigation/46_figures.R` |
| Fig. 4 typology | `4_typology/56_figures.R` |
| Fig. S1 influence analyses | `2_wells/36_figures.R` → `figure_si_influence.pdf` |
| Fig. S2 hydraulic connectivity | `2_wells/34_run_robustness.R` (estimation) + `36_figures.R` (plot) |
| Fig. S3 alternative outcomes | `3_irrigation/46_figures.R` |
| Fig. S4 economic asymmetry | `3_irrigation/49_asymmetry.R` |
| Fig. S5 matching balance | `3_irrigation/46_figures.R` |
| Fig. S6 river-border exclusion | `3_irrigation/46_figures.R` |
| Fig. S7 cooperation mosaics | `4_typology/56_figures.R` |
| Table S2 balance / overlap designs | `2_wells/37_tables.R` |
| Table S3 main wells contrasts | `2_wells/37_tables.R` → `table_main_results.csv` |
| Table S4 specification registry | `2_wells/34_run_robustness.R` + `37_tables.R` |
| Table S5 governance proxy | `2_wells/1_data/governance_proxy.csv` (an input; see below) |
| Table S6 segment stability | `2_wells/37_tables.R` → `table_si_aquifer_stability.csv` |
| Tables S7-S9 irrigation | `3_irrigation/47_tables.R` |
| Table S10 border zones | `2_wells/1_data/IGRAC_Properties.csv` via `34_run_robustness.R` |
| Tables S11-S13 typology | `4_typology/57_tables.R` |
| Datasets S1-S5 | `assemble_SI_datasets.R` |

**Two artefacts are inputs, not pipeline products.** `governance_proxy.csv`
(SI Table S5) is an ordinal 0-2 coding of national groundwater governance,
hand-coded from national statutes with a source recorded per country; it
cannot be derived from any public file and is therefore shipped here.
`IGRAC_Properties.csv` carries the transmissivity and storativity values
behind the border-zone radii of SI Table S10.

**One SI table is not reproducible from this code.** SI Table S1 (well-analysis
summary statistics) predates the current sample: it reports N = 637 domestic
segments against the 654 used everywhere else, and three of its columns are not
in the current configuration. No script in this repository regenerates it.

## Supplementary datasets

| file | rows | contents |
|---|---|---|
| `Dataset_S1.csv` | 44 | Monitored transboundary aquifer–country segments: well counts, border-distance shares, raw and shrunken excess depletion and borderward gradient, ranks, bootstrap and robustness selection frequencies |
| `Dataset_S2.csv` | 990 | Every IGRAC aquifer–country intersection with the visual-consistency verdict (929 retained for analysis) |
| `Dataset_S3.csv` | 929 | The transboundary segments of the irrigation analysis: every outcome, sample filter and matching/regression covariate the 41-series specifications use |
| `Dataset_S4.csv` | 568 | The dyads: typology class under current and future-demand conditions, SDG 6.5.2 cooperation level, river border, distance class, and the stage-by-stage build with the note explaining every reclassification |
| `Dataset_S5.csv` | 568 | External-audit census of all dyad classifications (evidence, verdicts, geometric validation) |

Column definitions are in `supplementary_datasets/CODEBOOK.md`.
`assemble_SI_datasets.R` asserts row counts and key uniqueness on every build,
so a dataset cannot silently change length.

## Input data (not redistributed)

Raw third-party inputs are excluded from this repository. All are publicly
available from the sources below; their preparation is documented in `1_gee/`
and in each module's `NN_config.R`.

| Input | Source | DOI / URL |
|---|---|---|
| Groundwater-level trends; aquifer delineations | Jasechko et al. (2024), *Nature* 625, 715–721 | https://doi.org/10.1038/s41586-023-06879-8 |
| Transboundary aquifer polygons | IGRAC, *Transboundary Aquifers of the World*, Sept 2025 update | https://www.un-igrac.org |
| Irrigated areas | Meier et al. (2018), *HESS* 22, 1119–1133 | https://doi.org/10.5194/hess-22-1119-2018 |
| Groundwater share of irrigation | Siebert et al. (2010), *HESS* 14, 1863–1880 | https://doi.org/10.5194/hess-14-1863-2010 |
| Cropland extent | Potapov et al. (2022), *Nature Food* 3, 19–28 | https://doi.org/10.1038/s43016-021-00429-z |
| Green/blue water scarcity (WATNEEDS) | Chiarelli et al. (2020), *Sci. Data* 7, 273 | https://doi.org/10.1038/s41597-020-00612-0 |
| Blue-water availability | PCR-GLOBWB 2, Sutanudjaja et al. (2018) | https://doi.org/10.5194/gmd-11-2429-2018 |
| Settlements | European Commission JRC, GHS-SMOD R2023A | https://human-settlement.emergency.copernicus.eu |
| Hydrological basins | HydroSHEDS, Lehner & Grill (2013) | https://www.hydrosheds.org |
| River network (border-crossing rivers) | GloRiC v1.0, Ouellet Dallaire et al. (2019), *ERL* 14, 024003 | https://doi.org/10.1088/1748-9326/aad8e9 |
| River-defined borders | GSRB, Popelka & Smith (2020) | https://doi.org/10.1080/1755876X.2020.1747183 |
| Country boundaries | GADM; geoBoundaries CGAZ ADM0 | https://gadm.org · https://www.geoboundaries.org |
| Crop suitability | FAO/IIASA GAEZ v4 | https://gaez.fao.org |
| Natural (pre-development) water-table depth | Fan et al. (2013), *Science* 339, 940–943 | https://doi.org/10.1126/science.1229881 |
| Maximum economic pumping depth | Bierkens et al. (2022), *ERL* 17, 104008 | https://doi.org/10.1088/1748-9326/ac8de1 |
| Transboundary cooperation (SDG 6.5.2) | Fraser et al. (2023), *Water Policy* 25, 736–753 | https://doi.org/10.2166/wp.2023.107 |
| Formal-arrangement corrections | Burchi (2018), *J. Hydrol. Reg. Stud.* 20, 15–20 | https://doi.org/10.1016/j.ejrh.2018.04.007 |
| Fuel prices | GlobalPetrolPrices via TheGlobalEconomy.com | https://www.theglobaleconomy.com |
| Agricultural value added, agricultural land | World Bank World Development Indicators | https://databank.worldbank.org |

The map figures of Fig. 2 D–E additionally need the Jasechko aquifer polygons
and geoBoundaries ADM0 shapefiles under `2_wells/map/`; that block stops with
an explicit message when they are absent, and the rest of the module runs.

## Requirements

R ≥ 4.3 with `tidyverse`, `sf`, `metafor`, `MatchIt`, `optmatch`, `lme4`,
`clubSandwich`, `sandwich`, `fwildclusterboot`, `vcd`, `modelsummary`,
`tinytable`, `patchwork`, `maps`, `ncdf4` (only for `get_wtd_fan.R`).
A Google Earth Engine account is needed for `1_gee/`.

## Citation

See `CITATION.cff`. 

## License

MIT — see `LICENSE`.
