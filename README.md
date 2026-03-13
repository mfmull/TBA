# Transboundary Aquifer Analysis Pipeline

This repository contains the code used to construct spatial datasets, process groundwater and irrigation indicators, generate basin typologies, and produce figures used in the associated research project.

The analysis combines:

- local spatial processing in **R**
- satellite-derived indicators processed using **Google Earth Engine (GEE)**
- vector spatial operations on aquifer and basin geometries

Large raw datasets are **not included in the repository** due to GitHub file size limits.

## Repository Structure

```text
0_PreprocessIGRAC/
1_Map/
2_Wells/
3_Irrig/
4_Typology/

figs/
SI1.csv
SI2.csv
```

Each directory corresponds to a stage in the analysis pipeline.

## Workflow Overview

```text
RawData
   │
   ▼
0_PreprocessIGRAC
   │  Prepare aquifer geometries
   │  Split aquifers by national borders
   │  Resolve discrepancies in country membership
   ▼
2_Wells
   │  Construct well datasets
   │  Aggregate well observations to basin level
   ▼
3_Irrig
   │  Construct irrigation indicators
   │  Build matched basin samples
   ▼
4_Typology
   │  Construct basin typologies
   │  Generate classification datasets
   ▼
1_Map
   ▼
Figures used in the manuscript
```

## Script Execution

Scripts are designed to be executed **from within their own directory**.

Before running a script, set the working directory to the folder containing the script. Then run the scripts in that folder **in the order indicated by their filenames**.

Example:

```r
setwd("0_PreprocessIGRAC")
source("0_Peprocess.R")
```

Run the pipeline folders in this order:

```text
0_PreprocessIGRAC
2_Wells
3_Irrig
4_Typology
1_Map
```

Within each folder, execute the scripts in sequence.

## Raw Data

Raw datasets are **not included in the repository**.

To reproduce the analysis, create the following local directory structure:

```text
RawData/

IGRAC/
  TBA2025/
    tba_map_2025.shp

geoBoundariesCGAZ_ADM0/
  geoBoundariesCGAZ_ADM0.shp

Hydrosheds/
  hybas_lake_*.shp

GSRB/
  GSRB_Level0.shp

jasechko_aquifs/
  jasechko_et_al_2024_aquifers.shp

Fraser2023_Agreements/
  _agreements.csv
```

These inputs include aquifer geometries, country boundaries, hydrological basin layers, and other external spatial/tabular datasets used throughout the pipeline.

## Google Earth Engine Components

Parts of the irrigation-data construction rely on **Google Earth Engine (GEE)**.

The typical workflow is:

```text
Google Earth Engine
       │
       ▼
Export raster or tabular datasets
       │
       ▼
RawData/
       │
       ▼
Local processing with R scripts
```

Outputs generated from GEE should be downloaded and placed in `RawData/` before running the relevant irrigation scripts.

## Script Inventory

### 0_PreprocessIGRAC

| Script | Description |
|---|---|
| `0_Peprocess.R` | Script in `0_PreprocessIGRAC` stage of the pipeline. |

### 1_Map

| Script | Description |
|---|---|
| `Fig1.R` | Script in `1_Map` stage of the pipeline. |

### 2_Wells

| Script | Description |
|---|---|
| `1_data/getWells.R` | Script in `2_Wells` stage of the pipeline. |
| `2_firststage/FirstStages.R` | Script in `2_Wells` stage of the pipeline. |
| `3_secondstage/1_preferred.R` | Script in `2_Wells` stage of the pipeline. |
| `3_secondstage/2_SI_AltCov.R` | Script in `2_Wells` stage of the pipeline. |
| `3_secondstage/3_SI_LMRobust.R` | Script in `2_Wells` stage of the pipeline. |
| `3_secondstage/4_SI_Umatched.R` | Script in `2_Wells` stage of the pipeline. |
| `3_secondstage/5_SI_dropTopSE.R` | Script in `2_Wells` stage of the pipeline. |
| `3_secondstage/6_SI_matchArchitecture.R` | Script in `2_Wells` stage of the pipeline. |
| `3_secondstage/7_SI_nMin.R` | Script in `2_Wells` stage of the pipeline. |
| `3_secondstage/8_SI_loo.R` | Script in `2_Wells` stage of the pipeline. |
| `3_secondstage/8_SI_looCountry.R` | Script in `2_Wells` stage of the pipeline. |
| `3_secondstage/9_SI_FE.R` | Script in `2_Wells` stage of the pipeline. |
| `3_secondstage/9_SI_dclust.R` | Script in `2_Wells` stage of the pipeline. |
| `3_secondstage/9_clustRobustSE.R` | Script in `2_Wells` stage of the pipeline. |
| `4_figures/Fig2A.R` | Script in `2_Wells` stage of the pipeline. |
| `4_figures/Fig2B.R` | Script in `2_Wells` stage of the pipeline. |
| `4_figures/FigS1S2.R` | Script in `2_Wells` stage of the pipeline. |
| `4_figures/FigS3.R` | Script in `2_Wells` stage of the pipeline. |
| `4_figures/FigS4.R` | Script in `2_Wells` stage of the pipeline. |
| `4_figures/FigS5.R` | Script in `2_Wells` stage of the pipeline. |
| `4_figures/FigS6.R` | Script in `2_Wells` stage of the pipeline. |
| `4_figures/FigS7.R` | Script in `2_Wells` stage of the pipeline. |
| `4_figures/TabS1.R` | Script in `2_Wells` stage of the pipeline. |
| `4_figures/TabS2S3.R` | Script in `2_Wells` stage of the pipeline. |

### 3_Irrig

| Script | Description |
|---|---|
| `0_prepControls/0_defMatchPop.R` | Script in `3_Irrig` stage of the pipeline. |
| `0_prepControls/1_getBasins.R` | Script in `3_Irrig` stage of the pipeline. |
| `1_buildData/2_buildDataset.R` | Script in `3_Irrig` stage of the pipeline. |
| `1_buildData/3_GenNonOverlappingCtrl.R` | Script in `3_Irrig` stage of the pipeline. |
| `2_RunAnalysis/4_Analyse.R` | Script in `3_Irrig` stage of the pipeline. |
| `3_Visualize/Fig3.R` | Script in `3_Irrig` stage of the pipeline. |
| `3_Visualize/FigSIViolin.R` | Script in `3_Irrig` stage of the pipeline. |
| `3_Visualize/TableBestMatch.R` | Script in `3_Irrig` stage of the pipeline. |
| `3_Visualize/TableSummary.R` | Script in `3_Irrig` stage of the pipeline. |
| `4_Robustness/JasechkoReplication/2_buildDataset_J.R` | Script in `3_Irrig` stage of the pipeline. |
| `4_Robustness/JasechkoReplication/Analyse_J.R` | Script in `3_Irrig` stage of the pipeline. |

### 4_Typology

| Script | Description |
|---|---|
| `1_makeDyad.R` | Script in `4_Typology` stage of the pipeline. |
| `Fig4_and_mosaic.R` | Script in `4_Typology` stage of the pipeline. |
| `FigSI_200.R` | Script in `4_Typology` stage of the pipeline. |
| `FigSI_5.R` | Script in `4_Typology` stage of the pipeline. |
| `FigSI_IR.R` | Script in `4_Typology` stage of the pipeline. |
| `FigSI_eps.R` | Script in `4_Typology` stage of the pipeline. |
| `FigSI_thresh.R` | Script in `4_Typology` stage of the pipeline. |

## Software Requirements

The repository is written in **R** and depends on the packages below.

This list was compiled by scanning all uploaded `.R` scripts for explicit `library()` and `require()` calls, so it is the exhaustive set of packages referenced directly in the codebase.

- `broom.mixed`
- `clubSandwich`
- `cobalt`
- `countrycode`
- `dplyr`
- `furrr`
- `future`
- `future.apply`
- `ggalluvial`
- `ggplot2`
- `ggrepel`
- `grid`
- `Hmisc`
- `kableExtra`
- `knitr`
- `lme4`
- `lmerTest`
- `lwgeom`
- `mapview`
- `MASS`
- `MatchIt`
- `metafor`
- `modelsummary`
- `multcompView`
- `nngeo`
- `parallel`
- `patchwork`
- `performance`
- `progress`
- `progressr`
- `purrr`
- `qs`
- `raster`
- `readr`
- `rnaturalearth`
- `rnaturalearthdata`
- `sandwich`
- `scales`
- `sf`
- `splines`
- `stringr`
- `terra`
- `tibble`
- `tidyr`
- `tidyverse`
- `units`
- `vcd`

A convenient installation command is:

```r
install.packages(c(
  "broom.mixed",
  "clubSandwich",
  "cobalt",
  "countrycode",
  "dplyr",
  "furrr",
  "future",
  "future.apply",
  "ggalluvial",
  "ggplot2",
  "ggrepel",
  "Hmisc",
  "kableExtra",
  "knitr",
  "lme4",
  "lmerTest",
  "lwgeom",
  "mapview",
  "MASS",
  "MatchIt",
  "metafor",
  "modelsummary",
  "multcompView",
  "nngeo",
  "patchwork",
  "performance",
  "progress",
  "progressr",
  "purrr",
  "qs",
  "raster",
  "readr",
  "rnaturalearth",
  "rnaturalearthdata",
  "sandwich",
  "scales",
  "sf",
  "stringr",
  "terra",
  "tibble",
  "tidyr",
  "tidyverse",
  "units",
  "vcd"
))
```

Notes:

- `parallel`, `grid`, and `splines` are part of base/recommended R and typically do not need separate installation.
- Some packages may pull in overlapping dependencies, for example `tidyverse` already includes several core tidyverse packages.

## Figures and Outputs

Figures generated by the analysis pipeline are stored in:

```text
figs/
```

Supporting tables used in the manuscript and supplementary material include:

```text
SI1.csv
SI2.csv
```

## Reproducing the Analysis

1. Clone the repository.
2. Download the required raw datasets.
3. Place them in `RawData/` using the structure above.
4. Open an R session.
5. Enter each pipeline folder in the order listed above.
6. Run the scripts in that folder in filename order.

## License

This repository is distributed under the license provided in the `LICENSE` file.
