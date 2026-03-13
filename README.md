
# Replication code for: Border-proximate irrigation and groundwater depletion in internationally shared aquifers

This repository contains the code used in the manuscript "Border-proximate irrigation and groundwater depletion in internationally shared aquifers" by Muller, Faure, Davis, Chiarelli and Hung (2025). 

The analysis combines:

- Local spatial processing in **R**
- Satellite-derived indicators processed using **Google Earth Engine (GEE)**
- Vector spatial operations on aquifer and basin geometries

Large raw datasets are **not included in the repository** due to GitHub file size limits.

---

# Repository Structure

```
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

---

# Workflow Overview

```
RawData
   │
   ▼
0_PreprocessIGRAC
   ▼
2_Wells
   ▼
3_Irrig
   ▼
4_Typology
   ▼
1_Map
   ▼
Figures
```

---

# Running the Analysis

Scripts are designed to be executed **from within their own directory**.

Example:

```r
setwd("0_PreprocessIGRAC")
source("script_name.R")
```

Run scripts **in each folder in the order indicated by their filenames**.

Pipeline order:

```
0_PreprocessIGRAC → 2_Wells → 3_Irrig → 4_Typology → 1_Map
```

---

# Raw Data

Raw datasets must be placed locally in:

```
RawData/
```

Example structure:

```
RawData/
  IGRAC/TBA2025/
  geoBoundariesCGAZ_ADM0/
  Hydrosheds/
  GSRB/
```

These files are not included in the repository.

---

# Script Inventory

The following scripts are included in the repository. Where relevant, the table also lists the Google Earth Engine (GEE) script(s) that need to be run beforehand or that provide the mapped/exported inputs used by the corresponding R script.

### 0_PreprocessIGRAC

| Script | Description | GEE links |
|---|---|---|
| `0_Peprocess.R` | Splits IGRAC transboundary aquifers by country boundaries, flags and resolves country-membership discrepancies, and exports cleaned aquifer-country polygons plus supplementary outputs used downstream. | — |

### 1_Map

| Script | Description | GEE links |
|---|---|---|
| `Fig1.R` | Builds Figure 1, a stacked bar chart that decomposes irrigated area from GEE outputs into surface-water vs. groundwater irrigation and non-TBA vs. TBA groundwater irrigation components. | Main figure map: `https://code.earthengine.google.com/fdf7e7e39ec8e3450c30f3ce5ea89ba2`  •  Data export: `https://code.earthengine.google.com/a83f29a381d3b721f5ab476ab5a06500?noload=1` |

### 2_Wells

| Script | Description | GEE links |
|---|---|---|
| `getWells.R` | Assembles the cleaned, spatially enriched well-level dataset from GEE exports; links wells to countries, nearest land borders, and aquifers; and exports both well-level and aquifer-level summary tables. | Wells export: `https://code.earthengine.google.com/6ed1434fbca51b7acaeb2658f048d5e4?noload=1`  •  Well–aquifer export: `https://code.earthengine.google.com/9661594f3ec29ae3f375764555334fea?noload=1` |
| `FirstStages.R` | Estimates first-stage well regressions for each aquifer-country unit, including baseline OLS, robust-regression, and spatial-declustering variants, and exports unit-level coefficients. | — |
| `1_preferred.R` | Runs the preferred second-stage matched meta-regression comparing transboundary and non-transboundary aquifer-country units. | — |
| `2_SI_AltCov.R` | Repeats the preferred second-stage analysis using an alternative set of matching covariates. | — |
| `3_SI_LMRobust.R` | Repeats the second-stage analysis using first-stage estimates derived from robust regression instead of OLS. | — |
| `4_SI_Umatched.R` | Estimates pooled second-stage effects without any covariate matching. | — |
| `5_SI_dropTopSE.R` | Re-runs the preferred second-stage analysis after excluding predefined high-leverage or top-heavy aquifers. | — |
| `6_SI_matchArchitecture.R` | Tests sensitivity of second-stage results to matching architecture, including nearest-neighbor matching and alternative control-to-treated ratios. | — |
| `7_SI_nMin.R` | Tests sensitivity of second-stage results to the minimum wells-per-unit threshold used to retain aquifer-country units. | — |
| `8_SI_loo.R` | Performs a leave-`k`-out robustness exercise by repeatedly dropping subsets of transboundary aquifers and re-estimating pooled effects. | — |
| `8_SI_looCountry.R` | Performs an overlap-country leave-one-out robustness check for the preferred second-stage analysis. | — |
| `9_SI_FE.R` | Re-estimates the second-stage model using country fixed effects instead of the multilevel random-effects specification. | — |
| `9_SI_dclust.R` | Summarizes how second-stage results change across alternative spatial declustering settings used in the first stage. | — |
| `9_clustRobustSE.R` | Re-estimates the preferred second-stage model while reporting country-cluster-robust inference for fixed effects. | — |
| `Fig2A.R` | Produces map-based visualizations of first-stage aquifer-country estimates and interactive close-up maps of well-level groundwater depletion patterns. | — |
| `Fig2B.R` | Creates the compact coefficient-style visualization used for the preferred second-stage meta-regression results. | — |
| `FigS1S2.R` | Produces matching diagnostics and covariate-balance figures for the second-stage analysis. | — |
| `FigS3.R` | Visualizes sensitivity of second-stage estimates to the nearest-neighbor matching ratio. | — |
| `FigS4.R` | Visualizes sensitivity of second-stage estimates to the minimum wells-per-unit threshold (`nMin`). | — |
| `FigS5.R` | Visualizes sensitivity of second-stage estimates to the leave-`k`-out transboundary-aquifer exercise. | — |
| `FigS6.R` | Visualizes how pooled second-stage estimates vary with the spatial declustering radius used in first-stage estimation. | — |
| `FigS7.R` | Visualizes the overlap-country leave-one-out results and summarizes how much identifying overlap is concentrated in countries containing both treated and control segments. | — |
| `TabS1.R` | Builds supplementary summary tables from the first-stage data, including matched-weighted comparisons between transboundary and non-transboundary units. | — |
| `TabS2S3.R` | Builds publication-ready regression tables comparing the preferred second-stage results with robustness and alternative specifications. | — |

### 3_Irrig

| Script | Description | GEE links |
|---|---|---|
| `0_defMatchPop.R` | Identifies candidate HydroSHEDS control basins that resemble each transboundary aquifer across basin levels and exports candidate-match files. | — |
| `1_getBasins.R` | Builds the control-basin polygon dataset matched to IGRAC transboundary aquifers for downstream analysis and for ingestion into GEE. | — |
| `2_buildDataset.R` | Builds the master irrigation analysis table by combining GEE-derived summaries for treated aquifers and matched control basins, including means tables and border-distance percentile curves. | Control basins export: `https://code.earthengine.google.com/4ea83dc49991c19cfdaea804b7f4de35?noload=1`  •  TBA export: `https://code.earthengine.google.com/74de9d6141ba8012ed146a1fe5615c64?noload=1` |
| `2_buildDataset_J.R` | Builds the irrigation analysis table for the Jasechko split-aquifer robustness exercise using GEE-exported summaries. | Jasechko robustness export: `https://code.earthengine.google.com/8ec4e37cf05eb8376a7416e2f098fc5d?noload=1` |
| `3_GenNonOverlappingCtrl.R` | Generates repeated random sets of non-overlapping HydroSHEDS control basins from the candidate-basin pool. | — |
| `4_Analyse.R` | Runs the main irrigation matching pipeline across many candidate control sets, estimates mixed-effects outcome models, and saves result objects for downstream comparison. | — |
| `Analyse_J.R` | Runs the irrigation matching and mixed-effects analysis on the Jasechko-based robustness dataset and produces the corresponding regression table. | — |
| `Fig3.R` | Produces the main irrigation figure showing mirrored density distributions of estimated treatment effects across candidate control sets. | Associated map / GEE reference: `https://code.earthengine.google.com/508d15eb659bdb3f4e44a4bdc3ab93cd` |
| `FigSIViolin.R` | Produces supplementary mirrored-density figures for additional irrigation outcomes and robustness specifications. | Associated map / GEE reference: `https://code.earthengine.google.com/508d15eb659bdb3f4e44a4bdc3ab93cd` |
| `TableBestMatch.R` | Generates the publication-ready regression table based on the top-ranked irrigation models from the matching analysis. | — |
| `TableSummary.R` | Generates the summary-statistics table for the aquifer-level irrigation dataset. | — |

### 4_Typology

| Script | Description | GEE links |
|---|---|---|
| `1_makeDyad.R` | Builds transboundary-aquifer dyads, attaches agreement scores and GEE-derived indicators for each side of the dyad, and exports the final dyad dataset used by the typology figures. | Dyad metrics export: `https://code.earthengine.google.com/e32121007cd578cd8c741d9a7815c02f?noload=1` |
| `Fig4_and_mosaic.R` | Produces the main dyad-typology figure by classifying aquifer dyads according to border-proximate irrigation patterns and summarizing the resulting classes graphically. | — |
| `FigSI_200.R` | Repeats the dyad-typology analysis using the 200 km buffer-based inputs as a sensitivity check. | — |
| `FigSI_5.R` | Repeats the dyad-typology analysis using the 5 km buffer-based inputs as a sensitivity check. | — |
| `FigSI_IR.R` | Repeats the dyad-typology analysis using the alternative IR-based classification setup as a sensitivity check. | — |
| `FigSI_eps.R` | Tests sensitivity of the dyad-typology classification to the epsilon parameter used in class assignment. | — |
| `FigSI_thresh.R` | Tests sensitivity of the dyad-typology classification to the threshold parameter used in class assignment. | — |

---

# Software Requirements

Main R packages used:

```
broom.mixed
clubSandwich
cobalt
countrycode
dplyr
furrr
future
future.apply
ggalluvial
ggplot2
ggrepel
Hmisc
kableExtra
knitr
lme4
lmerTest
lwgeom
mapview
MASS
MatchIt
metafor
modelsummary
multcompView
nngeo
patchwork
performance
progress
progressr
purrr
qs
raster
readr
rnaturalearth
rnaturalearthdata
sandwich
scales
sf
stringr
terra
tibble
tidyr
tidyverse
units
vcd
```

Install using:

```r
install.packages(c(
"broom.mixed","clubSandwich","cobalt","countrycode","dplyr","furrr","future",
"future.apply","ggalluvial","ggplot2","ggrepel","Hmisc","kableExtra","knitr",
"lme4","lmerTest","lwgeom","mapview","MASS","MatchIt","metafor","modelsummary",
"multcompView","nngeo","patchwork","performance","progress","progressr",
"purrr","qs","raster","readr","rnaturalearth","rnaturalearthdata","sandwich",
"scales","sf","stringr","terra","tibble","tidyr","tidyverse","units","vcd"
))
```

---

# License

This repository is distributed under the license provided in the `LICENSE` file.
