
# Replication code for: Border-proximate irrigation and groundwater depletion in internationally shared aquifers

This repository contains the code used in the manuscript "Border-proximate irrigation and groundwater depletion in internationally shared aquifers" by Muller, Faure, Davis, Chiarelli and Hung (2026). 

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
1_Map
   ▼
2_Wells
   ▼
3_Irrig
   ▼
4_Typology
   ▼
Figures
```

---

# Running the Analysis

Scripts are designed to be executed **from within their own directory**.

Example:

```r
setwd("0_PreprocessIGRAC")
source("0_Preprocess.R")
```

Run scripts **in each folder in the order indicated by their filenames**.

Pipeline order:

```
0_PreprocessIGRAC → 2_Wells → 3_Irrig → 4_Typology → 1_Map
```

---

# Input Data

Raw datasets are **not included in the repository** due to size constraints and licensing restrictions. To reproduce the full preprocessing workflow, the required raw data should be placed locally in the directory:
```
RawData/
```

A typical directory structure is:
```
RawData/
  IGRAC/TBA2025/
  geoBoundariesCGAZ_ADM0/
  Hydrosheds/
  GSRB/
```
These datasets are used only in the **preprocessing stage** of the pipeline.

Most analysis scripts operate on **processed datasets that are already included in this repository**, so running the full preprocessing workflow is **not required to reproduce the main results**. Scripts that rely on raw datasets are marked with an asterisk `*` in the script inventory below.

## Google Earth Engine preprocessing

Some large-scale raster reduction operations were performed using **Google Earth Engine (GEE)**. These operations generate intermediate datasets used by several scripts in the irrigation and wells analyses.

Outputs from these Earth Engine workflows are included in this repository, typically within directories named:
```
geeOut/
```
Each R script that relies on Earth Engine outputs includes a **static link to the corresponding Earth Engine script in its preamble**. Running these scripts requires:

- access to a web browser  
- a valid Google Earth Engine account

However, executing the Earth Engine scripts is **not necessary to reproduce the main analyses**, because the exported Earth Engine outputs used in the paper are already included in the repository.

---
# Script Inventory

The following scripts are included in the repository.

### 0_PreprocessIGRAC

| Script | Description |
|------|-------------|
| *`0_Peprocess.R` | Splits IGRAC transboundary aquifers by country boundaries, flags and resolves country-membership discrepancies, and exports cleaned aquifer-country polygons plus supporting SI outputs. |
| *`OverlapRatios_GEECode.rtf` | Links to GEE codes to carry out overlap analaysis, i.e. to determine the fraction irrigated cropland that covers (i) distinct IGRAC aquifers that overlap spatially and (ii) IGRAC aquifers where the depth to GW is not econmically viable|

### 1_Map

| Script | Description |
|------|-------------|
| `Fig1.R` | Builds Figure 1, a stacked bar chart that decomposes irrigated area from GEE outputs into surface-water vs groundwater irrigation and TBA vs non-TBA groundwater components. |

### 2_Wells

| Script | Description |
|------|-------------|
| *`getWells.R` | Assembles the cleaned well-level dataset from GEE exports, links wells to countries, nearest borders, and aquifers, and exports well- and aquifer-level summaries. |
| `FirstStages.R` | Estimates first-stage well regressions by aquifer-country unit under baseline OLS, robust-regression, and spatial-declustering specifications, then exports unit-level coefficients. |
| `1_preferred.R` | Runs the preferred second-stage matched meta-regression comparing transboundary and non-transboundary aquifer-country units. |
| `2_SI_AltCov.R` | Repeats the preferred second-stage analysis using an alternative set of matching covariates. |
| `3_SI_LMRobust.R` | Repeats the second-stage analysis using first-stage estimates derived from robust regression rather than OLS. |
| `4_SI_Umatched.R` | Estimates pooled second-stage effects without any covariate matching. |
| `5_SI_dropTopSE.R` | Re-runs the preferred second-stage analysis after excluding predefined high-leverage / high-uncertainty aquifers. |
| `6_SI_matchArchitecture.R` | Tests sensitivity of second-stage results to the matching architecture, including nearest-neighbor matching and different control-to-treated ratios. |
| `7_SI_nMin.R` | Tests sensitivity of second-stage results to the minimum wells-per-unit threshold used to retain aquifer-country units. |
| `8_SI_loo.R` | Performs a leave-k-out robustness exercise by repeatedly dropping subsets of transboundary aquifers and re-estimating pooled effects. |
| `8_SI_looCountry.R` | Performs an overlap-country leave-one-out robustness check for the preferred second-stage analysis. |
| `9_SI_FE.R` | Re-estimates the second-stage model using country fixed effects instead of the multilevel random-effects specification. |
| `9_SI_dclust.R` | Summarizes how second-stage results change across alternative spatial declustering settings used in the first stage. |
| `9_clustRobustSE.R` | Re-estimates the preferred second-stage model while reporting country-cluster-robust inference for fixed effects. |
| `Fig2A.R` | Produces map-based visualizations of first-stage aquifer-country estimates and interactive close-up maps of well-level groundwater depletion patterns. |
| `Fig2B.R` | Creates the compact coefficient-style visualization used for the preferred second-stage meta-regression results. |
| `FigS1S2.R` | Produces matching diagnostics and covariate-balance figures for the second-stage analysis. |
| `FigS3.R` | Visualizes sensitivity of second-stage estimates to the nearest-neighbor matching ratio. |
| `FigS4.R` | Visualizes sensitivity of second-stage estimates to the minimum wells-per-unit threshold (`nMin`). |
| `FigS5.R` | Visualizes sensitivity of second-stage estimates to the leave-k-out transboundary-aquifer exercise. |
| `FigS6.R` | Visualizes how pooled second-stage estimates vary with the spatial declustering radius used in first-stage estimation. |
| `FigS7.R` | Visualizes the overlap-country leave-one-out results and summarizes how much identifying overlap is concentrated in countries containing both treated and control segments. |
| `TabS1.R` | Builds supplementary summary tables from the first-stage data, including matched-weighted comparisons between transboundary and non-transboundary units. |
| `TabS2S3.R` | Builds publication-ready regression tables comparing the preferred second-stage results with robustness and alternative specifications. |

### 3_Irrig

| Script | Description |
|------|-------------|
| *`0_defMatchPop.R` | Identifies candidate HydroSHEDS control basins that resemble each transboundary aquifer across basin levels and exports candidate-match files. |
| *`1_getBasins.R` | Builds the control-basin polygon dataset matched to IGRAC transboundary aquifers for downstream analysis and GEE processing. |
| *`2_buildDataset.R` | Builds the master irrigation analysis table by combining GEE-derived summaries for treated aquifers and matched control basins. |
| *`2_buildDataset_J.R` | Builds the corresponding irrigation analysis table using the Jasechko split-aquifer dataset for robustness / replication exercises. |
| `3_GenNonOverlappingCtrl.R` | Generates repeated random sets of non-overlapping HydroSHEDS control basins from the candidate basin pool. |
| `4_Analyse.R` | Runs the main irrigation matching pipeline across many candidate control sets, estimates mixed-effects outcome models, and saves result objects for downstream comparison. |
| `Analyse_J.R` | Runs the irrigation matching and mixed-effects analysis on the Jasechko-based robustness dataset and produces the corresponding regression table. |
| `Fig3.R` | Produces the main irrigation figure showing mirrored density distributions of estimated treatment effects across candidate control sets. |
| `FigSIViolin.R` | Produces supplementary mirrored-density figures for additional irrigation outcomes and robustness specifications. |
| `TableBestMatch.R` | Generates the LaTeX regression table based on the top-ranked irrigation models from the matching analysis. |
| `TableSummary.R` | Generates the summary-statistics table for the aquifer-level irrigation dataset. |

### 4_Typology

| Script | Description |
|------|-------------|
| *`1_makeDyad.R` | Builds transboundary aquifer dyads, attaches agreement scores and GEE-derived indicators for each side of the dyad, and exports the final dyad dataset. |
| `Fig4_and_mosaic.R` | Produces the main dyad-typology figure by classifying aquifer dyads according to border-proximate irrigation patterns and summarizing the resulting classes graphically. |
| `FigSI_200.R` | Repeats the dyad-typology analysis using the 200 km buffer-based inputs as a sensitivity check. |
| `FigSI_5.R` | Repeats the dyad-typology analysis using the 5 km buffer-based inputs as a sensitivity check. |
| `FigSI_IR.R` | Repeats the dyad-typology analysis using the alternative IR-based classification setup as a sensitivity check. |
| `FigSI_eps.R` | Tests sensitivity of the dyad-typology classification to the epsilon parameter used in class assignment. |
| `FigSI_thresh.R` | Tests sensitivity of the dyad-typology classification to the threshold parameter used in class assignment. |

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
