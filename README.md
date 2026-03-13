
# Transboundary Aquifer Analysis Pipeline

This repository contains the code used to construct spatial datasets, process groundwater and irrigation indicators, generate basin typologies, and produce figures used in the associated research project.

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

The following scripts are included in the repository.


### 0_PreporcessIGRAC

| Script | Description |
|------|-------------|
| 0_Peprocess.R | Script in `0_PreporcessIGRAC` stage of pipeline |


### 1_Map

| Script | Description |
|------|-------------|
| Fig1.R | Script in `1_Map` stage of pipeline |


### 2_Wells

| Script | Description |
|------|-------------|
| 1_preferred.R | Script in `2_Wells` stage of pipeline |
| 2_SI_AltCov.R | Script in `2_Wells` stage of pipeline |
| 3_SI_LMRobust.R | Script in `2_Wells` stage of pipeline |
| 4_SI_Umatched.R | Script in `2_Wells` stage of pipeline |
| 5_SI_dropTopSE.R | Script in `2_Wells` stage of pipeline |
| 6_SI_matchArchitecture.R | Script in `2_Wells` stage of pipeline |
| 7_SI_nMin.R | Script in `2_Wells` stage of pipeline |
| 8_SI_loo.R | Script in `2_Wells` stage of pipeline |
| 8_SI_looCountry.R | Script in `2_Wells` stage of pipeline |
| 9_SI_FE.R | Script in `2_Wells` stage of pipeline |
| 9_SI_dclust.R | Script in `2_Wells` stage of pipeline |
| 9_clustRobustSE.R | Script in `2_Wells` stage of pipeline |
| Fig2A.R | Script in `2_Wells` stage of pipeline |
| Fig2B.R | Script in `2_Wells` stage of pipeline |
| FigS1S2.R | Script in `2_Wells` stage of pipeline |
| FigS3.R | Script in `2_Wells` stage of pipeline |
| FigS4.R | Script in `2_Wells` stage of pipeline |
| FigS5.R | Script in `2_Wells` stage of pipeline |
| FigS6.R | Script in `2_Wells` stage of pipeline |
| FigS7.R | Script in `2_Wells` stage of pipeline |
| FirstStages.R | Script in `2_Wells` stage of pipeline |
| TabS1.R | Script in `2_Wells` stage of pipeline |
| TabS2S3.R | Script in `2_Wells` stage of pipeline |
| getWells.R | Script in `2_Wells` stage of pipeline |


### 3_Irrig

| Script | Description |
|------|-------------|
| 0_defMatchPop.R | Script in `3_Irrig` stage of pipeline |
| 1_getBasins.R | Script in `3_Irrig` stage of pipeline |
| 2_buildDataset.R | Script in `3_Irrig` stage of pipeline |
| 2_buildDataset_J.R | Script in `3_Irrig` stage of pipeline |
| 3_GenNonOverlappingCtrl.R | Script in `3_Irrig` stage of pipeline |
| 4_Analyse.R | Script in `3_Irrig` stage of pipeline |
| Analyse_J.R | Script in `3_Irrig` stage of pipeline |
| Fig3.R | Script in `3_Irrig` stage of pipeline |
| FigSIViolin.R | Script in `3_Irrig` stage of pipeline |
| TableBestMatch.R | Script in `3_Irrig` stage of pipeline |
| TableSummary.R | Script in `3_Irrig` stage of pipeline |


### 4_Typology

| Script | Description |
|------|-------------|
| 1_makeDyad.R | Script in `4_Typology` stage of pipeline |
| Fig4_and_mosaic.R | Script in `4_Typology` stage of pipeline |
| FigSI_200.R | Script in `4_Typology` stage of pipeline |
| FigSI_5.R | Script in `4_Typology` stage of pipeline |
| FigSI_IR.R | Script in `4_Typology` stage of pipeline |
| FigSI_eps.R | Script in `4_Typology` stage of pipeline |
| FigSI_thresh.R | Script in `4_Typology` stage of pipeline |


---

# Software Requirements

Main R packages used:

```
sf
lwgeom
dplyr
stringr
future.apply
furrr
progressr
units
```

Install using:

```r
install.packages(c(
"sf","lwgeom","dplyr","stringr",
"future.apply","furrr","progressr","units"
))
```

---

# License

This repository is distributed under the license provided in the `LICENSE` file.
