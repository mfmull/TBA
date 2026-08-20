# Supplementary datasets — codebook

Numbering follows SI section order (wells → irrigation → typology), matching
the order of the `\dataset{}` blocks at the end of `tex/_SI_vFINAL.tex`.
Rebuild with `Rscript assemble_SI_datasets.R` from the repository root.

| | file | unit | rows | cols | produced by |
|---|---|---|---|---|---|
| S1 | `Dataset_S1.csv` | monitored TBA–country segment | 44 | 22 | `2_wells/output/table_si_aquifer_stability.csv` |
| S2 | `Dataset_S2.csv` | IGRAC aquifer–country intersection | 990 | 8 | `0_preprocess_igrac/igrac_segment_inventory.csv` (static) |
| S3 | `Dataset_S3.csv` | TBA–country segment (irrigation) | 929 | 20 | `3_irrigation/output/Dataset_S3.csv` |
| S4 | `Dataset_S4.csv` | dyad | 568 | 23 | `4_typology/output/Dataset_S4.csv` |
| S5 | `Dataset_S5.csv` | dyad | 568 | 29 | `4_typology/output/Data_S3_classification_audit.csv` |

Changes on 2026-08-19:

* **S3 is new** — the segment-level irrigation analysis sample.
* **S4 is the merge** of the former S3 (dyad classification) and the former S5
  (stage-by-stage build). Both were 568 rows on `(code, CC_1, CC_2)`, so the
  join is 1:1 and lossless.
* **S5 is the audit census**, formerly S4. Contents unchanged.

---

## Dataset S3 — irrigation analysis, segment level (929 rows)

One row per transboundary aquifer–country segment. Key: `aq_id` × `CountryName`.
Links to Dataset S2 on `aq_id` (S2's 990 intersections, of which these 929 are
the retained set).

**Inclusion rule.** A column is here only if the 41-series irrigation pipeline
uses it as an outcome, a sample filter, a matching covariate, a regression
covariate, or the grouping factor of the random intercept — plus the
identifiers needed to join and read the file. Controls (HydroSHEDS basins) are
not in this file; the control pool and the 500 non-overlapping realizations
are `3_irrigation/1_data/CtrlNoOverlapHYBAS_B.rds` in the code repository.

| column | pipeline name | role | notes |
|---|---|---|---|
| `aq_id` | `aq_id` | identifier | internal aquifer id; joins to Dataset S2 |
| `IGRAC_Code` | — | identifier | IGRAC aquifer code |
| `AquiferName` | — | identifier | IGRAC aquifer name |
| `CountryName` | `CntrName` | identifier **and** analysis | grouping factor of the country random intercept (143 levels) |
| `CountryCode` | — | identifier | ISO3 |
| `Region` | — | identifier | IGRAC region |
| `lat_c` | `lat_c` | matching covariate | centroid latitude; in `cov_int`, `cov_gini`, `cov_geo` |
| `lon_c` | `lon_c` | matching covariate | centroid longitude; same three sets |
| `RDS_mainDist` | `RDS_mainDist` | matching + regression covariate | km to nearest main road; also `cov_reg` in both double-robust specs |
| `CSI` | `CS_max` | matching covariate **and** outcome | crop-suitability index (0–10,000); in `cov_int`; outcome of `CropSuitInt` |
| `G_CSI` | `giniCSI` | matching + regression covariate | Gini of crop suitability by border distance; in `cov_gini` |
| `NearRiverBorder` | `near_river_border` | regression covariate + sample filter | TRUE = within 50 km of a GSRB river border. `river_cov` in the `\| Rivers` specs; `!near_river_border` defines the river-exclusion robustness sample |
| `Irrig` | `Ir` | outcome + filter + covariate | irrigated fraction of segment. Outcome of `IrIntens`, `IrRivsInt`, `IrIntDoubleRobust`, `IrIntens_exRiv`; defines the `Ir > 0` sample (562 segments); matching covariate of `Overdraft_Irrig` |
| `GWIrrig` | `GW` | sample filter | groundwater-irrigated fraction; defines the `GW > 0` sample (159 segments) |
| `IrrigNeed` | `IrNeed3` | outcome | irrigation need, ≥3 scarce months; outcome of `IrNeedInt` |
| `Overdraft` | `OverIR3` | outcome | outcome of `Overdraft` and `Overdraft_Irrig` |
| `OverIR6` | `OverIR6` | outcome | outcome of `Overdraft6` (≥6 scarce months) |
| `G_Irrig` | `giniIr` | outcome | borderward concentration of irrigation; outcome of `IrGini`, `IrRivsGini`, `IrGiniDoubleRobust`, `IrGini_exRiv` |
| `G_IrrigGW` | `giniGw` | outcome | borderward concentration of GW irrigation; outcome of `GWGini`, `GWRivsGini`, `GWGini_exRiv` |
| `gini_gwpct_gt0` | `gini_gwpct_gt0` | outcome | same under the permissive GW-fraction > 0 attribution; outcome of `GWShareGt0` |

NAs in the Gini columns are structural: a Gini is undefined where the
underlying land-cover band is absent over the segment.

**Dropped from the former `3_irrigation/output/Data_S2.csv`** (six columns, none of
which enters any specification):

| dropped | why |
|---|---|
| `type` | constant `"treat"` — the file is treated-only |
| `IrNeed0` | unused (the analysis uses the ≥3-month variant, `IrrigNeed`) |
| `gini_crpGWSpct_gt0` | unused **and** an exact duplicate of `G_IrrigNeed` (max abs difference 0 over all 929 rows) |
| `G_IrrigNeed` | descriptive only; no specification uses it |
| `area_km2` | descriptive only; the 0.25×–4× control-size caliper is applied upstream, in the control construction |
| `Countries` | riparian list; already in Dataset S2 on the same `aq_id` |

Dropping `area_km2` and `G_IrrigNeed` removes their two rows from the SI
summary-statistics table; `3_irrigation/47_tables.R` is updated accordingly and now
also reports `OverIR6` and `gini_gwpct_gt0`, which are outcomes.

---

## Dataset S4 — dyad classification and build, dyad level (568 rows)

One row per dyad. Key: `code` × `CC_1` × `CC_2`. Merge of the former Dataset S3
(classification + mosaic characteristics) and the former Dataset S5 (the
stage-by-stage build ledger).

| column | source | description |
|---|---|---|
| `code`, `name` | both | IGRAC aquifer code and common name |
| `CC_1`, `CC_2` | both | ISO3 codes of the two countries |
| `class` | old S3 | **finalised** class under current abstraction (IU0, IU1, UB, DA, BL) |
| `classF` | old S3 | finalised class under the future irrigation scenario |
| `cooperation_level` | old S3 | SDG 6.5.2 level, 0–3 (Methods) |
| `has_river_border` | old S3 | shared border is river-defined |
| `dist_class` | old S3 | distance class of the aquifer footprint |
| `class_stage1`, `classF_stage1` | old S5 (`class_1`, `classF_1`) | Stage 1: published segment-based classification |
| `class_stage2`, `classF_stage2` | old S5 (`class_2`, `classF_2`) | Stage 2: after the dyad-border correction |
| `measured_stage2` | old S5 | whether Stage 2 measured this dyad against its own border |
| `changed_stage2` | old S5 (`changed_2`) | Stage 2 changed the class |
| `changed_stage3` | old S5 (`changed_3`) | Stage 3 (audit adjudication) changed the class — 41 dyads |
| `changed_overall` | old S5 | class differs between Stage 1 and the finalised result — 78 dyads |
| `near_share_stage1`, `near_share_stage2` | old S5 (`near_share_1/2`) | near-zone share before and after the geometric correction |
| `audit_error_class`, `audit_tier` | old S5 | audit error class and evidence tier where one applies |
| `comment_stage2`, `comment_stage3` | old S5 | written note explaining each reclassification |

`class_3` and `classF_3` of the old Dataset S5 are **not** carried: Stage 3 *is*
the finalised classification, so they were identical to `class` and `classF`
row for row. `4_typology/55_summaries.R` asserts that identity on every run
rather than assuming it.

Class counts, finalised: IU0 87, IU1 111, UB 205, DA 125, BL 40.
Future scenario: IU0 74, IU1 91, UB 105, DA 211, BL 87.

---

## Dataset S5 — external audit census, dyad level (568 rows)

Unchanged from the former Dataset S4. One row per dyad: published class, audit
stratum, claim tested, verdict, proposed class, evidence tier, confidence,
Tier-1 re-citation outcome, and the geometric validation of every candidate
reclassification.
