# =============================================================================
# 51_config.R -- CONFIGURATION (transboundary aquifer dyad typology)
#
# Every substantive analytical choice lives here and nowhere else. The defaults
# ARE the FINALISED specification -- stage 3 of the four-stage build described
# below, i.e. dyad-specific border zones plus the adjudicated audit corrections.
# classify_dyads() with no overrides must reproduce reference_class /
# reference_classF, and the stage1_published variant must reproduce the
# published segment-based benchmark. 53_run_main.R asserts both.
#
# The pipeline is self-contained on 1_data/: FinalDiads.csv is the analysis
# input, and 1_data/geeOut/ holds the Google Earth Engine exports it was built
# from. Sourcing this file defines CFG and nothing else; it fits nothing and
# writes nothing.
# =============================================================================

.this_file <- (function() {
  for (fr in rev(sys.frames())) {
    f <- tryCatch(get("ofile", envir = fr), error = function(e) NULL)
    if (is.character(f) && grepl("51_config[.]R$", f)) return(normalizePath(f))
  }
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) && grepl("51_config[.]R$", f[1])) return(normalizePath(f[1]))
  if (file.exists("51_config.R")) return(normalizePath("51_config.R"))
  stop("Cannot locate 51_config.R; run from the 4_typology folder.")
})()
PROJECT_ROOT <- dirname(.this_file)


# =============================================================================
# THE FOUR-STAGE BUILD
#
#   1  Segment-based classification, as published, plus the distance and
#      parameter sensitivities. All of CFG$variants below except the three
#      stage_* entries are evaluated on THIS frame.
#   2  Dyad-border correction: zones re-measured against the dyad's own border
#      rather than against every international land border. CFG$dyad_zones.
#   3  External audit adjudication, re-derived on the stage-2 frame.
#      CFG$manual_class = AUDIT_FINAL.
#   4  Inference -- Fig. 4B/C, the mosaics, the honest test -- on the output of
#      stage 3, which is the FINALISED dataset and the preferred specification.
#
# The order is 1 -> 2 -> re-derive from/fromF -> 3 -> 4, and it is not
# interchangeable: the audit validated its findings against the DYAD'S OWN
# border, so its corrections are expressed in the post-stage-2 frame. Applying
# them to segment-based data would justify a correction in a frame the data is
# not in.
#
# =============================================================================
# EXTERNAL AUDIT RECLASSIFICATIONS (SI Section S5)
#
# Each row is a dyad whose classification is contradicted by external evidence,
# with the class the rule currently produces (`from`) and the class the rule
# would produce from corrected inputs (`to`). `tier` records the strength of the
# evidence: 1 = peer-reviewed or government / national / international agency
# report; 2 = any other source. Tier 1 alone is retained as the stage3_tier1
# variant, so a reader can see exactly how much rests on weaker evidence: three
# dyads on `class`.
#
# These tables ARE applied in the preferred specification (author decision,
# 12 Aug 2026). The finalised dataset is stage 3, and the Results, Fig. 4 and
# every inference in 53-57 are computed on it. The published segment-based
# classification remains reproducible and asserted -- see
# reference_class_published and the stage1_published variant.
#
# Coverage: a CENSUS. All 568 dyads on all 409 aquifers were audited against
# external sources, producing 134 candidate reclassifications. Every candidate was
# then checked against the aquifer polygons themselves
# (0_PreprocessIGRAC/TBA_Split_2025.shp): the feature the finding rests on must
# fall inside the polygon of the side the claim declares empty AND inside the
# 100 km band of that dyad's own border. The auditors had no GIS and judged this
# from names and descriptions; often they were wrong.
#
# Each of the 134 falls in exactly one ERROR CLASS, and each class was adjudicated
# by the authors as a class rather than dyad by dyad (SI Section S5.3):
#
#   class                                           n   ->DA/BL  classF  status
#   A settlement below the urban-cluster floor     31     0        22    EXCLUDED
#   B settlement above the floor, yet zero         23     0        20    APPLIED
#   C GW irrigation absent from the whole polygon  14    13         4    APPLIED
#   D irrigation on the polygon, outside the band   2     2         1    APPLIED
#   E interior cropland understated, BL share       3     3         1    APPLIED
#   F undecidable even after the buffer             2     2         -    EXCLUDED
#   G rejected: feature outside the polygon        56     4         -    EXCLUDED
#   H rejected: inside polygon, beyond the band     2     0         -    EXCLUDED
#   I rejected: evidence on the side not in question 1    1         -    EXCLUDED
#
# COORDINATE-UNCERTAINTY BUFFER (author decision, 12 Aug 2026). A strict
# point-in-polygon test is not meaningful when BOTH the feature coordinate and the
# IGRAC polygon boundary carry locational error. A finding is therefore accepted
# when its distance outside the polygon is within
#     (the feature's own stated coordinate uncertainty) + 2 km
# where the 2 km is the polygon delineation tolerance. Applied uniformly to every
# rejected and undecidable finding, not only to the marginal ones. Chosen because
# 2 km is the smallest tolerance at which the undecidable bucket resolves (26 of
# 28) while disturbing only 2 of the 58 polygon-rejected findings; at 3 km the
# latter jumps to 7 with no further gain. Without the polygon term the rule is
# perverse -- a town geocoded to +/-1 km at 2 km out fails while one geocoded to
# +/-5 km at 3 km out passes, penalising better geocoding.
# The buffer accepted 28 findings (26 from F, 2 from G); 10 of them are
# sub-threshold settlements and so fall in class A and stay excluded, leaving 18
# new corrections. It is what finally admits S018 Caplina-Concordia, the dyad that
# prompted the audit (0.4 km outside the Chilean polygon, +/-5 km coordinate).
#
# A is excluded as DEFINITIONAL, not mis-measured: UR is GHS-SMOD 2005 thresholded
# at .gt(20), i.e. urban cluster/centre only (>=~5,000 people at >=300/km2), so a
# village is a RURAL class and scores zero by construction. Its 21 corrections
# move dyads only within IU0/IU1/UB and never touch DA or BL. The caveat lives in
# the methods (SI S3, dyad typology). Applying it as well would give
# 71/114/213/127/43 for class.
#
# D and E ARE APPLIED (author decision, 12 Aug 2026). They were the two classes
# that turn on quantities the dyad-border correction recomputes -- the 100 km band
# for D, the near-zone share for E -- and the re-check they required has now been
# done. Outcome: of the five rows, one is withdrawn as corroborated (AS002
# EGY-PSE, class E: stage 2 reaches DA on geometry alone), one is weakened
# (AS130 IRQ-SAU, class D: the baseline moves UB -> IU1 and the audited
# irrigation now lands UB rather than DA), and the remaining three stand.
# Classes B and C, by contrast, are provably invariant to the re-banding: both
# were selected on a WHOLE-POLYGON zero (max 0.84 ha of groundwater-irrigated
# cropland across all corrected sides), and re-banding redistributes area within
# a polygon without creating area the polygon does not contain. Their `from`
# values still move where stage 2 changes the baseline, which is what the
# manual_class guard is for: it stops the run rather than applying a stale
# correction, and it is the reason the tables below carry re-derived values.
#
# F now holds only the two findings that remain undecidable after the buffer, both
# far outside the polygon on a diffuse feature: N082 MEX-USA (57 km) and
# EU086/EU088-EU089 BLR-UKR (43 km), each with a 25 km locational uncertainty.
#
# G, H and I are rejected on the geometry and need no further work; G's features
# lie a median 26 km outside the polygon and up to 144 km.
#
# Can stage 2 ADMIT a correction the segment frame rejected? No. H's two rows
# were rejected because the feature sits inside the polygon but beyond the band
# of the dyad's own border -- 163.6 km for AF005 NAM-ZAF, 116.9 km for AF059
# DJI-ERI, both already measured in the stage-2 frame -- and I's single row
# rests on evidence for the side not in question. F's two survivors are
# bilateral aquifers and so are not measured by the correction at all.
#
# The 42 adjudicated corrections (B, C, D and E) are applied to BOTH class and
# classF, and 41 survive re-derivation. classF propagation is monotone by
# construction: classF takes pmax(current, future) per zone, so a correction
# establishing presence for `class` also establishes it for `classF`. A number
# of rows are guarded no-ops on classF (toF == fromF) because the water-stressed
# cropland proxy already supplied presence in that zone; they are retained so
# the guard still checks them.
#
# WITHDRAWN and not to be reinstated: that the urban layer is unclipped (it is
# clipped; max urban/segment-area = 0.97) and that many segments carry no
# land-use values (all 201 all-zero segments are present with literal zeros; 128
# have land use over the full polygon but none in the 100 km band, which is the
# intended signal). Both rested on a 100x unit error.
#
# UNIT NOTE -- RESOLVED EMPIRICALLY, AND THE TWO LAYERS DIFFER BY TEN.
# getTypologyOG.js divides every band by 1e7 and its header calls the result
# thousands of hectares. That is right for the URBAN band, which is a binary
# mask, and wrong by a factor of ten for the CROPLAND bands, because the
# Crop1k_0307 raster is not a 0-1 fraction. Established from the data, not the
# code: (i) `Area` is raw pixel area in m2 -- exported country areas reproduce
# the real ones to 1-4% (Aruba 182 vs 180 km2); (ii) at 1e8 m2/unit the urban
# band would give Bangladesh 5.9x its own area, while the cropland bands reach
# at most 0.43 of theirs; (iii) at 1e8 m2/unit the cropland bands reproduce
# national statistics (India 37.1 Mha of GW-irrigated cropland against ~39
# reported), and at 1e7 every national total is exactly 10x too low.
#
#   one unit of CR / IR / GW / CR3 / GW3 / IR3 = 1e8 m2 = 10,000 ha
#   one unit of UR                             = 1e7 m2 =  1,000 ha
#
# So eps = 1e-3 is 10 ha of groundwater-irrigated cropland AND 1 ha of urban
# land; eps_lo = 1e-5 is 0.1 ha and 0.01 ha; bl_min_ha = 1e-2 is 100 ha, and
# applies to cropland only. AUTHOR DECISION (12 Aug 2026): keep these settings.
# NOTHING in the classification depends on the conversion -- the rule is
# evaluated in raw units throughout -- but every hectare figure quoted in the
# manuscript and the SI does, and the urban limb of the interior test is ten
# times more sensitive than the cropland limb.
#
# CR3 MASKING -- OPEN, DOCUMENTED, NOT RESOLVED HERE. In getTypologyOG.js, CR3
# is masked on GWS.gt(2) while GW3 and IR3 use BWS.gt(2), although all three are
# named *_GWS3 and the script header calls all three "groundwater stress". The
# inline comments are the accurate ones. This matters beyond Fig. 4: future_var
# is CR3 and classF is pmax(current, CR3) per zone, and FBW -- the Overdraft
# axis -- is (GW3_1 + GW3_2)/(GW_1 + GW_2), so an unintended mask would move both
# the future classification and the overdraft metric. AUTHOR DECISION (12 Aug
# 2026): document it as a caveat (SI S3) rather than re-export. Every classF and
# FBW number in this pipeline is computed against CR3 as currently masked.
# =============================================================================
# The tables below are the corrections AS RE-DERIVED ON THE STAGE-2 FRAME
# (see SI Section S5.8). `from` and `fromF` are the classes the rule produces
# after the dyad-border correction, not before it, so the manual_class guard
# checks the corrections against the data they are actually applied to. Where
# stage 2 already reached the audit's target on both classifications the row is
# withdrawn as corroborated rather than applied -- one row, AS002 EGY-PSE.
# AUDIT_FINAL therefore carries 41 rows, not 42, and AUDIT_FINAL_T1 carries 35.
#
# Six rows changed under re-derivation, all documented in
# output/table_audit_rederivation.csv:
#   AF009 BWA-ZWE, AF009 ZAF-ZWE, EB067-EB069 BGR-GRC and BGR-MKD -- stage 2
#     dissolved the BL future class these rows asserted; toF is re-inferred DA.
#   AF056 DZA-NER -- stage 2 moves the baseline IU1 -> IU0, and re-running the
#     rule with the audited urban presence lands IU1, not UB.
#   AS130 IRQ-SAU (class D) -- stage 2 moves the baseline UB -> IU1, and the
#     audited irrigation now lands UB, not DA. This is the class-D re-check the
#     build order required, and it weakens the correction rather than voiding it.
#
# `to` is re-derived from the RULE for error classes B and D, whose assertion
# names a single quantity on a named side (urban land, groundwater-irrigated
# cropland, in the far zone): raise it to the detection threshold and let
# classify_dyads() decide. That method reproduces the author's recorded target
# on 23 of the 25 class-B/D rows on the stage-1 frame, which is what licenses
# it on the stage-2 frame. Classes C and E assert presence without locating it,
# so their targets remain the author's adjudication, carried forward unchanged.
#
# classF is INFERRED throughout, never asserted:
#     toF = to                  when fromF is BL and to is DA
#         = max(fromF, to)      in the IU0 < IU1 < UB < DA < BL order otherwise
# The exception is not ad hoc -- classes C and E add land use to the interior,
# which enlarges the denominator of the near-zone share, so a BL classification
# resting on that share collapses in the future frame as it does in the current
# one. The rule reproduces all 42 of the author's recorded toF values.

AUDIT_FINAL_T1 <- data.frame(
  code  = c(
    "AF002", "AF003", "AF005", "AF005", "AF009", "AF009", "AF019", "AF042",
    "AF056", "AF059", "C009", "EB023", "EB052", "EB067 - EB069", "EB077",
    "EB078", "EB090", "S004", "AF001", "AF009", "AF052", "AF052", "AF052",
    "AF054", "AF058", "AF067", "AS128", "AS161", "N002", "N013", "N025",
    "S018", "AS130", "AF074; AF078", "AS089"),
  CC_1  = c(
    "NAM", "MOZ", "BWA", "BWA", "BWA", "ZAF", "MWI", "CMR", "DZA", "DJI",
    "BLZ", "BIH", "MNE", "BGR", "BGR", "BGR", "MDA", "BRA", "LSO", "BWA",
    "CMR", "CMR", "SDN", "GHA", "GMB", "DZA", "IRQ", "ISR", "CAN", "MEX",
    "MEX", "CHL", "IRQ", "DZA", "KHM"),
  CC_2  = c(
    "ZAF", "ZAF", "NAM", "ZAF", "ZWE", "ZWE", "ZMB", "NGA", "NER", "ETH",
    "GTM", "HRV", "XKX", "GRC", "GRC", "GRC", "ROU", "GUY", "ZAF", "ZAF",
    "NGA", "TCD", "TCD", "TGO", "SEN", "LBY", "KWT", "JOR", "USA", "USA",
    "USA", "PER", "SAU", "MAR", "VNM"),
  from  = c(
    "IU0", "IU1", "IU0", "IU0", "IU1", "IU1", "IU1", "IU1", "IU0", "IU1",
    "IU0", "IU1", "IU0", "IU1", "IU0", "IU0", "IU1", "IU0", "UB", "UB",
    "UB", "UB", "UB", "UB", "UB", "UB", "UB", "IU1", "UB", "UB", "IU0",
    "IU1", "IU1", "BL", "BL"),
  to    = c(
    "UB", "UB", "IU1", "IU1", "UB", "UB", "UB", "UB", "IU1", "UB", "UB",
    "UB", "IU1", "UB", "IU1", "IU1", "UB", "IU1", "DA", "BL", "DA", "DA",
    "DA", "DA", "DA", "DA", "DA", "BL", "DA", "DA", "IU1", "BL", "UB", "DA",
    "DA"),
  fromF = c(
    "IU0", "IU1", "IU0", "IU0", "DA", "DA", "IU1", "IU1", "IU0", "IU1",
    "IU0", "IU1", "IU0", "DA", "BL", "IU0", "BL", "IU0", "BL", "BL", "DA",
    "DA", "DA", "DA", "DA", "UB", "DA", "BL", "DA", "DA", "IU0", "IU1",
    "IU1", "DA", "DA"),
  toF   = c(
    "UB", "UB", "IU1", "IU1", "DA", "DA", "UB", "UB", "IU1", "UB", "UB",
    "UB", "IU1", "DA", "BL", "IU1", "BL", "IU1", "DA", "BL", "DA", "DA",
    "DA", "DA", "DA", "DA", "DA", "BL", "DA", "DA", "IU1", "BL", "UB", "DA",
    "DA"),
  class = c(
    "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B",
    "B", "B", "B", "B", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C",
    "C", "C", "C", "C", "D", "E", "E"),
  tier  = 1L,
  stringsAsFactors = FALSE)

AUDIT_FINAL <- data.frame(
  code  = c(
    "AF002", "AF003", "AF005", "AF005", "AF009", "AF009", "AF019", "AF042",
    "AF056", "AF059", "C009", "EB023", "EB052", "EB067 - EB069", "EB077",
    "EB078", "EB090", "S004", "AF001", "AF009", "AF052", "AF052", "AF052",
    "AF054", "AF058", "AF067", "AS128", "AS161", "N002", "N013", "N025",
    "S018", "AS130", "AF074; AF078", "AS089", "AF077", "AS008", "CB004",
    "EB051", "EB067 - EB069", "AF063"),
  CC_1  = c(
    "NAM", "MOZ", "BWA", "BWA", "BWA", "ZAF", "MWI", "CMR", "DZA", "DJI",
    "BLZ", "BIH", "MNE", "BGR", "BGR", "BGR", "MDA", "BRA", "LSO", "BWA",
    "CMR", "CMR", "SDN", "GHA", "GMB", "DZA", "IRQ", "ISR", "CAN", "MEX",
    "MEX", "CHL", "IRQ", "DZA", "KHM", "DZA", "GEO", "DOM", "MNE", "BGR",
    "EGY"),
  CC_2  = c(
    "ZAF", "ZAF", "NAM", "ZAF", "ZWE", "ZWE", "ZMB", "NGA", "NER", "ETH",
    "GTM", "HRV", "XKX", "GRC", "GRC", "GRC", "ROU", "GUY", "ZAF", "ZAF",
    "NGA", "TCD", "TCD", "TGO", "SEN", "LBY", "KWT", "JOR", "USA", "USA",
    "USA", "PER", "SAU", "MAR", "VNM", "MAR", "RUS", "HTI", "SRB", "MKD",
    "LBY"),
  from  = c(
    "IU0", "IU1", "IU0", "IU0", "IU1", "IU1", "IU1", "IU1", "IU0", "IU1",
    "IU0", "IU1", "IU0", "IU1", "IU0", "IU0", "IU1", "IU0", "UB", "UB",
    "UB", "UB", "UB", "UB", "UB", "UB", "UB", "IU1", "UB", "UB", "IU0",
    "IU1", "IU1", "BL", "BL", "IU1", "IU1", "IU1", "IU0", "IU1", "UB"),
  to    = c(
    "UB", "UB", "IU1", "IU1", "UB", "UB", "UB", "UB", "IU1", "UB", "UB",
    "UB", "IU1", "UB", "IU1", "IU1", "UB", "IU1", "DA", "BL", "DA", "DA",
    "DA", "DA", "DA", "DA", "DA", "BL", "DA", "DA", "IU1", "BL", "UB", "DA",
    "DA", "UB", "UB", "UB", "IU1", "UB", "DA"),
  fromF = c(
    "IU0", "IU1", "IU0", "IU0", "DA", "DA", "IU1", "IU1", "IU0", "IU1",
    "IU0", "IU1", "IU0", "DA", "BL", "IU0", "BL", "IU0", "BL", "BL", "DA",
    "DA", "DA", "DA", "DA", "UB", "DA", "BL", "DA", "DA", "IU0", "IU1",
    "IU1", "DA", "DA", "IU1", "IU1", "IU1", "IU1", "DA", "UB"),
  toF   = c(
    "UB", "UB", "IU1", "IU1", "DA", "DA", "UB", "UB", "IU1", "UB", "UB",
    "UB", "IU1", "DA", "BL", "IU1", "BL", "IU1", "DA", "BL", "DA", "DA",
    "DA", "DA", "DA", "DA", "DA", "BL", "DA", "DA", "IU1", "BL", "UB", "DA",
    "DA", "UB", "UB", "UB", "IU1", "DA", "DA"),
  class = c(
    "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B",
    "B", "B", "B", "B", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C",
    "C", "C", "C", "C", "D", "E", "E", "B", "B", "B", "B", "B", "D"),
  tier  = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2),
  stringsAsFactors = FALSE)

CFG <- list(

  ## ---- paths (all relative to PROJECT_ROOT) --------------------------------
  root        = PROJECT_ROOT,
  data_dir    = file.path(PROJECT_ROOT, "1_data"),
  dyad_file   = file.path(PROJECT_ROOT, "1_data", "FinalDiads.csv"),
  gee_dir     = file.path(PROJECT_ROOT, "1_data", "geeOut"),
  # 2026-08-19. The stage-by-stage ledger produced by the four-stage audit
  # build (archived in derived/audit/typology_derived_tables.tar.gz as
  # Data_S4_dyad_stages.csv). It is a FIXED INPUT here, not a pipeline
  # product: stages 1-3 are a completed historical record of how the
  # published classification was corrected, and re-running 53_run_main.R
  # reproduces the stage-3 result, not the path taken to it. 55_summaries.R
  # joins it onto the classification so that the SI ships ONE dyad-level
  # dataset (Dataset S4) instead of the two it used to (old S3 + old S5),
  # and asserts stage 3 == the classification the pipeline just computed.
  stages_file = file.path(PROJECT_ROOT, "1_data", "dyad_stages.csv"),
  out_dir     = file.path(PROJECT_ROOT, "output"),
  fig_dir     = file.path(PROJECT_ROOT, "figure"),
  cache_dir   = file.path(PROJECT_ROOT, "derived", "cache"),
  audit_dir   = file.path(PROJECT_ROOT, "derived", "audit"),
  log_file    = file.path(PROJECT_ROOT, "derived", "audit", "run_log.txt"),

  ## ---- STAGE 2: dyad-specific border zones ---------------------------------
  # Every zone column in FinalDiads.csv measures distance to `lb`, the global
  # land-border layer, over an aquifer-country SEGMENT. A segment facing three
  # riparians therefore carries the same near-zone area in all three dyads.
  # gee/dyad_border_corrections.js re-measures each zone against the stretch of
  # THIS dyad's border lying under the aquifer, one row per (code, CC, CC_other).
  #
  # TRUE in the preferred specification: the corrected geometry is what the
  # analysis is run on (SI Section S5.8). Set FALSE to reproduce the published
  # segment-based classification, which the stage1_published variant does.
  dyad_zones     = TRUE,
  dyad_zone_dir  = file.path(PROJECT_ROOT, "1_data", "geeOut"),
  # GLOBS, not filenames: the export is partitioned by aquifer and
  # gee/finish_dyad_10k.js re-runs the ones that time out, so a zone arrives in
  # many parts. read_dyad_zone() concatenates every match. Names are the zone
  # suffixes in FinalDiads.csv, so _B1 is the 10 km near zone and _B2 the 100 km
  # far zone. dyadB5k_* is deliberately absent: that export is incomplete and
  # the 5 km zone (_B11) is a stage-1 robustness variant, not a corrected one.
  dyad_zone_files = list("_B1" = "dyadB10k_*.csv",
                         "_B2" = "dyadB100k_*.csv"),
  dyad_zone_vars  = c("Area", "GW", "GW3", "IR", "IR3", "UR", "CR3"),
  # Relative tolerance when a key appears in more than one part. Earth Engine
  # re-serialises the reduction on a re-export, so repeated keys agree to the
  # last one or two ulps and NOT byte-identically; anything larger is a genuine
  # disagreement and stops the run.
  dyad_zone_dup_tol = 1e-6,
  # A corrected value may legitimately exceed the segment value it replaces:
  # `lb` exists only where two countries share a LAND boundary, while the dyad
  # zone is the counterpart segment buffered and can cross water. 105 values on
  # 40 dyads do so (median +1.7%), all on genuinely adjacent dyads (gap_m = 0).
  # The run stops only above this relative excess, which catches a join error
  # without tripping on the real effect. See table_dyad_border_excess.csv.
  dyad_zone_tol         = 1e-6,
  dyad_zone_excess_stop = 5,
  dyad_zone_strict      = FALSE,

  ## ---- optional dyad rebuild (52_core.R::build_dyads) ----------------------
  # FinalDiads.csv ships as the frozen analysis input. Rebuilding it from
  # geometry needs three files that are NOT part of this folder (and are not
  # currently anywhere in the repository): the IGRAC aquifer split, the GSRB
  # river-border layer, and the Fraser agreement table. Point these at real
  # files to enable the rebuild; build_dyads() stops with a clear message
  # naming whichever is missing. Nothing else in the pipeline needs them.
  tba_shapefile   = file.path(PROJECT_ROOT, "..", "0_PreprocessIGRAC",
                              "TBA_Split_2025.shp"),
  rivers_shapefile = NULL,   # e.g. ".../GSRB/GSRB_Level0.shp"
  agreements_file  = NULL,   # e.g. ".../Fraser2023_Agreements/_agreements.csv"
  river_border_tol_m = 5000, # buffer around the shared boundary, metres
  # Burchi (2018) https://doi.org/10.1016/j.ejrh.2018.04.007 -- aquifers whose
  # riparians are parties to a formal arrangement, overriding the Fraser score.
  # NOTE (audit 4.5): this override is applied aquifer-wide. For NSAS, NWSAS,
  # Guarani and the bilateral cases every riparian pair is a party, but for
  # Iullemmeden and Taoudeni coverage depends on the MoU vintage. It feeds only
  # the cooperation mosaic.
  burchi_codes = c("AF056", "AS126", "AF069", "AF063", "S021",
                   "EU024", "AF064", "N015", "N001"),
  burchi_score = 3,

  ## ---- classification ------------------------------------------------------
  # The typology is a deterministic, ordered rule. Every variant in the
  # robustness registry is one of these five settings changed; no variant
  # re-implements the rule.
  #
  #   IU0  neither side has land use within the far zone      (interior, both)
  #   IU1  one side has neither GW irrigation nor urban land  (interior, one)
  #   UB   at least one side lacks GW irrigation in the far zone, but has
  #        urban land -- bilateral agricultural interaction precluded
  #   BL   >thresh of a side's total is inside the NEAR zone  (border-localised)
  #   DA   otherwise: bilateral, distributed across the aquifer
  #
  # NOTE (audit 2.2): UB as coded does not require urban land on BOTH sides;
  # the Methods wording is being aligned to the code, not the other way round.
  land_var   = "GW",      # GW = groundwater-irrigated cropland; IR = all irrigated
  near_zone  = "_B1",     # column suffix for the near buffer  (_B1 = 10 km, _B11 = 5 km)
  far_zone   = "_B2",     # column suffix for the far buffer   (_B2 = 100 km, _B22 = 200 km)
  eps        = 1e-3,      # detection threshold = 1 ha. One unit is 1e7 m2 =
                          # 1,000 ha (getTypologyOG.js divides pixel area by
                          # 1e7), so 1e-3 units = 1 ha. Earlier comments said
                          # 10 ha; that was a 10x unit misreading, not a change
                          # of specification. See the unit note above.
  thresh     = 0.99,      # concentration threshold for BL
  # Absolute floor on the side that triggers BL, in the same units as eps
  # (1 unit = 1,000 ha), so bl_min_ha = 1e-2 is 10 ha.
  # ZERO in the preferred specification -- this changes no published number.
  # BL is defined as a share, so it can fire on a side holding a trivial area:
  # 24 of the 45 BL dyads bind on under 10 ha and 16 on under 5 ha, and in
  # several the trivial side sits on a dyad of thousands of hectares (Sarakhas
  # binds on 14.7 ha out of 13,511; Umm er Radhuma ARE-SAU on 27.9 out of
  # 21,345). The bl_floor variant below tests whether BL survives a floor.
  bl_min_ha  = 0,
  # STAGE 3: the adjudicated audit corrections, Tier 1 and Tier 2, re-derived
  # against the stage-2 frame. APPLIED in the preferred specification -- the
  # finalised dataset is stages 1 -> 2 -> 3, and that is what the Results and
  # the inference in 53-57 are run on. Set NULL (with dyad_zones = FALSE) to
  # reproduce the published classification; stage1_published does exactly that.
  manual_class = AUDIT_FINAL,
  # Future scenario: per zone, land use is taken as the maximum of current
  # extent and cropland under >=3-month green-water scarcity, holding urban
  # extent constant. pmax avoids double counting and is a lower bound on the
  # union, so the scenario is conservative (audit 4.6).
  future_var = "CR3",

  ## ---- class levels and colours -------------------------------------------
  class_levels = c("IU0", "IU1", "UB", "DA", "BL"),
  # "published" reproduces Fig. 4 as printed. "cvd" is an accessible
  # alternative -- see the note below; switching is a one-word change and
  # affects nothing but the rendering.
  palette = "published",
  # The published palette. NOTE: checked with a colour-vision validator, the
  # adjacent pair IU0 #ABABAB and IU1 #A6CEE3 separates by only dE 10.2 for
  # NORMAL vision (a floor of 15 is the usual target) and 9.4 under deuteranopia
  # -- these two categories are hard to tell apart for everyone, not only for
  # colourblind readers, and they sit side by side in the alluvial and the
  # legend. Kept as the default so the figure reproduces as published.
  class_cols_published = c(IU0 = "#ABABAB", IU1 = "#A6CEE3", UB = "#6A3D9A",
                           DA  = "#FF7F00", BL  = "#E31A1C"),
  # Validated alternative: passes adjacent-pair separation for normal vision
  # (worst dE 16.6) and for deuteranopia, protanopia and tritanopia (worst
  # dE 13.2). IU0 stays neutral grey deliberately -- it is the "no land use on
  # either side" category, and grey carries that meaning; the greens and reds
  # are Paul Tol / Okabe-Ito steps. Every mark is also directly labelled, so
  # identity never rests on colour alone.
  class_cols_cvd = c(IU0 = "#A0A0A0", IU1 = "#117733", UB = "#6A3D9A",
                     DA  = "#E08214", BL  = "#CC3311"),

  ## ---- Fig. 4B/C statistics ------------------------------------------------
  # FIX (audit 2.3): the code used conf.level = 0.90 / threshold = 0.10 while
  # the Fig. 4 caption said 95%. The Tukey letter groupings are identical at
  # both levels, so the figure does not change -- but code and caption now
  # agree. Set these back to 0.90 / 0.10 to reproduce the legacy run exactly.
  tukey_conf   = 0.95,
  tukey_alpha  = 0.05,
  # FIX (audit 3.2): standard errors divided by sqrt(n()), which counts rows
  # whose value is NA. TRUE uses sqrt(sum(!is.na(x))).
  se_drop_na   = TRUE,
  # FIX (audit 4.3): 11 of 20 cells in the class x cooperation table have
  # expected counts below 5, where the chi-square approximation is unreliable.
  chisq_simulate = TRUE,
  chisq_B        = 10000L,
  # Reported alongside the ANOVA as a distribution-free check (audit 4.1):
  # FGW and FBW are shares in [0,1] with group SDs spanning two orders of
  # magnitude, which strains the homoscedasticity assumption behind Tukey.
  report_kruskal = TRUE,
  # Dyads are not independent -- a country appears in up to 39 of them, an
  # aquifer in up to 10 (audit 4.2). Fig. 4B/C is framed as descriptive; this
  # adds a country-clustered bootstrap of the group means as a check.
  cluster_boot   = TRUE,
  cluster_boot_reps = 2000L,
  cluster_boot_seed = 20260808L,

  ## ---- distance class (mosaic) ---------------------------------------------
  # Computed from aquifer EXTENT within each buffer, not from observed use:
  # the question is whether inland exploitation is geometrically feasible
  # (audit 2.4b).
  dist_class_levels = c("0-10km", "10-100km", "100+km"),

  ## ---- robustness registry -------------------------------------------------
  # Each entry is the preferred specification with ONE setting changed. `si`
  # names the SI figure the variant corresponds to in v3.
  #
  # TWO KINDS OF ENTRY, and the distinction matters.
  #
  #   STAGE variants reproduce the intermediate datasets of the build:
  #   stage1_published is the segment-based classification exactly as published,
  #   stage2_dyad_border adds the dyad-border correction, stage3_tier1 restricts
  #   the audit to Tier 1 evidence. The preferred specification IS stage 3 with
  #   both tiers, so these show what each correction did.
  #
  #   PARAMETER variants -- buffer widths, detection threshold, land variable,
  #   concentration threshold, BL floor -- are all evaluated ON THE STAGE-1
  #   FRAME, with dyad_zones = FALSE and manual_class = NULL. That is deliberate
  #   (author decision, 12 Aug 2026). These sensitivities are properties of the
  #   segment-based measurement the published analysis rests on; running them on
  #   the corrected frame would confound a measurement sensitivity with a
  #   correction, and the reader could not tell which moved a number.
  variants = list(
    stage1_published = list(
      label = "Stage 1: segment-based zones, no audit corrections (published)",
      si = "S18", dyad_zones = FALSE, manual_class = NULL),
    stage2_dyad_border = list(
      label = "Stage 2: dyad-specific border zones, no audit corrections",
      si = "S19", dyad_zones = TRUE, manual_class = NULL),
    stage3_tier1 = list(
      label = "Stage 3: dyad borders + Tier 1 audit corrections only",
      si = "S20", dyad_zones = TRUE, manual_class = AUDIT_FINAL_T1),
    near5 = list(label = "Near zone 5 km instead of 10 km",
                 si = "S12", near_zone = "_B11",
                 dyad_zones = FALSE, manual_class = NULL),
    far200 = list(label = "Far zone 200 km instead of 100 km",
                  si = "S13", far_zone = "_B22",
                  dyad_zones = FALSE, manual_class = NULL),
    eps_lo = list(label = "Detection threshold 0.01 ha instead of 1 ha",
                  si = "S10", eps = 1e-5,
                  dyad_zones = FALSE, manual_class = NULL),
    irrig = list(label = "All irrigated cropland instead of groundwater-irrigated",
                 si = "S11", land_var = "IR",
                 dyad_zones = FALSE, manual_class = NULL),
    # Not cited in v3. Retained as a canonical output: BL grows but keeps its
    # high-overdraft profile, which is the point worth reporting (audit 4.4).
    thresh50 = list(label = "Concentration threshold 50% instead of 99%",
                    si = "S14", thresh = 0.50,
                    dyad_zones = FALSE, manual_class = NULL),
    # BL is a share and can be triggered by a side holding a trivial area.
    # This variant requires the concentrated side to hold at least 10 ha,
    # ten times the 1 ha detection threshold. It removes precisely the 24 BL
    # dyads that bind on a side under 10 ha, which is why BL lands on 21.
    # 10 ha is a round number and
    # not a natural constant; the SI reports the 50 ha and 500 ha results
    # alongside, and BL keeps the highest overdraft exposure at all three.
    bl_floor = list(label = "BL requires the concentrated side to exceed 10 ha",
                    si = "S15", bl_min_ha = 1e-2,
                    dyad_zones = FALSE, manual_class = NULL)),

  ## ---- figures -------------------------------------------------------------
  # R1.7 asked for "water scarcity on the x-axis and overdraft on the y-axis
  # then show groupings of their aquifer typologies". The scatter is the main
  # Fig. 4B/C; the bar panels move to the SI so the Tukey letters stay reported.
  fig4_main_panel = "scatter",

  ## ---- parallel / caching --------------------------------------------------
  # Nothing here is expensive enough to parallelise: the whole pipeline is a
  # deterministic classification over 568 rows. Caching exists so figures and
  # tables can be re-rendered without recomputing, and is content-stamped.
  cache_check  = TRUE,
  source_files = file.path(PROJECT_ROOT, c("51_config.R", "52_core.R")),

  ## ---- frozen benchmark ----------------------------------------------------
  # TWO benchmarks, and 53_run_main.R asserts both.
  #
  #   reference_class / reference_classF   the FINALISED dataset (stage 3).
  #     This is what the preferred settings above produce and what the Results
  #     report.
  #   reference_class_published / ..._F    the PUBLISHED segment-based
  #     classification, locked against 4_Typology/SI2.csv, which an independent
  #     re-implementation reproduces 568/568 on both `class` and `classF`.
  #     53_run_main.R rebuilds it with dyad_zones = FALSE and manual_class =
  #     NULL and asserts it, so the corrections can never quietly become the
  #     only thing the pipeline can produce.
  #
  # A run that fails to reproduce either STOPS the pipeline.
  reference_n_dyads     = 568L,
  reference_n_countries = 142L,
  # ---- SI supplementary dataset (Dataset S4) --------------------------------
  # 2026-08-19. The consolidated dyad-level dataset: the old Dataset S3
  # (classification + the dyad characteristics used in the mosaics) merged
  # with the old Dataset S5 (the stage-by-stage build) on (code, CC_1, CC_2).
  # Both were 568 rows on the same key, so the merge is 1:1 and lossless. The
  # ledger's class_3 / classF_3 are dropped: they are identical to `class` /
  # `classF` by construction, which 55_summaries.R asserts rather than assumes.
  # 58_run_all.R asserts this column list, so a new column cannot appear in
  # the published dataset without being declared here.
  si_dataset_name = "Dataset_S4.csv",
  si_dataset_cols = c(
    # -- dyad identity --------------------------------------------------------
    "code", "name", "CC_1", "CC_2",
    # -- published classification (= stage 3, the finalised result) -----------
    "class", "classF",
    # -- dyad characteristics used in the Fig. S mosaics ----------------------
    "cooperation_level", "has_river_border", "dist_class",
    # -- the build: stage 1 published, stage 2 dyad-border, stage 3 = class ---
    "class_stage1", "classF_stage1", "class_stage2", "classF_stage2",
    "measured_stage2", "changed_stage2", "changed_stage3", "changed_overall",
    "near_share_stage1", "near_share_stage2",
    # -- audit provenance and written justification ---------------------------
    "audit_error_class", "audit_tier", "comment_stage2", "comment_stage3"),
  reference_n_aquifers  = 409L,
  reference_class  = c(IU0 = 87L, IU1 = 111L, UB = 205L, DA = 125L, BL = 40L),
  reference_classF = c(IU0 = 74L, IU1 = 91L,  UB = 105L, DA = 211L, BL = 87L),
  reference_class_published  =
    c(IU0 = 85L, IU1 = 118L, UB = 207L, DA = 113L, BL = 45L),
  reference_classF_published =
    c(IU0 = 70L, IU1 = 95L,  UB = 99L,  DA = 204L, BL = 100L),
  # Per-variant counts, same source and same order (IU0, IU1, UB, DA, BL).
  # 54_run_robustness.R asserts each one. These make the audit's central
  # robustness claim checkable rather than asserted: the interior/urban/
  # bilateral partition barely moves (UB 207/207/203/207 across the buffer and
  # threshold variants), while the BL-DA boundary inside the bilateral group
  # shifts with the detection and concentration thresholds exactly as a
  # tail-defined category should.
  reference_variants = list(
    stage1_published   = list(class  = c(85L, 118L, 207L, 113L, 45L),
                              classF = c(70L, 95L,  99L,  204L, 100L)),
    stage2_dyad_border = list(class  = c(98L, 118L, 201L, 112L, 39L),
                              classF = c(83L, 94L,  96L,  208L, 87L)),
    stage3_tier1       = list(class  = c(88L, 114L, 202L, 124L, 40L),
                              classF = c(74L, 94L,  103L, 210L, 87L)),
    near5    = list(class  = c(85L, 118L, 207L, 139L, 19L),
                    classF = c(70L, 95L,  99L,  238L, 66L)),
    far200   = list(class  = c(81L, 119L, 203L, 119L, 46L),
                    classF = c(66L, 96L,  95L,  210L, 101L)),
    eps_lo   = list(class  = c(76L, 104L, 153L, 141L, 94L),
                    classF = c(64L, 75L,  75L,  217L, 137L)),
    irrig    = list(class  = c(76L, 111L, 154L, 158L, 69L),
                    classF = c(65L, 95L,  91L,  216L, 101L)),
    thresh50 = list(class  = c(85L, 118L, 207L, 62L,  96L),
                    classF = c(70L, 95L,  99L,  113L, 191L)),
    bl_floor = list(class  = c(85L, 118L, 207L, 137L, 21L),
                    classF = c(70L, 95L,  99L,  224L, 80L))),
  reproduce_tol = 1e-8
)

# Resolve the active palette. Kept out of the CFG literal so the two candidate
# palettes stay side by side and switching is a one-word edit above.
CFG$class_cols <- if (identical(CFG$palette, "cvd"))
  CFG$class_cols_cvd else CFG$class_cols_published
