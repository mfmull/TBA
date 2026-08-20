###############################################################################
# Script: Split Transboundary Aquifers by Country
#
# Purpose:
# This script processes the IGRAC Transboundary Aquifer (TBA) dataset and
# splits each aquifer polygon by national boundaries. The resulting segments
# represent the portion of each aquifer located within each country.
#
# Workflow overview:
# 1. Load aquifer and country boundary datasets.
# 2. Ensure geometries are valid for spatial operations.
# 3. Identify which countries intersect each aquifer polygon.
# 4. Perform parallelized intersection to split aquifers by country.
# 5. Export segments with potential inconsistencies for manual inspection.
# 6. Integrate manual decisions on disputed segments.
# 7. Dissolve segments and produce:
#       - SI1: documentation table of modifications
#       - tbaFinal: cleaned aquifer-country spatial dataset
#
# Key outputs:
# - discrepancies.shp: segments requiring manual inspection
# - SI1.csv: documentation table describing segment decisions
# - TBA_Split_2025N: final aquifer-country shapefile
#
# Author: [Add author]
# Date: [Add date]
###############################################################################

# --- Load required packages --------------------------------------------------

library(sf)
library(lwgeom)         # robust geometry operations (e.g., st_geod_area)
library(dplyr)
library(future.apply)   # parallel map utilities
library(units)
library(stringr)
library(furrr)
library(progressr)

# --- Input data paths (edit if necessary) ------------------------------------

# IGRAC Transboundary Aquifer shapefile
aq <- st_read("../RawData/IGRAC/TBA2025/tba_map_2025.shp")

# Ensure aquifer geometries are valid before spatial operations
aq <- st_make_valid(aq)

# Country boundaries from geoBoundaries
cn <- st_read("../RawData/geoBoundariesCGAZ_ADM0/geoBoundariesCGAZ_ADM0.shp") %>%
  rename(CC = shapeGroup, CountryName = shapeName)

# Ensure country geometries are valid
cn <- st_make_valid(cn)


###############################################################################
# STEP 1 — Split aquifers by countries
###############################################################################

# Assign stable ID to each aquifer polygon
aq <- aq %>% mutate(aq_id = row_number())

# Determine which countries intersect each aquifer
idx <- st_intersects(aq, cn, sparse = TRUE)

# Setup parallel backend using available CPU cores
plan(multisession, workers = parallel::detectCores() - 1)

# Enable progress bar
handlers(global = TRUE)

# Perform parallel intersection
with_progress({
  
  # Progress tracker
  p <- progressor(along = seq_len(nrow(aq)))
  
  aq_split <- future_map_dfr(seq_len(nrow(aq)), function(i) {
    
    p()  # update progress bar
    
    # Countries intersecting aquifer i
    hits <- idx[[i]]
    
    # Skip if no overlap
    if (!length(hits)) return(NULL)
    
    # Subset relevant country attributes
    csub <- cn[hits, c("CC", "CountryName"), drop = FALSE]
    
    # Intersect aquifer polygon with country boundaries
    piece <- st_intersection(aq[i, , drop = FALSE], csub)
    
    # Return intersection results only if geometries exist
    if (nrow(piece)) piece else NULL
  })
})

# Ensure resulting geometries remain valid
aq_split <- st_make_valid(aq_split)


###############################################################################
# STEP 2 — Quality control of aquifer-country segments
###############################################################################
# Identify inconsistencies between:
#   - Country name from spatial overlay
#   - Country membership listed in IGRAC attributes
#
# These cases are exported for manual inspection.

discrepancies <- aq_split %>%
  filter(!str_detect(countries, fixed(CountryName)))

# Export problematic segments for manual review
st_write(
  discrepancies,
  "./Discrepancies/discrepancies0.shp",
  delete_layer = FALSE
)

###############################################################################
# MANUAL STEP REQUIRED
###############################################################################
# Review the exported shapefile and edit it manually.
#
# Required fields in the edited file:
#   CC        : country code
#   code      : aquifer code
#   comment   : decision on the segment
#
# Example values for "comment":
#   kept              -> segment should be included in dataset
#   <text reason>     -> explanation for discarding segment
###############################################################################


###############################################################################
# STEP 3 — Integrate manual decisions
###############################################################################

library(sf)
library(dplyr)
library(stringr)

# Read edited discrepancy file
discr <- st_read("./Discrepancies/discrepancies.shp", quiet = TRUE)

# Remove geometry and keep decision attributes
discr_tbl <- discr %>%
  st_drop_geometry() %>%
  select(CC, code, comment)

# Join decisions back to the spatial segments
tba <- aq_split %>%
  left_join(discr_tbl, by = c("CC", "code"))

# Recode disputed territories to Palestine (PSE)
tba <- tba %>%
  mutate(CC = if_else(CC %in% c("118", "129"), "PSE", CC))

# Dissolve duplicate aquifer-country segments
tba <- tba %>%
  group_by(code, CC) %>%
  summarise(
    across(-geometry, dplyr::first),
    geometry = sf::st_union(geometry),
    .groups = "drop"
  ) %>%
  st_as_sf()


###############################################################################
# STEP 4 — Produce SI1 documentation table
###############################################################################
# SI1 records how segments were handled:
#   - Original IGRAC entry
#   - Added to IGRAC
#   - Discarded segments with explanation

SI1 <- tba %>%
  st_drop_geometry() %>%
  mutate(
    Comment = case_when(
      comment == "kept" ~ "Added to IGRAC",
      !is.na(comment) & comment != "kept" ~ paste("Discarded:", comment),
      TRUE ~ "Original IGRAC entry"
    )
  ) %>%
  rename(IGRAC.region = region) %>%
  select(aq_id, CC, CountryName, countries, code, name, IGRAC.region, Comment)

# Export documentation table
write.csv(SI1, "./SI1.csv", row.names = FALSE)


###############################################################################
# STEP 5 — Create final spatial dataset
###############################################################################
# Keep:
#   - Original IGRAC entries
#   - Segments explicitly marked as "kept"

tbaFinal <- tba %>%
  filter(is.na(comment) | comment == "kept") %>%
  select(aq_id, CC, CountryName, countries, code, name, region, geometry) %>%
  rename(
    CntryNm = CountryName,
    countrs = countries,
    Region = region
  ) %>%
  st_make_valid()

# Export final aquifer-country dataset
st_write(tbaFinal, "TBA_Split_2025N", delete_layer = TRUE)