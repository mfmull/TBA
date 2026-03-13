# ==============================================================================
# Script purpose (preamble)
# ==============================================================================
# This script builds a *control-basin* polygon dataset (HydroSHEDS HYBAS basins)
# matched to each IGRAC transboundary aquifer, for downstream comparisons and for
# ingestion in Google Earth Engine (GEE).
#
# In brief, it:
#   (1) Loads precomputed aquifer→candidate-basin matches for multiple HYBAS
#       levels (one CSV per level).
#   (2) Chooses a single “best” HYBAS level per aquifer using an area-similarity
#       criterion (aquifer area vs median candidate basin area at each level).
#   (3) Loads HydroSHEDS basin geometries only for the selected levels.
#   (4) For each aquifer, extracts candidate basins at its selected level and
#       keeps up to NMax basins, prioritizing same-country candidates and then
#       the closest remaining candidates by centroid distance.
#   (5) Exports:
#       - corrTable.csv: aquifer↔basin correspondence table
#       - controlBas2025.shp: deduplicated basin polygons (unique HYBAS_ID)
#
# IMPORTANT:
#   The exported shapefile `matched_basins_out/controlBas2025.shp` must be
#   uploaded to GEE. Your GEE scripts then use this control-basin layer as `bas`.
#
# Inputs
#   A) Candidate match tables (one per HYBAS level):
#      - matched_basins_out/matched_basins_lev{L}.csv
#        Required fields (used downstream in this script):
#          aq_id, HYBAS_ID, area_km2, dist_km, name (country label), ...
#
#   B) Aquifer polygons:
#      - ../../RawData/IGRAC/TBA_Split_2025.shp
#
#   C) HydroSHEDS basin polygons (by level and continent):
#      - ../../RawData/Hydrosheds/hybas_lake_{cont}_lev0{L}_v1c/*.shp
#        cont ∈ {eu, af, as, na, sa, au}
#
#   D) Land-with-borders mask (filters out basins not intersecting relevant land):
#      - ../../RawData/geoBoundariesCGAZ_ADM0/LandWithBorders_rough.shp
#
#   E) Countries (loaded for consistent CRS preparation; not essential for outputs):
#      - rnaturalearth::ne_countries()
#
# Outputs
#   1) matched_basins_out/corrTable.csv
#      - Mapping from aquifer id → retained HYBAS_IDs (plus selected level and distance).
#
#   2) matched_basins_out/controlBas2025.shp
#      - Deduplicated (unique HYBAS_ID) basin polygons in EPSG:4326 with centroid lon/lat.
#      - This is the file to upload to GEE.
#
# CRS conventions
#   - EPSG:6933 (equal-area, meters) is used for area and distance computations.
#   - EPSG:4326 (WGS84) is used for the final shapefile to ensure GEE compatibility.
#
# Workflow outline
#   1) Read and stack all matched_basins_lev*.csv; infer `level` from filename.
#   2) Read aquifers, compute aquifer areas (km²) in EPSG:6933.
#   3) For each aquifer×level, compute median candidate-basin area; choose the
#      level minimizing |median_basin_area − aquifer_area|.
#   4) Read HydroSHEDS basins only for the selected levels (all continents) and
#      reproject/validate.
#   5) For each aquifer:
#        - subset candidate HYBAS_IDs at its chosen level
#        - compute centroid distance to aquifer (km)
#        - keep up to NMax, prioritizing same-country (if `name` exists) then closest
#   6) Join candidate metadata, write corrTable.csv, deduplicate HYBAS_IDs, add lon/lat,
#      filter to LandWithBorders, and write controlBas2025.shp
# ==============================================================================
############################################################
# 0. Load packages
############################################################
library(dplyr)
library(readr)
library(purrr)
library(sf)
library(progressr)


####


############################################################
# 1. Load all matched CSVs (levels)
############################################################
# Each CSV corresponds to matched basins at one hydroshed level
csv_files <- list.files(
  path = "matched_basins_out",
  pattern = "^matched_basins_lev\\d+\\.csv$",
  full.names = TRUE
)

# Combine all levels into one table
# read and combine
matches_all <- csv_files %>%
  set_names() %>%
  map_dfr(readr::read_csv, .id = "source") %>%
  mutate(level = as.integer(sub(".*lev([0-9]+)\\.csv$", "\\1", basename(source)))) %>%
  select(-source)



############################################################
# 2. Load aquifers and countries
############################################################
aq <- st_read("../../RawData/IGRAC/TBA_Split_2025.shp", quiet = TRUE) %>%
  rename(aqName = name) %>%
  mutate(id = 1:nrow(.)) %>%
  st_make_valid()

countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Reproject everything to a projected CRS (meters)
crs_proj <- 6933
aq <- st_transform(aq, crs_proj)
countries <- st_transform(countries, crs_proj)


############################################################
# 3. Compute area per aquifer and identify optimal basin level
############################################################
# Compute aquifer areas (km²)
aq_area <- aq %>%
  mutate(area_km2 = as.numeric(st_area(geometry) / 1e6)) %>%
  st_drop_geometry() %>%
  select(id, aq_area_km2 = area_km2)

# Median basin area by aquifer × level
level_summary <- matches_all %>%
  group_by(aq_id, level) %>%
  summarise(med_basin_area = median(area_km2, na.rm = TRUE), .groups = "drop")

# Select the hydroshed level with basin areas closest to each aquifer’s area
aq_best_level <- level_summary %>%
  left_join(aq_area, by = c("aq_id" = "id")) %>%
  mutate(diff = abs(med_basin_area - aq_area_km2)) %>%
  group_by(aq_id) %>%
  slice_min(diff, n = 1, with_ties = FALSE) %>%
  ungroup()


############################################################
# 4. Load Hydrosheds basins for relevant levels
############################################################
get_bas_path <- function(level, continent) {
  sprintf("../../RawData/Hydrosheds/hybas_lake_%s_lev0%d_v1c/hybas_lake_%s_lev0%d_v1c.shp",
          continent, level, continent, level)
}

# Read basins for one level (with console feedback)
read_bas_level <- function(level) {
  conts <- c("eu", "af", "as", "na", "sa", "au")
  message("\n🗺️  Reading basins for level ", level)
  start <- Sys.time()
  
  bas_list <- vector("list", length(conts))
  for (i in seq_along(conts)) {
    cont <- conts[i]
    f <- get_bas_path(level, cont)
    message("  [", i, "/", length(conts), "] ", cont, ": ", basename(f))
    if (!file.exists(f)) {
      warning("     → file missing, skipping")
      next
    }
    bas_list[[i]] <- st_read(f, quiet = FALSE)
  }
  
  bas <- bind_rows(bas_list)
  bas <- st_make_valid(st_transform(bas, crs_proj))
  message("✅ Finished level ", level, " in ",
          round(difftime(Sys.time(), start, units = "secs"), 1), " sec")
  bas
}

# Read only the needed levels
levels_needed <- sort(unique(aq_best_level$level))

handlers(global = TRUE)
handlers("txtprogressbar")

with_progress({
  p <- progressor(along = levels_needed)
  bas_by_level <- lapply(levels_needed, function(lv) {
    p(sprintf("Reading Hydrosheds level %d", lv))
    read_bas_level(lv)
  })
})
names(bas_by_level) <- levels_needed


############################################################
# 5. Match aquifers to their selected basins (one level per aquifer)
############################################################
n <- nrow(aq_best_level)
bas_list <- vector("list", n)
NMax=40
pb <- txtProgressBar(min = 0, max = n, style = 3)

for (i in seq_len(n)) {
  aq_id_i <- aq_best_level$aq_id[i]
  level_i <- aq_best_level$level[i]
  this_level <- as.character(level_i)
  
  bas <- bas_by_level[[this_level]]
  
  # HYBAS_IDs matched to this aquifer and level
  target_ids <- matches_all$HYBAS_ID[
    matches_all$aq_id == aq_id_i &
      matches_all$level == level_i
  ]
  
  if (length(target_ids) > 0) {
    # Subset basins and annotate with aquifer info
    bas_sub <- bas %>%
      filter(HYBAS_ID %in% target_ids) %>%
      mutate(aq_id = aq_id_i,
             lev   = level_i)
    
    # Get aquifer geometry and main country
    aq_geom <- aq[aq$id == aq_id_i, ]
    aq_country <- aq_geom$name[1]
    
    # Compute centroid distances (in meters)
    dist_vec <- as.numeric(
      st_distance(st_centroid(bas_sub), st_centroid(aq_geom))
    )
    bas_sub$dist_km <- dist_vec / 1000
    
    # Filter by country first (if available)
    if ("name" %in% names(bas_sub)) {
      same_country <- bas_sub %>%
        filter(name == aq_country)
    } else {
      same_country <- bas_sub[0, ]
    }
    
    # Prioritize same-country basins, then closest others
    if (nrow(same_country) >= NMax) {
      bas_keep <- same_country %>%
        slice_min(dist_km, n = NMax)
    } else {
      remaining <- bas_sub %>%
        filter(!(HYBAS_ID %in% same_country$HYBAS_ID)) %>%
        slice_min(dist_km, n = NMax - nrow(same_country))
      bas_keep <- bind_rows(same_country, remaining)
    }
    
    bas_list[[i]] <- bas_keep
  } else {
    bas_list[[i]] <- NULL
  }
  
  print(paste("Aquifer", aq_id_i, ":", nrow(bas_list[[i]]), "basins kept"))
  setTxtProgressBar(pb, i)
}

close(pb)

# Combine all basins into a single sf object
bas_country <- bas_list %>%
  compact() %>%
  bind_rows()%>%
    left_join(matches_all %>% select(HYBAS_ID, aq_id, level, dist_km, area_km2, name),
              by = c("HYBAS_ID", "aq_id", "lev" = "level"))%>%
  select(HYBAS_ID,aq_id,lev,dist_km.x,area_km2,name,geometry)%>%rename(CntrName=name,dist_km=dist_km.x)

#That's the correspondance table if needed
corrTable=bas_country%>%as.data.frame()%>%select(HYBAS_ID,aq_id,lev,dist_km)
write.csv(corrTable,'matched_basins_out/corrTable.csv')
  


# Ensure CRS is WGS84 (EPSG:4326) for Earth Engine
#That's the shp without duplicate, which we can perhaps use as input for our matching
bas_country=bas_country%>%select(HYBAS_ID,lev,CntrName,area_km2,geometry)
bas_country=  bas_country[!duplicated(bas_country$HYBAS_ID), ]
bas_country <- st_transform(bas_country, 4326)%>%st_make_valid()
bas_country <- bas_country %>%
  mutate(HYBAS_ID = as.character(HYBAS_ID),
         lon = st_coordinates(st_centroid(geometry))[, "X"],
         lat = st_coordinates(st_centroid(geometry))[, "Y"])
# Write shapefile (creates .shp, .dbf, .shx, .prj)



#Filter out Hybas that are not on a land-border containing land
landWithBorders=st_transform(st_read("../../RawData/geoBoundariesCGAZ_ADM0/LandWithBorders_rough.shp"), st_crs(bas_country))
bas_country <- bas_country[ lengths(st_intersects(bas_country, landWithBorders)) > 0, ]



st_write(bas_country, "matched_basins_out/controlBas2025.shp", delete_dsn = TRUE)




