library(sf)
library(dplyr)
library(units)
library(rnaturalearth)
library(furrr)
library(progressr)
library(parallel)
# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script builds a “look-alike basin” candidate set for each transboundary
# aquifer polygon by searching HydroSHEDS HYBAS basins (multiple nesting levels)
# and exporting one CSV per HYBAS level with candidate matches.
# THIS IS TAKES A LONG TIME TO RUN (~45 min on my laptop)
#
# Main idea
# - Input aquifers: IGRAC transboundary aquifer polygons (TBA_Split_2025.shp).
# - Candidate units: HydroSHEDS HYBAS basins (hybas_lake_*_levXX), read for all
#   continents and filtered to exclude lake basins (LAKE == 1).
# - Exclusions: any basin that intersects (i) any aquifer polygon or (ii) the
#   Lebanon boundary polygon (LB.shp) is removed, so candidates come only from
#   “non-aquifer” areas and avoid Lebanon.
#
# Processing steps
# 1) Load and validate geometries, and transform everything to an equal-area
#    projection (EPSG:6933) for robust area and distance calculations.
# 2) For each aquifer, identify its “main country” by intersecting with Natural
#    Earth country polygons and selecting the country with the largest overlap.
# 3) For each HYBAS level requested (lev = 2…10):
#    - Read HYBAS basin shapefiles across all continents and bind them together.
#    - Compute basin areas (km^2) and assign each basin to a country using the
#      same “largest overlap” rule.
#    - For each aquifer polygon, find candidate basins that are:
#        (a) similar in size: basin area within [0.25×, 4×] aquifer area
#        (b) spatially proximate: centroid-to-centroid distance <= (1000 km + r),
#            where r is the aquifer’s equivalent-circle radius from its area.
#      The output retains basin HYBAS_ID, basin area, assigned country name, the
#      aquifer id, and the centroid distance in km.
# 4) Parallelization: matching is run in parallel over aquifers using furrr,
#    with a progress bar via progressr, then written to CSV.
#
# Outputs
# - matched_basins_out/matched_basins_lev{L}.csv for each L in {2,…,10},
#   containing one row per (aquifer, candidate basin) match plus a column `lev`.
# ------------------------------------------------------------------------------


# local matching helper
find_matches <- function(aq_row) {
  a_geom <- st_geometry(aq_row)
  a_area <- aq_row$area_km2
  a_id   <- aq_row$id
  
  cand <- bas_country %>%
    filter(area_km2 >= 0.25 * a_area, area_km2 <= 4 * a_area)
  if (nrow(cand) == 0) return(NULL)
  
  a_radius <- sqrt(a_area * 1e6 / pi)
  dist_thr <- 1e6 + a_radius
  
  # centroids
  a_cent <- st_centroid(a_geom)
  cand_cent <- st_centroid(cand)
  
  # filter by distance threshold
  within <- st_is_within_distance(cand_cent, a_cent, dist = dist_thr, sparse = FALSE)
  cand <- cand[within, ]
  if (nrow(cand) == 0) return(NULL)
  
  # compute actual centroid-to-centroid distance in km
  dist_vals <- drop_units(set_units(st_distance(cand_cent[within, ], a_cent), "km"))
  
  cand %>%
    st_drop_geometry() %>%
    select(HYBAS_ID, area_km2, name) %>%
    mutate(
      aq_id = a_id,
      dist_km = dist_vals
    )
}
sf_use_s2(FALSE)

# Read and prepare once
aq <- st_read("../../0_PreporcessIGRAC/TBA_Split_2025.shp", quiet = TRUE) %>%
  rename(aqName = name) %>%
  mutate(id = 1:nrow(.))

lb <- st_read("../../RawData/geoBoundariesCGAZ_ADM0/LB.shp", quiet = TRUE)


countries <- ne_countries(scale = "medium", returnclass = "sf")
countries <- st_transform(st_make_valid(countries), 6933)

lb <- st_transform(st_make_valid(lb), st_crs(countries))

aq <- st_transform(st_make_valid(aq), st_crs(countries))
aq_proj <- st_transform(aq, 6933)
aq_buf <- st_transform(st_buffer(aq_proj, dist = set_units(1000, "km")), st_crs(aq))

# Identify main country for each aquifer
overlaps <- st_intersection(st_make_valid(aq), st_make_valid(countries)) %>%
  mutate(intersect_km2 = as.numeric(st_area(.) / 1e6)) %>%
  group_by(id) %>%
  slice_max(intersect_km2, n = 1) %>%
  select(id, name) %>%
  as.data.frame()

aq_country <- aq %>%
  left_join(overlaps, by = "id") %>%
  mutate(area_km2 = as.numeric(st_area(.) / 1e6))

# --------------------------
# MAIN FUNCTION
# --------------------------
run <- function(lev = 6) {
  message("Running level ", lev)
  
  # read basins for all continents
  read_bas <- function(cont) {
    f <- sprintf("../../RawData/Hydrosheds/hybas_lake_%s_lev0%d_v1c/hybas_lake_%s_lev0%d_v1c.shp",
                 cont, lev, cont, lev)
    st_read(f, quiet = TRUE) %>% filter(LAKE != 1)
  }
  
  bas <- bind_rows(
    read_bas("eu"), read_bas("af"), read_bas("as"),
    read_bas("na"), read_bas("sa"), read_bas("au")
  )
  
  bas <- bas %>%
    st_transform(6933) %>%
    st_make_valid() %>%
    st_buffer(0)
    
  # drop basins intersecting aquifers
  bas <- bas[!apply(st_intersects(bas, aq, sparse = FALSE), 1, any), ]
  bas <- bas[!apply(st_intersects(bas, lb, sparse = FALSE), 1, any), ]
  
  # link to country
  overlaps <- st_intersection(st_make_valid(bas), st_make_valid(countries)) %>%
    mutate(intersect_km2 = as.numeric(st_area(.) / 1e6)) %>%
    group_by(HYBAS_ID) %>%
    slice_max(intersect_km2, n = 1) %>%
    select(HYBAS_ID, name, intersect_km2) %>%
    as.data.frame()
  
  bas_country <- bas %>%
    left_join(overlaps, by = "HYBAS_ID") %>%
    mutate(area_km2 = as.numeric(st_area(.) / 1e6))
  n_cores <- max(1, detectCores() - 1)
  plan(multisession, workers = n_cores)
  handlers(global = TRUE)
  handlers("txtprogressbar")
  
  with_progress({
    p <- progressor(along = 1:nrow(aq_country))
    matches <- future_map_dfr(
      1:nrow(aq_country),
      ~{
        # load required packages inside the worker
        suppressMessages({
          library(sf)
          library(dplyr)
          library(units)
        })
        p(sprintf("Matching %d / %d", .x, nrow(aq_country)))
        find_matches(aq_country[.x, ])
      },
      .options = furrr_options(
        globals = c("aq_country", "bas_country", "find_matches"),
        seed = TRUE
      )
    )
  })
  matches$lev=lev
  plan(sequential)
  gc()
  
  return(matches)
}
write.csv(run(4), sprintf("matched_basins_out/matched_basins_lev4.csv"), row.names = FALSE)

write.csv(run(5), sprintf("matched_basins_out/matched_basins_lev5.csv"), row.names = FALSE)

write.csv(run(6), sprintf("matched_basins_out/matched_basins_lev6.csv"), row.names = FALSE)

write.csv(run(2), sprintf("matched_basins_out/matched_basins_lev2.csv"), row.names = FALSE)
write.csv(run(3), sprintf("matched_basins_out/matched_basins_lev3.csv"), row.names = FALSE)

write.csv(run(7), sprintf("matched_basins_out/matched_basins_lev7.csv"), row.names = FALSE)
# NOT RUN
# write.csv(run(8), sprintf("matched_basins_out/matched_basins_lev8.csv"), row.names = FALSE)
# write.csv(run(9), sprintf("matched_basins_out/matched_basins_lev9.csv"), row.names = FALSE)
# write.csv(run(10), sprintf("matched_basins_out/matched_basins_lev10.csv"), row.names = FALSE)
