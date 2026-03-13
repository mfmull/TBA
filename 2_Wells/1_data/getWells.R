# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script assembles a cleaned, spatially-enriched well dataset (from Google
# Earth Engine exports) and links each well to (i) a country, (ii) the nearest
# land-border segment (classified as river border vs non-river border), and
# (iii) a “Nature” aquifer polygon (Jasechko et al. 2024), then exports both
# well-level and aquifer-level summary tables.
#
#GEE SCRIPTS TO RUN BEFORE
#     Wells: https://code.earthengine.google.com/6ed1434fbca51b7acaeb2658f048d5e4?noload=1
#     Well Aquifers: https://code.earthengine.google.com/9661594f3ec29ae3f375764555334fea?noload=1
#
# Inputs
# - GEE well attribute exports split across three CSV parts:
#     GEE_Out/Wells_withAttr_part1.csv
#     GEE_Out/Wells_withAttr_part2.csv
#     GEE_Out/Wells_withAttr_part3.csv
# - Country boundaries (geoBoundaries ADM0)
# - Land border lines (LB.shp)
# - River-border lines (GSRB_Level0.shp) used to label which border segments are rivers
# - Aquifer polygons (jasechko_et_al_2024_aquifers.shp) + GEE aquifer attributes:
#     GEE_Out/AqfWells_withAttr.csv
#
# Main processing steps
# 1) Merge the three well CSV parts into one table, create an id, and convert to
#    an sf point object (WGS84 lon/lat). Rename the groundwater depletion slope
#    variable and later convert it to mm/year.
# 2) Assign each well to a country (CC) by spatial intersection with ADM0
#    polygons; for wells that fall outside polygons (CC is NA), assign the
#    nearest country polygon instead.
# 3) Classify border segments as “river border” or not:
#    - Reproject land borders (LB) and river borders (GSRB) to EPSG:3857.
#    - Buffer river borders by 5 km and split LB into two parts:
#        * segments within 5 km of a river border (LB_river = 1)
#        * remaining segments (LB_river = 0)
# 4) Compute well distance to the nearest border segment:
#    - Reproject wells and LB_split to EPSG:3857.
#    - For each well, find the nearest border segment, store:
#        * dist_LB_km: nearest distance in km
#        * LB_river: whether that nearest segment is a river border
# 5) Link wells to Nature aquifers (Jasechko et al. 2024):
#    - Read and validate aquifer polygons; mark aquifers as transboundary (TB)
#      if they intersect land borders.
#    - Join additional aquifer attributes from the GEE export.
#    - Spatially join wells to aquifer polygons (keeping only intersecting wells).
#      A fallback nearest-feature assignment is attempted for unmatched wells.
# 6) Export well-level products:
#    - Create a flat CSV (wellsData.csv) with lon/lat columns and no geometry.
#    - Write a variable-name “readme” CSV (ReadMeWells.csv).
# 7) Create aquifer-level summaries:
#    - Count wells per aquifer and compute mean dist_LB_km and mean GW slope.
#    - Intersect aquifers with country polygons to create aquifer-by-country
#      segments and compute summary statistics per segment.
#    - Compute each segment centroid distance to the nearest land border.
#    - Export the final aquifer-by-country attribute table (aqfData.csv) and a
#      variable-name “readme” (ReadMeAqf.csv).
#
# Outputs
# - wellsData.csv: well-level dataset with country code, aquifer name, TB flag,
#   nearest-border distance (km), river-border indicator, and covariates.
# - ReadMeWells.csv: list of well-level variable names.
# - aqfData.csv: aquifer-by-country segment summary dataset (attributes averaged
#   within each Aquifer×CC segment) plus centroid distance to land border.
# - ReadMeAqf.csv: list of aquifer-segment variable names.
# ------------------------------------------------------------------------------

library(tidyverse)
library(sf)
library(mapview)
library(terra)
library(raster)
library(nngeo)
#Those are wells with sampled values for covariates from here:
#Wells: https://code.earthengine.google.com/6ed1434fbca51b7acaeb2658f048d5e4?noload=1
#Well Aquifers: https://code.earthengine.google.com/9661594f3ec29ae3f375764555334fea?noload=1
#Open all wells run once
wells=read.csv('GEE_Out/Wells_withAttr_part1.csv')
wells2=read.csv('GEE_Out/Wells_withAttr_part2.csv')
wells3=read.csv('GEE_Out/Wells_withAttr_part3.csv')
wells=wells%>%dplyr::select(-.geo)%>%
  left_join(wells2%>%dplyr::select(-.geo,-lat,-lon))%>%
  left_join(wells3%>%dplyr::select(-.geo,-lat,-lon))%>%dplyr::select(-system.index)
# save(wells, file='wellsRAW.rdata')
# load('wellsRAW.rdata')

wells=wells%>%
  mutate(id=1:nrow(.))%>%rename(GWSlp=Slope_mPerY_Since2000)%>%
  st_as_sf(coords=c('lon','lat'),crs=4326)

#Overlap with country
cntr <- st_read("../../RawData/geoBoundariesCGAZ_ADM0/geoBoundariesCGAZ_ADM0.shp")%>%
  st_make_valid(.)
wells <- st_join(
  wells,
  cntr %>% dplyr::select(CC = shapeGroup),
  join = st_intersects,
  left = TRUE
)

# 2. Find wells without intersection (CC = NA)
no_hit <- which(is.na(wells$CC))

if (length(no_hit) > 0){
  nearest_idx <- st_nearest_feature(wells[no_hit, ], cntr)
  # assign nearest info
  wells$CC[no_hit] <- cntr$shapeGroup[nearest_idx]
}



#Label borders are rivers or not:
# 2️⃣ Reproject both to a metric CRS (e.g., EPSG:3857 for global data)
LB <- st_read("../../RawData/geoBoundariesCGAZ_ADM0/LB.shp")
LBr <- st_read("../../RawData/GSRB/GSRB_Level0.shp")#river borders
LB  <- st_transform(LB, 3857)
LBr <- st_transform(LBr, 3857)
LBr_buf <- st_buffer(LBr, dist = 5000)   # 5000 m = 5 km
LB_cut <- st_intersection(LB, st_union(LBr_buf))
LB_cut <- LB_cut %>%
  mutate(LB_river = 1)
LB_diff <- st_difference(LB, st_union(LBr_buf)) %>%
  mutate(LB_river = 0)
LB_split <- bind_rows(LB_cut, LB_diff)


#Distance to LB
# sf_use_s2(FALSE)  # disables S2 and frees its cache # 
LB_split=LB_split%>%
  st_transform(3857) %>%          # meters, safe globally
  st_simplify(dTolerance = 100)%>%   # 10
  st_make_valid(.)#%>%st_simplify(dTolerance = 0.001)
wells=st_transform(wells, 3857)
nearest_id <- st_nearest_feature(wells, LB_split)
wells$dist_LB_km <- as.numeric(st_distance(wells, LB_split[nearest_id, ], by_element = TRUE))/1000
wells$LB_river <- as.data.frame(LB_split)[nearest_id,'LB_river' ]



#Determine Overlap with Nature Aquifers
###Open Nature Aquifers, determine overlap with land borders
aqf <- st_read("../../RawData/jasechko_aquifs/jasechko_et_al_2024_aquifers.shp")%>%
  st_transform(3857)%>%
  st_make_valid(.)
aqf$TB <- lengths(st_intersects(aqf, LB)) > 0

aqfDat=read.csv('GEE_Out/AqfWells_withAttr.csv')%>%dplyr::select(-system.index,-.geo)
aqf=aqf%>%left_join(aqfDat)


#Overlap with wells: Add AQF, TB; duplicate if multiple aquifers. There are no duplicated.
wells <- st_join(
  wells, aqf[, c("Aquifer", "TB")],
  join = st_intersects,
  left = FALSE
)

# 2. Find wells without intersection (CC = NA)
no_hit <- which(is.na(wells$Aquifer))

if (length(no_hit) > 0){
  nearest_idx <- st_nearest_feature(wells[no_hit, ], cntr)
  # assign nearest info
  wells$Aquifer[no_hit] <- aqf$Aquifer[nearest_idx]
  wells$TB[no_hit] <- aqf$TB[nearest_idx]
}



wells_df <- wells %>%  st_transform(4326) %>%  # project to WGS84 lat/lon
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()


wells_df=wells_df%>%mutate(GWSlp=GWSlp*1000)#mm per year

write.csv(wells_df,'wellsData.csv')
cnames <- names(wells_df)
write.csv(data.frame(Variable=cnames), "ReadMeWells.csv", row.names = FALSE)



#AQF: count wells,
aqf=aqf%>%
  left_join(wells_df%>%group_by(Aquifer)%>%summarize(nW=n(),mean_dist_LB_km=mean(dist_LB_km),mean_GWSlp=mean(GWSlp)))


#Intersect with cntr to make segment
# Ensure both are in the same CRS
aqf_cntr <- st_intersection(
  aqf,
  cntr %>%
    st_transform(st_crs(aqf)) %>%
    st_make_valid() %>%
    dplyr::select(CC = shapeGroup)
) %>%
  group_by(Aquifer, CC) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    across(!where(is.numeric), ~ first(.x)),
    .groups = "drop"
  ) %>%
  st_make_valid()

#centroid dist to border
aqf_cent <- st_centroid(aqf_cntr)
nearest_id <- st_nearest_feature(aqf_cent, LB)
aqf_cntr$dist_to_LB_km <- as.numeric(
  st_distance(aqf_cent, LB[nearest_id, ], by_element = TRUE)
) / 1000


#Final cleanup
aqf_cntr_db=aqf_cntr%>%
  st_drop_geometry()%>%dplyr::select(-constant)


write.csv(aqf_cntr_db,'aqfData.csv')


cnames <- names(aqf_cntr_db)
write.csv(data.frame(Variable=cnames), "ReadMeAqf.csv", row.names = FALSE)


