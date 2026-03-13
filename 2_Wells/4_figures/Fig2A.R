
library(ggplot2)
library(dplyr)
library(ggrepel)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script produces spatial visualizations of (i) Aquifer × Country first-stage
# results (global point map) and (ii) well-level groundwater depletion patterns
# within transboundary aquifers (interactive close-up maps). It combines summary
# statistics from the first-stage regressions with contextual aquifer and well
# datasets, then renders both static ggplot maps and interactive mapview layers.
#
# Part A: Global map of Aquifer × Country first-stage “distance trend” results
# Inputs
# - ../2_firststage/firstStageMain.csv
#     First-stage Aquifer × Country results (beta_1, se_1, p_1, n, TB, lat_c/lon_c)
# - ../1_data/aqfData.csv
#     Aquifer × Country reference table (ensures full coverage of Aquifer × CC units)
# - Natural Earth country polygons (rnaturalearth::ne_countries)
#
# Steps
# 1) Build a mapping dataset (d, then d_map):
#    - Select key fields from first-stage results and compute vi_1 = se_1^2.
#    - Right-join aqfData.csv to keep aquifer-country units even if some are missing
#      first-stage estimates; set TB = NA where estimates are missing.
#    - Filter to units with n > 20 wells.
# 2) Encode visual attributes:
#    - sig indicator: p_1 < 0.1
#    - color/fill categories based on sign and significance of beta_1
#    - shape categories for TB vs non-TB (and missing TB)
# 3) Plot a world basemap and overlay points at (lon_c, lat_c):
#    - A white semi-transparent “halo” layer reduces overlap ambiguity.
#    - Non-TB units plotted with size ∝ n.
#    - TB units emphasized by inflating point size (size = 5*n) and higher stroke.
#    - Coordinate window is restricted to x ∈ [-130, 130], y ∈ [-40, 60].
#    - Legends are suppressed in the final theme settings (legend.position = "none").
#
# Output (Part A)
# - A static global map showing where Aquifer × Country units have positive/negative
#   and significant/non-significant distance-to-border trends (beta_1), with TB units
#   visually emphasized.
#
# Part B: Interactive close-up maps for transboundary aquifers (mapview)
# Inputs
# - ../1_data/jasechko_aquifs/jasechko_et_al_2024_aquifers.shp
#     Aquifer polygons (filtered to only those flagged TB in firstStageMain.csv)
# - ../1_data/wellsData.csv
#     Well-level dataset with GWSlp and coordinates (lon/lat), filtered to remove
#     selected island nations
# - ../1_data/geoBoundariesCGAZ_ADM0/LB.shp
#     Land border polyline layer
# - Natural Earth country polygons (as background)
#
# Steps
# 1) Read aquifer polygons and keep only aquifers with TB == TRUE in the first-stage file.
# 2) Read wells, filter to TB aquifers, and convert to sf points (WGS84).
# 3) Winsorize and normalize depletion values per aquifer:
#    - winsorize_vec trims GWSlp to the 5th–95th percentiles within each aquifer
#      to reduce the influence of extreme wells.
#    - rescale() maps trimmed GWSlp into [0, 1] per aquifer for comparable coloring.
# 4) Build an interactive map stack (mapview):
#    - Countries as an outline-only background
#    - Aquifer polygons with dashed outlines and semi-transparent fill
#    - Land border polylines
#    - Wells colored by normalized depletion (GWSlp_norm) using a custom palette
#      and rendered with large symbols for visibility
#
# Output (Part B)
# - An interactive mapview widget allowing pan/zoom inspection of:
#    * transboundary aquifer polygons,
#    * land borders, and
#    * well points colored by within-aquifer normalized depletion.
# ------------------------------------------------------------------------------

###################
aquifersFS=read.csv('../2_firststage/firstStageMain.csv')
d=aquifersFS%>%dplyr::select(lat_c,lon_c,n,se_1,beta_1,p_1,TB,Aquifer, CC)%>%mutate(vi_1=se_1^2)%>%
  right_join(read.csv('../1_data/aqfData.csv')%>%dplyr::select(Aquifer,CC,TB, lat_c, lon_c))%>%mutate(TB=ifelse(is.na(vi_1),NA,TB))%>%
  filter(n>20)


# Assume your data are in `d_map` and already include:
# lon_c, lat_c, n, weights, color, shape

world <- ne_countries(scale = "medium", returnclass = "sf")
# --- Prepare point data ---
d_map <- d %>%
  mutate(
    sig = p_1 < 0.1,
    color = case_when(
      sig & beta_1 < 0  ~ "red",
      sig & beta_1 > 0  ~ "blue",
      TRUE              ~ "grey70"
    ),
    shape = case_when(
      is.na(TB) ~ "cross",
      TB        ~ "TBA",
      !TB       ~ "non-TBA"
    )
  )

ggplot() +
  # 1. Clean world background
  geom_sf(data = world, fill = "grey97", color = "grey90", linewidth = 0.15) +
  
  # 2. Semi-transparent white halo around points to separate overlaps
  geom_point(
    data = d_map,
    aes(x = lon_c, y = lat_c, size = n),
    color = "white", alpha = 0.25, stroke = 1.5
  ) +
  
  # 3. Main points
  geom_point(
    data = d_map%>%filter(TB==FALSE),
    aes(
      x = lon_c, y = lat_c,
      size = n,
      # alpha = w,
      fill = case_when(
        p_1 < 0.1 & beta_1 < 0 ~ "#d95f02",
        p_1 < 0.1 & beta_1 > 0 ~ "#1b9e77",
        TRUE                     ~ "white"
      ),
      shape = case_when(
        is.na(TB)  ~ "cross",
        TB         ~ "TBA",
        TRUE       ~ "non-TBA"
      )
    ),
    color = "black",
    stroke = 0.3
  ) +
  geom_point(
    data = d_map%>%filter(TB==TRUE),
    aes(
      x = lon_c, y = lat_c,
      size = 5*n,
      alpha = 0.9,
      fill = case_when(
        p_1 < 0.1 & beta_1 < 0 ~ "#d95f02",  # Teal (negative)
        p_1 < 0.1 & beta_1 > 0 ~ "#1b9e77",  # Orange (positive)
        TRUE                     ~ "white"
      ),
      shape = case_when(
        is.na(TB)  ~ "cross",
        TB         ~ "TBA",
        TRUE       ~ "non-TBA"
      )
    ),
    color = "black",
    stroke = 0.6
  )+
  # # 4. Scales and legends
  scale_fill_identity(guide = "legend", name = "Effect",
                      labels = c("Positive (p<0.1)" = "#d95f02",
                                 "Negative (p<0.1)" = "#1b9e77",
                                 "Non-significant" = "white")) +

    scale_shape_manual(
    name = "Border type",
    values = c("cross" = 4, "non-TBA" = 21, "TBA" = 22)
  ) +
  scale_size_continuous(range = c(0.8, 4), name = "Number of wells") +
  scale_alpha_continuous(range = c(0.4, 0.9), guide = "none") +
  
  # 5. Coordinate and theme
  coord_sf(xlim = c(-130, 130), ylim = c(-40, 60), expand = FALSE) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "none",
    
  )


###################Closeup maps
library(sf)
library(dplyr)
library(mapview)
library(rnaturalearth)
library(scales)

# --- 1. Read inputs
aqfsh <- st_read("../1_data/jasechko_aquifs/jasechko_et_al_2024_aquifers.shp") %>%
  filter(Aquifer %in% aquifersFS$Aquifer[aquifersFS$TB])
df=read.csv('../1_data/wellsData.csv')%>%
  mutate(LB_river=(LB_river==1))%>%filter(!CC%in%c("AUS", "NZL","TWN","JPN", "IRL", "GBR"))
dfxx <- df %>% filter(Aquifer %in% aquifersFS$Aquifer[aquifersFS$TB])
wells_sf <- st_as_sf(dfxx, coords = c("lon", "lat"), crs = 4326)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
LB <- st_read("../1_data/geoBoundariesCGAZ_ADM0/LB.shp")



winsorize_vec <- function(x, probs = c(0.05, 0.95)) {
  qs <- quantile(x, probs = probs, na.rm = TRUE)
  x[x < qs[1]] <- qs[1]
  x[x > qs[2]] <- qs[2]
  return(x)
}

# Apply per aquifer before normalization
wells_sf <- wells_sf %>%
  group_by(Aquifer) %>%
  mutate(
    GWSlp_wins = winsorize_vec(GWSlp, probs = c(0.05, 0.95)),  # trim top/bottom 5%
    GWSlp_norm = rescale(GWSlp_wins, to = c(0, 1), na.rm = TRUE)
  ) %>%
  ungroup()


pal <- colorRampPalette(c("#2b83ba", "#998ec3", "#fdae61", "#d73027"))(100)
# --- 4. Construct map


mapview(world, color = "black", lwd = 1, alpha.regions = 0, 
        layer.name = "Countries") +
  mapview(aqfsh,
          color = "grey40",       # dashed outline
          lwd = 1.2,
          lty = 2,
          alpha.regions = 0.3,    # light transparency for fill
          col.regions = "#B3CDE3",# light purple-blue fill
          layer.name = "Aquifers (light blue-purple)") +
  mapview(LB, color = "black", lwd = 2, layer.name = "Land Border") +
  mapview(
    wells_sf,
    zcol = "GWSlp_norm",
    cex =2.5,
    col.regions = pal,
    alpha = 1,
    alpha.regions = 1,  # enforce full opacity for points
    lwd = 0.1,               # thin outline (set to 0 to remove border)
    layer.name = "Wells (normalized GWSlp)"
  )

