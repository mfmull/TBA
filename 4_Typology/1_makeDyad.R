
############################################################
# Script: Build transboundary aquifer dyads and attach metrics
#
# Purpose:
#   1. Identify country pairs that share a border segment within the
#      same transboundary aquifer (TBA)
#   2. Attach groundwater agreement scores
#   3. Join country-, aquifer-, and buffer-level indicators exported
#      from Google Earth Engine
#
# Inputs:
#   - TBA_Split_2025.shp
#   - Fraser agreement table
#   - GEE summary tables:
#       outCntr_Typol.csv
#       outTBA_Typol.csv
#       outB5k_Typol.csv
#       outB10k_Typol.csv
#       outB100k_Typol.csv
#       outB200k_Typol.csv
#       GEE SOURCE: https://code.earthengine.google.com/e32121007cd578cd8c741d9a7815c02f?noload=1
#
# Output:
#   - FinalDiadsN.csv
#
# Notes:
#   - Dyads are built only when two country polygons within the same
#     aquifer share a non-zero common boundary segment.
#   - Agreement scores are merged by aquifer-country pair and reduced
#     to the minimum score across the two sides of the dyad.
#   - GEE-derived variables are joined separately for side 1 and side 2
#     of each dyad.
############################################################

############################################################
# 0. Libraries and parallel setup
############################################################

library(sf)
library(dplyr)
library(furrr)
library(progressr)
library(countrycode)
library(tidyr)
future::plan(sequential)
plan(multisession)

############################################################
# 1. Load aquifer geometry
############################################################

tba <- st_read("../0_PreporcessIGRAC/TBA_Split_2025.shp", quiet = TRUE)

# Unique aquifer codes
codes <- unique(tba$code)

############################################################
# 2. Load GEE outputs
############################################################

buf1  <- read.csv("./geeOut/outB10k_Typol.csv")
buf2  <- read.csv("./geeOut/outB100k_Typol.csv")
buf11 <- read.csv("./geeOut/outB5k_Typol.csv")
buf22 <- read.csv("./geeOut/outB200k_Typol.csv")
cntr  <- read.csv("./geeOut/outCntr_Typol.csv")
TBAx   <- read.csv("./geeOut/outTBA_Typol.csv")

############################################################
# 3. Function to detect bordering country segments within an aquifer
############################################################

process_aquifer <- function(cd) {
  
  g <- tba %>% filter(code == cd)
  if (nrow(g) < 2) return(NULL)
  
  nb <- st_touches(g)
  res <- list()
  k <- 1
  
  for (i in seq_len(nrow(g))) {
    for (j in nb[[i]]) {
      if (i < j) {
        
        # Shared boundary length between the two country segments
        inter_len <- sum(as.numeric(
          st_length(
            st_intersection(
              st_boundary(g$geometry[i]),
              st_boundary(g$geometry[j])
            )
          )
        ))
        
        if (inter_len > 0) {
          res[[k]] <- data.frame(
            code = cd,
            CC_1 = min(g$CC[i], g$CC[j]),
            CC_2 = max(g$CC[i], g$CC[j])
          )
          k <- k + 1
        }
      }
    }
  }
  
  if (length(res) == 0) return(NULL)
  
  distinct(bind_rows(res))
}

############################################################
# 4. Build dyad table in parallel
############################################################

handlers(global = TRUE)

with_progress({
  
  p <- progressor(along = codes)
  
  diads <- future_map_dfr(
    codes,
    function(cd) {
      p(sprintf("Aquifer %s", cd))
      process_aquifer(cd)
    },
    .options = furrr_options(seed = TRUE)
  )
  
})

##attach back names
diads=diads%>%left_join(tba%>%select(code, name)%>%st_drop_geometry()%>%unique())

############################################################
# 5. Load and attach agreement scores
############################################################

agr <- read.csv("../RawData/Fraser2023_Agreements/_agreements.csv") %>%
  select(CC, code, ArrangementScore) %>%
  filter(!is.na(CC), !is.na(code), !is.na(ArrangementScore)) %>%
  group_by(CC, code) %>%
  summarise(score = max(ArrangementScore), .groups = "drop")

diads <- diads %>%
  left_join(
    agr %>% rename(CC_1 = CC, score_1 = score),
    by = c("code", "CC_1")
  ) %>%
  left_join(
    agr %>% rename(CC_2 = CC, score_2 = score),
    by = c("code", "CC_2")
  ) %>%
  mutate(ArrangementScore = coalesce(pmin(score_1, score_2), 0)) %>%
  select(-score_1, -score_2)

############################################################
# 6. Manual agreement-score corrections
# Source: Burchi 2018
# https://doi.org/10.1016/j.ejrh.2018.04.007
############################################################

diads[
  diads$code %in% c(
    "AF056", "AS126", "AF069", "AF063", "S021",
    "EU024", "AF064", "N015", "N001"
  ),
  "ArrangementScore"
] <- 3

############################################################
# 7. Helper to join side-specific GEE summaries
############################################################

join_side <- function(diadsx, df, side, src = NULL, by_code = TRUE,
                      vars = c("CR3", "GW", "GW3", "UR", "Area", 'IR','IR3')) {
  
  keys <- if (by_code) c("code", "CC") else "CC"
  
  x <- df %>%
    { if (!by_code) rename(., CC = GID_0) else . } %>%
    select(any_of(keys), any_of(vars)) %>%
    group_by(across(all_of(keys))) %>%
    summarise(across(all_of(vars), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
    rename_with(
      ~ if (is.null(src)) paste0(.x, "_", side) else paste0(.x, src, "_", side),
      all_of(vars)
    )
  
  by <- if (by_code) {
    c("code", setNames("CC", paste0("CC_", side)))
  } else {
    setNames("CC", paste0("CC_", side))
  }
  
  left_join(diadsx, x, by = by)
}

############################################################
# 8. Join GEE indicators to dyads
############################################################

diadsOut <- diads %>%
  # 10 km buffer
  join_side(buf1,  side = 1, src = "_B1",  by_code = TRUE) %>%
  join_side(buf1,  side = 2, src = "_B1",  by_code = TRUE) %>%
  
  # 100 km buffer
  join_side(buf2,  side = 1, src = "_B2",  by_code = TRUE) %>%
  join_side(buf2,  side = 2, src = "_B2",  by_code = TRUE) %>%
  
  # Aquifer-country totals
  join_side(TBAx,   side = 1, src = NULL,   by_code = TRUE) %>%
  join_side(TBAx,   side = 2, src = NULL,   by_code = TRUE) %>%
  
  # Country totals
  join_side(cntr,  side = 1, src = "_C",   by_code = FALSE) %>%
  join_side(cntr,  side = 2, src = "_C",   by_code = FALSE) %>%
  
  # 5 km buffer
  join_side(buf11, side = 1, src = "_B11", by_code = TRUE) %>%
  join_side(buf11, side = 2, src = "_B11", by_code = TRUE) %>%
  
  # 200 km buffer
  join_side(buf22, side = 1, src = "_B22", by_code = TRUE) %>%
  join_side(buf22, side = 2, src = "_B22", by_code = TRUE)

############################################################
# 9. Flag dyads associated with river borders
############################################################


LBr <- st_read("../RawData/GSRB/GSRB_Level0.shp", quiet = TRUE) |>
  st_transform(st_crs(tba)) |>
  st_make_valid()

sf::sf_use_s2(TRUE)
# tolerance in meters because sf:_use_s2 is TRUE
tol <-5000# 0.05

has_river_border <- function(code_i, cc1_i, cc2_i, segs, rivers, tol = 100) {
  g1 <- segs |>
    filter(code == code_i, CC == cc1_i)
  
  g2 <- segs |>
    filter(code == code_i, CC == cc2_i)
  
  if (nrow(g1) == 0 || nrow(g2) == 0) return(FALSE)
  
  b1 <- st_boundary(st_union(g1$geometry))
  b2 <- st_boundary(st_union(g2$geometry))
  shared <- st_intersection(b1, b2)
  
  if (length(shared) == 0 || all(st_is_empty(shared))) return(FALSE)
  
  shared_buf <- st_buffer(shared, tol)
  
  any(lengths(st_intersects(shared_buf, rivers)) > 0)
}

handlers(global = TRUE)

with_progress({
  p <- progressor(steps = nrow(diadsOut))
  
  river_flag <- furrr::future_pmap_lgl(
    list(diadsOut$code, diadsOut$CC_1, diadsOut$CC_2),
    function(code, CC_1, CC_2) {
      p(sprintf("Dyad %s: %s-%s", code, CC_1, CC_2))
      has_river_border(code, CC_1, CC_2, segs = tba, rivers = LBr, tol = tol)
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
})

diadsOut <- diadsOut |>
  mutate(has_river_border = river_flag)

############################################################
# 9. Export final dyad table
############################################################

write.csv(diadsOut, "FinalDiads.csv", row.names = FALSE)
