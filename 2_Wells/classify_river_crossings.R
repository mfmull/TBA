#https://code.earthengine.google.com/ac70c49545f88173f735714d7105c215
###

# =============================================================================
# classify_river_crossings.R  --  SIMPLIFIED (2026-08-18, rev. 2)
# -----------------------------------------------------------------------------
# One question, one parameter:
#
#   Does at least one river with long-term mean discharge >= Q_MIN_CMS cross
#   into this aquifer segment EXACTLY ONCE?
#
# "Exactly once" is the whole meander/border-parallel exclusion: a river that
# runs along the border and re-crosses it repeatedly enters the segment two or
# more times and is classified as "meandering", not as a crossing.
#
# CHANGES IN REV. 2 (all in sections 4-6; sections 1-3 unchanged):
#   (a) Entries are counted on the reach NETWORK, not on merged geometry.
#       st_union() nodes a river at every confluence and st_line_merge() will
#       not merge across the junction, so a river with ONE tributary inside a
#       segment came back as 3 pieces and was misclassified as meandering.
#       Verified: a Y-shaped river with no re-entry at all returned 3 pieces.
#   (b) Upstream/downstream now compares this segment against the COUNTERPART,
#       not against the river's maximum. The maximum usually sits outside the
#       aquifer, so every side read "upstream" (both sides of Elvas did).
#   (c) n_reaches_window distinguishes "no river here" from "outside the GEE
#       export footprint" -- a FALSE with zero candidate reaches is not a
#       finding.
#   (d) Reaches are prefiltered by bbox in lon/lat before being reprojected,
#       instead of reprojecting all ~15k reaches once per aquifer.
#   (f) "Crosses once" is now measured as depth of penetration on each side
#       rather than as a piece count -- see the note above side_stats().
#   (e) The sweep writes one CSV per Q, so a change in the flagged count is
#       diagnosable.
#
# INPUT   1_data/river_crossings_raw/river_reaches.geojson
#         map/jasechko_aquifs/jasechko_CountrySplit.shp
#
# OUTPUT  1_data/river_crossings_by_unit.csv
#         1_data/river_crossings_sweep.csv   (if Q_SWEEP set)
#         1_data/river_crossings_by_unit_Q<q>.csv (one per sweep value)
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
})

sf::sf_use_s2(FALSE)


# ---- PARAMETERS -------------------------------------------------------------

Q_MIN_CMS    <- 10
Q_SWEEP      <- c(10, 30, 100, 300)   # NULL to skip

MIN_SIDE_KM  <- 5       # require at least this much channel on each side
MIN_DEPTH_KM <- 10      # ... reaching at least this far from the border on each side
SNAP_TOL_M   <- 500     # tolerate small polygon slivers/gaps
PREFILTER_M  <- 5000    # candidate-reach buffer around whole aquifer

LOGQ_SCALE   <- NA      # NA = auto-detect


ROOT <- "."

REACHES <- file.path(
  ROOT, "1_data", "river_crossings_raw", "river_reaches.geojson"
)

SEG_SHP <- file.path(
  ROOT, "map", "jasechko_aquifs", "jasechko_CountrySplit.shp"
)

OUT_UNIT  <- file.path(ROOT, "1_data", "river_crossings_by_unit.csv")
OUT_SWEEP <- file.path(ROOT, "1_data", "river_crossings_sweep.csv")


# ---- 1. segments ------------------------------------------------------------

segs <- st_read(
  SEG_SHP,
  quiet = TRUE
) %>%
  mutate(
    unit_id = paste0(Aquifer, "_", shapeGroup)
  ) %>%
  st_make_valid() %>%
  
  # Dissolve multipart pieces belonging to the same country-side segment.
  group_by(Aquifer, shapeGroup, unit_id) %>%
  summarise(.groups = "drop") %>%
  
  # Only aquifers split across >= 2 countries.
  group_by(Aquifer) %>%
  filter(n() >= 2) %>%
  ungroup()


cat(sprintf(
  "segments: %d in %d shared aquifers\n",
  nrow(segs),
  n_distinct(segs$Aquifer)
))


# ---- 2. reaches and discharge scale -----------------------------------------

reaches <- st_read(
  REACHES,
  quiet = TRUE
)

stopifnot(
  all(c("Reach_ID", "Next_down", "Log_Q_avg") %in% names(reaches))
)


# GloRiC stores Log_Q_avg log10-transformed.
# Some distributions store the log value multiplied by 10.
if (is.na(LOGQ_SCALE)) {
  
  LOGQ_SCALE <- if (max(abs(reaches$Log_Q_avg), na.rm = TRUE) > 12) 10 else 1
  
  cat(sprintf("Log_Q_avg scale auto-detected as /%d\n", LOGQ_SCALE))
}


reaches$q_cms <- 10^(reaches$Log_Q_avg / LOGQ_SCALE)

q_available_min <- min(reaches$q_cms, na.rm = TRUE)


cat(sprintf(
  "reaches: %d | implied discharge min/median/max = %.2g / %.2g / %.2g m3/s\n",
  nrow(reaches),
  min(reaches$q_cms, na.rm = TRUE),
  median(reaches$q_cms, na.rm = TRUE),
  max(reaches$q_cms, na.rm = TRUE)
))


cat("  reaches surviving each candidate floor:\n")

for (q in c(0.1, 1, 3, 10, 30, 100, 1000)) {
  
  cat(sprintf(
    "    >= %7.1f m3/s : %6d\n",
    q,
    sum(reaches$q_cms >= q, na.rm = TRUE)
  ))
}


if (Q_MIN_CMS < q_available_min) {
  
  warning(sprintf(
    paste0(
      "Q_MIN_CMS = %g is below the minimum discharge present in the reach ",
      "export (%.2g m3/s). Thresholds below %.2g m3/s therefore cannot add ",
      "any reaches."
    ),
    Q_MIN_CMS, q_available_min, q_available_min
  ))
}


# ---- 3. group reaches into rivers -------------------------------------------
#
# A "river" is an undirected connected component of the retained network.
#
# Next_down is used only to say:
#   reach A is connected to reach B.
#
# If B is absent because it was removed by the discharge threshold, that
# connection simply does not exist at this threshold.

river_id_of <- function(id, next_id) {
  
  idx <- match(next_id, id)
  
  parent <- seq_along(id)
  
  
  find <- function(i) {
    
    # Find root.
    r <- i
    
    while (parent[r] != r) {
      r <- parent[r]
    }
    
    
    # Path compression.
    #
    # find() is nested inside river_id_of(), so <<- deliberately modifies
    # river_id_of()'s local parent vector.
    while (parent[i] != r) {
      
      p <- parent[i]
      
      parent[i] <<- r
      
      i <- p
    }
    
    r
  }
  
  
  for (i in seq_along(id)) {
    
    if (is.na(idx[i])) {
      next
    }
    
    a <- find(i)
    b <- find(idx[i])
    
    
    # This assignment occurs directly inside river_id_of(), therefore it must
    # use <- rather than <<-.
    if (a != b) {
      parent[b] <- a
    }
  }
  
  
  vapply(seq_along(id), find, integer(1))
}


# ---- 4. how a river sits relative to the border ------------------------------
#
# For one connected river and one segment, measure on each side of the shared
# border: how much channel there is, and how FAR from the border it reaches.
#
#   crossing   : substantial channel reaching well inland on BOTH sides
#   following  : interacts with the border but stays in the border corridor,
#                i.e. the meandering / border-parallel case
#
# Depth of penetration is the discriminator, not a count of pieces or of
# crossing points. Two counting schemes were tried and both fail:
#
#   - merged LINESTRING pieces: st_union() nodes the network at every
#     confluence and st_line_merge() will not merge across it, so a river with
#     one tributary inside the segment returns 3 pieces and is thrown out as
#     meandering (verified on a Y with no re-entry at all);
#   - connected components of the reach network: a river meandering along the
#     border has reaches that each straddle the line, so the chain never
#     breaks and an oscillating river reads as a single clean entry (verified
#     on a 5-reach zigzag).
#
# A meander belt is a few km wide; a river that reaches MIN_DEPTH_KM into both
# national parts is connecting them whatever its planform does at the line.
# n_border_points is carried as a diagnostic only.

side_stats <- function(sub, poly, border) {
  
  clip <- suppressWarnings(st_intersection(sub, poly))
  
  if (nrow(clip) == 0) {
    return(list(len = 0, depth = 0, q = NA_real_))
  }
  
  gt <- as.character(st_geometry_type(clip))
  
  if (any(gt == "GEOMETRYCOLLECTION")) {
    clip <- suppressWarnings(st_collection_extract(clip, "LINESTRING"))
  }
  
  clip <- clip[
    as.character(st_geometry_type(clip)) %in% c("LINESTRING", "MULTILINESTRING"),
  ]
  
  if (nrow(clip) == 0) {
    return(list(len = 0, depth = 0, q = NA_real_))
  }
  
  clip$len_km <- as.numeric(st_length(clip)) / 1000
  
  len <- sum(clip$len_km)
  
  # Furthest the channel gets from the border on this side.
  depth <- if (length(border) == 0 || all(st_is_empty(border))) {
    NA_real_
  } else {
    pts <- suppressWarnings(st_cast(st_geometry(clip), "POINT"))
    max(as.numeric(st_distance(pts, border))) / 1000
  }
  
  # Length-weighted MEAN discharge, not the maximum: the reach that straddles
  # the border belongs to both sides, and taking maxima makes it the maximum on
  # both, so the sides tie and the upstream/downstream label goes missing.
  ok <- is.finite(clip$q_cms) & is.finite(clip$len_km) & clip$len_km > 0
  
  list(
    len   = len,
    depth = depth,
    q     = if (any(ok)) {
      sum(clip$q_cms[ok] * clip$len_km[ok]) / sum(clip$len_km[ok])
    } else {
      NA_real_
    }
  )
}

# ---- 5. crossing check ------------------------------------------------------

crossings_for_Q <- function(q_min) {
  
  r <- reaches[!is.na(reaches$q_cms) & reaches$q_cms >= q_min, ]
  
  
  if (nrow(r) == 0) {
    stop("no reaches at Q >= ", q_min)
  }
  
  
  # Assign retained reaches to connected river components.
  r$river <- river_id_of(r$Reach_ID, r$Next_down)
  
  
  # Work aquifer by aquifer.
  aquifer_groups <- split(seq_len(nrow(segs)), segs$Aquifer)
  
  
  out <- lapply(
    aquifer_groups,
    function(ii) {
      
      a <- segs[ii, ]
      
      
      # -----------------------------------------------------------------------
      # Cheap lon/lat bbox prefilter, so that only nearby reaches are
      # reprojected. Reprojecting the whole network once per aquifer is the
      # dominant cost otherwise.
      # -----------------------------------------------------------------------
      
      bb <- st_bbox(a)
      
      pad <- PREFILTER_M / 80000   # degrees; generous at any latitude
      
      bb[c("xmin", "ymin")] <- bb[c("xmin", "ymin")] - pad
      bb[c("xmax", "ymax")] <- bb[c("xmax", "ymax")] + pad
      
      near <- r[
        lengths(
          suppressWarnings(st_intersects(r, st_as_sfc(bb)))
        ) > 0,
      ]
      
      
      # -----------------------------------------------------------------------
      # Local metric CRS. Everything from here onward is in metres.
      # -----------------------------------------------------------------------
      
      lon0 <- mean(c(bb["xmin"], bb["xmax"]))
      lat0 <- mean(c(bb["ymin"], bb["ymax"]))
      
      crs_a <- sprintf(
        "+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m +no_defs",
        lat0, lon0
      )
      
      ap <- st_transform(a, crs_a)
      
      
      none <- tibble(
        unit_id            = ap$unit_id,
        n_reaches_window   = 0L,
        n_crossing_rivers  = 0L,
        n_meander_rivers   = 0L,
        max_q_crossing_cms = NA_real_,
        side               = NA_character_,
        min_depth_km       = NA_real_,
        n_border_points    = NA_integer_
      )
      
      if (nrow(near) == 0) {
        return(none)
      }
      
      
      rp <- st_transform(near, crs_a)
      
      rp <- rp[
        lengths(
          st_intersects(rp, st_buffer(st_union(st_geometry(ap)), PREFILTER_M))
        ) > 0,
      ]
      
      
      if (nrow(rp) == 0) {
        return(none)
      }
      
      
      # Each list element is one connected river inside this aquifer's window.
      rivers <- split(seq_len(nrow(rp)), rp$river)
      
      
      # -----------------------------------------------------------------------
      # Test each country-side segment.
      # -----------------------------------------------------------------------
      
      segment_results <- lapply(
        seq_len(nrow(ap)),
        function(k) {
          
          this <- st_geometry(ap[k, ])
          
          
          # "Other side" means all other country segments of this aquifer.
          # The small buffer prevents tiny polygon gaps/slivers from making a
          # genuine crossing look disconnected.
          other <- st_buffer(
            st_union(st_geometry(ap[-k, ])),
            SNAP_TOL_M
          )
          
          
          # The shared border: the part of this segment's boundary that adjoins
          # the counterpart. Inside an aquifer this IS the international
          # border, by construction.
          border <- suppressWarnings(
            st_intersection(st_boundary(this), other)
          )
          
          
          res <- lapply(
            rivers,
            function(jj) {
              
              sub <- rp[jj, ]
              
              
              # No shared boundary: these two parts of the aquifer do not meet,
              # so there is no border here to cross.
              if (length(border) == 0 || all(st_is_empty(border))) {
                return(NULL)
              }
              
              
              In <- side_stats(sub, this, border)
              
              if (In$len < MIN_SIDE_KM) {
                return(NULL)
              }
              
              
              Out <- side_stats(sub, other, border)
              
              if (Out$len < MIN_SIDE_KM) {
                return(NULL)
              }
              
              
              # ---------------------------------------------------------------
              # Classification.
              # ---------------------------------------------------------------
              
              is_cross <- is.finite(In$depth) &&
                is.finite(Out$depth) &&
                In$depth >= MIN_DEPTH_KM &&
                Out$depth >= MIN_DEPTH_KM
              
              
              # Upstream/downstream: discharge grows downstream, so compare
              # this segment against the counterpart. Comparing against the
              # river's own maximum is wrong -- the maximum normally lies
              # outside the aquifer altogether, which makes every side read
              # "upstream".
              side <- if (!is.finite(In$q) || !is.finite(Out$q) || In$q == Out$q) {
                NA_character_
              } else if (In$q > Out$q) {
                "downstream"
              } else {
                "upstream"
              }
              
              
              # Diagnostic only: how often this river meets the border.
              npts <- suppressWarnings(
                st_intersection(st_union(st_geometry(sub)), border)
              )
              
              n_pts <- if (length(npts) == 0) {
                0L
              } else {
                length(suppressWarnings(st_cast(npts, "POINT")))
              }
              
              
              list(
                cross   = is_cross,
                q       = max(rp$q_cms[jj], na.rm = TRUE),
                side    = side,
                depth   = min(In$depth, Out$depth),
                n_pts   = n_pts
              )
            }
          )
          
          
          # Rivers without enough channel on both sides returned NULL.
          res <- Filter(Negate(is.null), res)
          
          cr <- Filter(function(z) isTRUE(z$cross), res)
          
          
          if (length(cr) > 0) {
            
            q_cross <- vapply(cr, `[[`, numeric(1), "q")
            
            strongest <- which.max(q_cross)
            
            max_q <- max(q_cross, na.rm = TRUE)
            side  <- cr[[strongest]]$side
            depth <- cr[[strongest]]$depth
            n_pts <- cr[[strongest]]$n_pts
            
          } else {
            
            max_q <- NA_real_
            side  <- NA_character_
            depth <- NA_real_
            n_pts <- NA_integer_
          }
          
          
          tibble(
            unit_id            = ap$unit_id[k],
            n_reaches_window   = nrow(rp),
            n_crossing_rivers  = length(cr),
            n_meander_rivers   = length(res) - length(cr),
            max_q_crossing_cms = max_q,
            side               = side,
            min_depth_km       = depth,
            n_border_points    = n_pts
          )
        }
      )
      
      
      bind_rows(segment_results)
    }
  )
  
  
  bind_rows(out) %>%
    
    mutate(
      has_crossing_river = n_crossing_rivers > 0,
      q_min_cms = q_min
    ) %>%
    
    left_join(
      st_drop_geometry(segs) %>%
        select(unit_id, Aquifer, shapeGroup),
      by = "unit_id"
    ) %>%
    
    select(
      unit_id, Aquifer, shapeGroup, q_min_cms,
      has_crossing_river, n_crossing_rivers, n_meander_rivers,
      max_q_crossing_cms, side, min_depth_km, n_border_points,
      n_reaches_window
    ) %>%
    
    arrange(unit_id)
}


# ---- 6. main run ------------------------------------------------------------

cat(sprintf("\n== Q >= %g m3/s ==\n", Q_MIN_CMS))

main <- crossings_for_Q(Q_MIN_CMS)

write.csv(main, OUT_UNIT, row.names = FALSE)


cat(sprintf(
  "%d of %d segments flagged has_crossing_river (%d meandering rivers ruled out)\n",
  sum(main$has_crossing_river),
  nrow(main),
  sum(main$n_meander_rivers)
))


# A FALSE with no candidate reaches means the segment is outside the exported
# river network, not that it has no crossing river. Keep the two apart.
cat(sprintf(
  "%d segment(s) had no candidate reach in the window (not classified).\n",
  sum(main$n_reaches_window == 0)
))


cat(sprintf("wrote %s\n", OUT_UNIT))


print(
  as.data.frame(
    main[
      main$has_crossing_river,
      c("unit_id", "n_crossing_rivers", "max_q_crossing_cms", "side")
    ]
  )
)


# ---- 7. threshold sensitivity -----------------------------------------------

if (!is.null(Q_SWEEP)) {
  
  if (any(Q_SWEEP < q_available_min)) {
    
    cat(sprintf(
      paste0(
        "\nNOTE: reach export begins at %.2g m3/s; Q sweep values below this ",
        "cannot add lower-discharge reaches.\n"
      ),
      q_available_min
    ))
  }
  
  
  sweep <- bind_rows(
    lapply(
      Q_SWEEP,
      function(q) {
        
        cat(sprintf("  sweep Q >= %g m3/s\n", q))
        
        x <- crossings_for_Q(q)
        
        # One file per threshold, so that a change in the flagged count can be
        # traced to a specific aquifer rather than guessed at.
        write.csv(
          x,
          sub("\\.csv$", sprintf("_Q%g.csv", q), OUT_UNIT),
          row.names = FALSE
        )
        
        tibble(
          q_min_cms  = q,
          n_flagged  = sum(x$has_crossing_river),
          n_segments = nrow(x),
          n_meander  = sum(x$n_meander_rivers)
        )
      }
    )
  )
  
  
  write.csv(sweep, OUT_SWEEP, row.names = FALSE)
  
  cat("\nsensitivity to Q:\n")
  
  print(as.data.frame(sweep))
  
  cat(sprintf("wrote %s\n", OUT_SWEEP))
}