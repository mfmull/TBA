library(tibble)
library(dplyr)
library(ggplot2)
###################Data sources:
#GEE Map (for maing part of figure 1)
#https://code.earthengine.google.com/fdf7e7e39ec8e3450c30f3ce5ea89ba2
#GEE data
#https://code.earthengine.google.com/a83f29a381d3b721f5ab476ab5a06500?noload=1
# inputs
####################
# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script constructs a simple stacked bar chart summarizing irrigation area
# decomposition from Google Earth Engine (GEE) outputs. It splits:
#   (1) total irrigation area into surface-water (Sw) vs groundwater (Gw), and
#   (2) groundwater irrigation area into non-transboundary-aquifer (non TBA) vs
#       transboundary-aquifer (TBA) components,
# then plots the resulting areas in million km².
#
# Inputs
# - MnIr:     mean fraction (or share) of total area under irrigation
# - MnIrGw:   mean fraction of total area under groundwater irrigation
# - MnIrGwTba:mean fraction of total area under groundwater irrigation in TBAs
# - PolygonArea: total study polygon area, scaled to million km²
#   (here computed from an m² value divided by 1e12)
#
# Consistency checks
# - stopifnot(MnIr >= MnIrGw, MnIrGw >= MnIrGwTba)
#   Ensures nested shares are logically ordered:
#     total irrigation ≥ groundwater irrigation ≥ groundwater-in-TBA
#
# Data construction
# - Build a tidy tibble (plot_df) with two bars, each decomposed into two stacks:
#   Bar 1: "Total irrigation (Ir)"
#     * Sw area  = (MnIr - MnIrGw) * PolygonArea
#     * Gw area  = MnIrGw * PolygonArea
#   Bar 2: "Groundwater irrigation (Gw)"
#     * non TBA  = (MnIrGw - MnIrGwTba) * PolygonArea
#     * TBA      = MnIrGwTba * PolygonArea
# - Convert both Bar and Component to factors to control ordering in the plot.
#
# Plot
# - ggplot2 stacked columns (geom_col) with:
#   * x-axis: bar category (total irrigation vs groundwater irrigation)
#   * y-axis: area in million km²
#   * fill: component category (Sw/Gw or non TBA/TBA)
# - Manual colors distinguish components:
#   * Sw in light grey
#   * Gw and non TBA in blue
#   * TBA in green
#
# Output
# - A stacked bar chart (printed as `p`) showing the irrigation area breakdown in
#   million km², based on the provided GEE-derived shares and total polygon area.
# ------------------------------------------------------------------------------



MnIr=0.010395152045752865
MnIrGw=0.004230504800487288
MnIrGwTba=0.0014572515707868263
PolygonArea=364141029582227.75/1e12  # million km^2

stopifnot(MnIr >= MnIrGw, MnIrGw >= MnIrGwTba)

plot_df <- tibble(
  Bar = c("Total irrigation (Ir)", "Total irrigation (Ir)",
          "Groundwater irrigation (Gw)", "Groundwater irrigation (Gw)"),
  Component = c("Sw", "Gw",
                "non TBA", "TBA"),
  Area = c(MnIr - MnIrGw, MnIrGw,
           MnIrGw - MnIrGwTba, MnIrGwTba) * PolygonArea
) %>%
  mutate(
    Bar = factor(Bar, levels = c("Total irrigation (Ir)", "Groundwater irrigation (Gw)")),
    Component = factor(Component, levels = c("Sw", "Gw", "non TBA", "TBA"))
  )

p=ggplot(plot_df, aes(x = Bar, y = Area, fill = Component)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c(
    "Sw"       = "grey80",
    "Gw"                = "#1f78b4",
    "non TBA"      = "#1f78b4",
    "TBA"= "#33a02c"
  )) +
  labs(x = NULL, y = "Area (million km²)", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

print(p)