# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script prepares supplementary-information (SI) data products from the
# Aquifer × Country first-stage results and generates a LaTeX summary-statistics
# table comparing transboundary versus non-transboundary units, including a
# matched-weighted non-TB comparison group.
#
# Inputs
# - ../2_firststage/firstStageMain.csv
#     Aquifer × Country first-stage dataset containing regression outputs
#     (beta/se/p for two models), aquifer/country descriptors, and matching
#     covariates/weights.
#
# Outputs (data export)
# 1) SI_Data_Wells.csv
#    - A curated subset of columns from firstStageMain.csv, renamed to clearer,
#      publication-facing variable names (e.g., nWells, beta_AveDepletion,
#      beta_distTrend, CropSuitIndex, probBorderRiver, TransboundaryAquifer).
#    - Intended as SI data accompanying the analysis.
#
# Helper functions
# - classify_slope(beta, p):
#     Categorizes a coefficient by sign and significance at p < 0.1
#     (positive/negative × significant/non-significant).
# - weighted_mean(), weighted_se():
#     Convenience functions for weighted summaries (not used in the final table
#     below; the table uses weighted quantiles instead).
#
# Summary-statistics table
# 1) Load the first-stage dataset and create slope-class labels (class_0/class_1)
#    based on beta_1 and p_1 for later grouping/diagnostics.
# 2) Define a set of variables to summarize (summarize_vars), spanning:
#    - covariates used for matching (e.g., CS_max, urbkHaKm2, LB_river, STS),
#    - geography/context (e.g., RDS_mainDist, dist_to_LB_km),
#    - outcome-related quantities (e.g., m_per_year),
#    - unit composition (e.g., nW).
# 3) summary_stats(data, w = NULL):
#    - Computes median and interquartile range [Q1, Q3] for each variable.
#    - If weights are provided, uses Hmisc::wtd.quantile to compute weighted
#      median and weighted quartiles; otherwise uses unweighted quantiles.
#    - Returns formatted strings "median [Q1, Q3]" for each variable.
# 4) Build a three-column comparison table:
#    - Non-Transboundary (unweighted)
#    - Matched non-TB (weighted by matching weights)
#    - Transboundary (unweighted)
# 5) Transpose the summary output so variables are rows and groups are columns,
#    then render as a LaTeX table using knitr::kable and kableExtra styling.
#
# Final output
# - A LaTeX-formatted table titled:
#   "Summary statistics by transboundary status, including matched non-TB sample:
#    Median and [Quartiles]"
#   suitable for inclusion in SI/appendix material.
# ------------------------------------------------------------------------------
# ============================================================
# Setup
# ============================================================
library(dplyr)
library(ggplot2)
library(knitr)
library(kableExtra)


####That's SI data
xx=read.csv('../2_firststage/firstStageMain.csv')%>%
  select(Aquifer,CC,area_km2,n,beta_0,se_0,p_0,beta_1,se_1,p_1,lat_c,lon_c,urbkHaKm2,CS_max,LB_river,RDS_mainDist,STS,prec_mm,pet_mm,dist_to_LB_km,m_per_year,TB)%>%
  rename(Country=CC,nWells=n,beta_AveDepletion=beta_0, se_AveDepletion=se_0, p_AveDepletion=p_0,beta_distTrend=beta_1,se_distTrend=se_1,p_distTrend=p_1,CropSuitIndex=CS_max,probBorderRiver=LB_river,SoitTerrainSuit=STS,depletion_m_yr=m_per_year,TransboundaryAquifer=TB)

write.csv(xx,'SI_Data_Wells.csv')


# ============================================================
# Helper functions
# ============================================================
classify_slope <- function(beta, p) {
  case_when(
    beta > 0 & p < 0.1  ~ "pos_sig",
    beta > 0 & p >= 0.1 ~ "pos_nonsig",
    beta < 0 & p >= 0.1 ~ "neg_nonsig",
    beta < 0 & p < 0.1  ~ "neg_sig",
    TRUE ~ NA_character_
  )
}

weighted_mean <- function(x, w) weighted.mean(x, w, na.rm = TRUE)
weighted_se <- function(x, w) {
  m <- weighted_mean(x, w)
  sqrt(sum(w * (x - m)^2, na.rm = TRUE) / (sum(w, na.rm = TRUE)^2))
}

# ============================================================
# Summary statistics
# ============================================================
df <- read.csv('../2_firststage/firstStageMain.csv')%>%
  mutate(
    class_1 = classify_slope(beta_1, p_1),
    class_0  = classify_slope(beta_1, p_1)
  )



summarize_vars <- c(
  "CS_max", "area_km2", "urbkHaKm2", "RDS_mainDist",'prec_mm','pet_mm','STS','LB_river',
  "m_per_year", "dist_to_LB_km", "nW"
)

summary_stats <- function(data, w = NULL) {
  library(dplyr)
  library(Hmisc)  # for wtd.quantile
  
  data %>%
    summarise(across(all_of(summarize_vars), ~ {
      if (is.null(w)) {
        median_val <- median(.x, na.rm = TRUE)
        q <- quantile(.x, probs = c(0.25, 0.75), na.rm = TRUE)
      } else {
        median_val <- Hmisc::wtd.quantile(.x, weights = w, probs = 0.5, na.rm = TRUE)
        q <- Hmisc::wtd.quantile(.x, weights = w, probs = c(0.25, 0.75), na.rm = TRUE)
      }
      sprintf("%.2f [%.2f, %.2f]", median_val, q[1], q[2])
    }))
}
sum_df <- list(
  df %>% filter(!TB) %>% summary_stats() %>% mutate(TB = "Non-Transboundary"),
  df %>% filter(!TB) %>% summary_stats(w = df$weights[df$TB == FALSE]) %>% mutate(TB = "Matched non-TB"),
  df %>% filter(TB) %>% summary_stats() %>% mutate(TB = "Transboundary")
) %>%
  bind_rows() %>%
  dplyr::select(TB, everything())

# Transpose for display
tab <- sum_df %>%
  column_to_rownames("TB") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable")

colnames(tab) <- c("Variable", "Non-Transboundary", "Matched non-TB", "Transboundary")

kable(tab,
      format = "latex",
      booktabs = TRUE,
      caption = "Summary statistics by transboundary status, including matched non-TB sample: Median and [Quartiles]",
      align = "lccc") %>%
  kable_styling(latex_options = c("hold_position"))


