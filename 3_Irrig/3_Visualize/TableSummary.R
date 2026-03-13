
library(dplyr)
library(knitr)
library(kableExtra)
generate_summary_table <- function(dout) {
  library(dplyr)
  library(knitr)
  library(kableExtra)
  
  # --- Continuous vars ---
  cont_vars <- c("CSI","RDS_mainDist", "area_km2", "GWIrrig", "Irrig", "IrrigNeed", "Overdraft")
  cont_summary <- lapply(cont_vars, function(v) {
    x <- dout[[v]] * 100
    zeros <- sum(x == 0, na.rm = TRUE)
    x_pos <- x[x > 0 & !is.na(x)]
    tibble(
      Variable = v,
      Zeros = zeros,
      Q25 = if (length(x_pos)) quantile(x_pos, 0.25, na.rm = TRUE) else NA,
      Median = if (length(x_pos)) median(x_pos, na.rm = TRUE) else NA,
      Q75 = if (length(x_pos)) quantile(x_pos, 0.75, na.rm = TRUE) else NA
    )
  }) %>% bind_rows()
  
  # --- Gini vars ---
  gini_vars <- c("G_Irrig", "G_IrrigGW", "G_IrrigNeed", "G_CSI")
  gini_summary <- lapply(gini_vars, function(v) {
    x <- dout[[v]]
    tibble(
      Variable = v,
      NAs = sum(is.na(x)),
      Q25 = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Q75 = quantile(x, 0.75, na.rm = TRUE)
    )
  }) %>% bind_rows()
  
  # --- Categorical vars ---
  river_tab <- dout %>%
    count(RiverBorder) %>%
    mutate(Variable = paste0("RiverBorder=", RiverBorder)) %>%
    select(Variable, Count = n)
  
  region_tab <- dout %>%
    count(Region) %>%
    mutate(Variable = paste0("Region=", Region)) %>%
    select(Variable, Count = n)
  
  cat_summary <- tibble(
    Variable = c("AquiferName (unique categories)", "CountryName (unique categories)"),
    Count = c(n_distinct(dout$AquiferName), n_distinct(dout$CountryName))
  )
  
  total_row <- tibble(Variable = "Total rows", Count = nrow(dout))
  
  # --- Combine ---
  final_table <- bind_rows(cont_summary, gini_summary, river_tab, region_tab, cat_summary, total_row)
  
  # --- Format numeric values ---
  format_num <- function(x) {
    ifelse(is.na(x) | x == "", "",
           ifelse(abs(x) < 0.01 & abs(x) > 0, formatC(x, format = "e", digits = 1),
                  formatC(x, format = "f", digits = 2)))
  }
  
  final_table <- final_table %>%
    mutate(across(where(is.numeric), format_num))
  
  # --- LaTeX table ---
  kbl(
    final_table,
    format = "latex",
    booktabs = TRUE,
    caption = "Summary statistics for aquifer-level dataset",
    align = "lcccccc"
  ) %>%
    kable_styling(latex_options = c("hold_position"), font_size = 9) %>%
    pack_rows("Continuous variables", 1, length(cont_vars)) %>%
    pack_rows("Gini indices", length(cont_vars) + 1, length(cont_vars) + length(gini_vars)) %>%
    pack_rows("Categorical variables", length(cont_vars) + length(gini_vars) + 1,
              nrow(final_table) - 1)
}

dout=read.csv('../1_buildData/Data_S2.csv')
generate_summary_table(dout)

