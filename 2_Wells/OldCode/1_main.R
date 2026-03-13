source('_Fct2.R')

#MAIN, also generratess plots, including the bar plot to be included instead of the bins. 
a2=run(plot=T)[[2]]
a2$lab='Preferred'
regressplot(a2)


##Regression table with robustness tests
#Unmatched
a1=run(matchWg=F)[[2]]
a1$lab='UnMatched'


#Robust lm #barely ns.
a3=run(rob=T)[[2]]
a3$lab='Robust LM'
#alt covs
a4=run(cov_match =c("lat_c","lon_c",'RDS_mainDist','STS','prec_mm','pet_mm','LB_river'))[[2]]
a4$lab='Alt Cov'
# #double robust# barely ns
a5=run(tVar=c('urbkHaKm2','CS_max','LB_river'),cov_match =c("lat_c","lon_c",'urbkHaKm2','CS_max','LB_river'))[[2]]
a5$lab='Double Robust'

a6=run(dropTopSE=T)[[2]]
a6$lab='Drop Heavy'

# dclust = NULL,         # NEW: radius in km, or NULL
# kde = FALSE,           # NEW: if TRUE use KDE weights; default fixed-radius
# lon_col = "lon",       # NEW: coords for weights

# a7=run(dclust=50)[[2]]
# a7$lab='Spatial Clust'

mods=list(a1,a2,a3,a4,a6)
make_regtable(mods,0,'Effect of TB status on mean depletion. (1) No matching weight; (2) Preferred specification: full matching weight, with LatLon, urban cover, Crop Suitability and Propensity of river border as matching covariates; (3) Aquifer-level robust linear models; (4) Alternative matching covariates:  Soil and terrain suitability, precipitations, potential ET; (5) Double Robust specification where matching covariates are also used as aquifer-level covariates; (6) Exclusion of 10 control aquifers with highest matching weights')
make_regtable(mods,1,'Effect of TB status on distance to border trend. (1) No matching weight; (2) Preferred specification: full matching weight, with LatLon, urban cover, Crop Suitability and Propensity of river border as matching covariates; (3) Aquifer-level robust linear models; (4) Alternative matching covariates:  Soil and terrain suitability, precipitations, potential ET; (5) Double Robust specification where matching covariates are also used as aquifer-level covariates; (6) Exclusion of 10 control aquifers with highest matching weights')


#####################################################################################################################
#######################################Robustness to alternate matching architecture#######################################
#####################################################################################################################
#nearest match, with varrying ratios: 
r_seq <- 1:10
results_ratio <- map_dfr(r_seq, function(r) {
  cat("Running ratio =", r, "\n")
  res <- run(
    match = 'nearest',   # nearest neighbor with ratio control
    r = r
  )[[2]]
  res$ratio <- r
  res
})
future::plan(sequential)
save(results_ratio, file = "results_ratio_sensitivity_nearest.rdata")
load("results_ratio_sensitivity_nearest.rdata")
plot_df <- bind_rows(
  # Regression 0: "Mean depletion"
  results_ratio %>%
    transmute(
      ratio,
      term  = nm_0,
      model = "Mean depletion",
      est   = b_0,
      se    = se_0
    ),
  # Regression 1: "Distance Trend"
  results_ratio %>%
    transmute(
      ratio,
      term  = nm_1,
      model = "Distance Trend",
      est   = b_1,
      se    = se_1
    )
) %>%
  filter(term %in% c("intrcpt", "TBTRUE")) %>%
  mutate(
    term   = recode(term, intrcpt = "(Intercept)", TBTRUE = "TBTRUE"),
    ci_low = est - 1.64 * se,
    ci_high = est + 1.64 * se
  )%>%mutate(model=factor(model,levels=c('Mean depletion', 'Distance Trend')))

ggplot(plot_df, aes(x = ratio, y = est, color = term)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = term),
              alpha = 0.25, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ model, scales = "free_y") +
  scale_color_manual(values = c("(Intercept)" = "darkorange", "TBTRUE" = "#1f78b4")) +
  scale_fill_manual(values = c("(Intercept)" = "darkorange", "TBTRUE" = "#1f78b4")) +
  labs(
    x = "Matching ratio (controls per treated)",
    y = "Estimated effect"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank())


#####################################################################################################################
#######################################Robustness to data cutoff#######################################
#####################################################################################################################
# nMin_seq <- c(10,20,30,40,50,60,70,80,90,100)
# results_nMin <- map_dfr(nMin_seq, function(nm) {
#   cat("Running nMin =", nm, "\n")
#   res <- run(
#     nMin = nm
#   )[[2]]  # extract your second element
#   res$nMin <- nm
#   res
# })
# save(results_nMin,file='results_nMin.rdata')

# Load results
load("results_nMin.rdata")

library(dplyr)
library(ggplot2)
library(patchwork)

# --- Prepare data from numeric tibble ---
plot_df <- results_nMin %>%
  transmute(
    nMin,
    term = nm_0,
    model = "Mean depletion",
    est = b_0,
    se = se_0
  ) %>%
  bind_rows(
    results_nMin %>%
      transmute(
        nMin,
        term = nm_1,
        model = "Distance Trend",
        est = b_1,
        se = se_1
      )
  ) %>%
  mutate(
    ci_low = est - 1.64 * se,
    ci_high = est + 1.64 * se
  ) %>%
  filter(term %in% c("intrcpt", "TBTRUE")) %>%
  mutate(
    term = recode(term, intrcpt = "(Intercept)", TBTRUE = "TBTRUE")
  )

# --- y-axis limits (adjustable) ---
ylims <- data.frame(
  model = c("Distance Trend", "Mean depletion"),
  ymin = c(-1, -220),
  ymax = c(1, 200)
)

# --- Plotting function ---
make_plot <- function(model_name) {
  lims <- ylims %>% filter(model == model_name)
  
  ggplot(filter(plot_df, model == model_name),
         aes(x = nMin, y = est, color = term, fill = term)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.25, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("(Intercept)" = "darkorange", "TBTRUE" = "#1f78b4")) +
    scale_fill_manual(values = c("(Intercept)" = "darkorange", "TBTRUE" = "#1f78b4")) +
    coord_cartesian(ylim = c(lims$ymin, lims$ymax)) +
    labs(
      x = "Minimum Obs. cutoff (nMin)",
      y = "Estimated effect",
      title = model_name
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.title = element_blank())
}

# --- Generate plots and combine ---
p1 <- make_plot("Mean depletion") + theme(legend.position = "none")
p2 <- make_plot("Distance Trend")

p1 + p2

#####################################################################################################################
#######################################Leave one out robustness check#######################################
#####################################################################################################################
# 
# # #leave x out robustness to evaluate external validity. There are 40 tbas, so this is not really flexible
# library(dplyr)
# library(progressr)
# library(furrr)
# library(progressr)
# library(purrr)
# library(dplyr)
# library(progressr)
# library(purrr)
# library(dplyr)
# 
# x <- run()
# 
# # Vector of aquifers to sample from
# TBA <- x[[1]] %>% filter(TB) %>% pull(Aquifer)
# 
# # Parameters
# M <- 250                               # Number of repetitions
# k_values <- c(1, 5, 10, 15, 20, 25, 30, 35, 40)
# 
# handlers(global = TRUE)
# handlers("progress")
# 
# resultsLOO <- with_progress({
#   total_steps <- M * length(k_values)
#   p <- progressor(steps = total_steps)
#   
#   results_all <- vector("list", M)
#   
#   for (m in seq_len(M)) {
#     results_m <- vector("list", length(k_values))
#     
#     for (ki in seq_along(k_values)) {
#       k <- k_values[ki]
#       p(message = sprintf("Iteration m=%d, dropping k=%d aquifers", m, k))
#       
#       set.seed(1000 + k * 100 + m)
#       TBA_sel <- sample(TBA, k)
#       
#       res <- tryCatch(run(loo = TBA_sel)[[2]], error = function(e) NULL)
#       if (!is.null(res)) {
#         res$looaqf <- paste(TBA_sel, collapse = "; ")
#         res$k <- k
#         res$iter <- m
#       }
#       
#       results_m[[ki]] <- res
#     }
#     
#     # Combine all results from this repetition
#     results_all[[m]] <- bind_rows(results_m)
#     
#     # Temporary save after each repetition
#     save(results_all, file = 'resultsLOO_final.rdata')
#   }
#   
#   bind_rows(results_all)
# })
# 
# # Final save
# save(resultsLOO, file = "resultsLOO_final.rdata")

#remove duplicate in romoved aquifers
load("resultsLOO_final.rdata")
p1=plot_loo_robustness(resultsLOO, reg = 1, coef = "treatment")
p2=plot_loo_robustness(resultsLOO, reg = 1, coef = "intercept")
p3=plot_loo_robustness(resultsLOO, reg = 0, coef = "treatment")
p3+p1+p2

##Robustness to clustering dist
dclust_seq <- c(10,20,30,40,50,60,70,80,90,100)
results_dclust <- map_dfr(dclust_seq, function(nm) {
  cat("Running nMin =", nm, "\n")
  res <- run(
    dclust = nm,kde=T
  )[[2]]  # extract your second element
  res$dclust <- nm
  res
})
results_dclust=rbind(results_dclust,data.frame(run(dclust=NULL)[[2]],dclust=0))

save(results_dclust,file='results_dclust.rdata')

# Load results
load("results_dclust.rdata")

library(dplyr)
library(ggplot2)
library(patchwork)

# --- Prepare data from numeric tibble ---
plot_df <- results_dclust %>%
  transmute(
    dclust,
    term = nm_0,
    model = "Mean depletion",
    est = b_0,
    se = se_0
  ) %>%
  bind_rows(
    results_dclust %>%
      transmute(
        dclust,
        term = nm_1,
        model = "Distance Trend",
        est = b_1,
        se = se_1
      )
  ) %>%
  mutate(
    ci_low = est - 1.64 * se,
    ci_high = est + 1.64 * se
  ) %>%
  filter(term %in% c("intrcpt", "TBTRUE")) %>%
  mutate(
    term = recode(term, intrcpt = "(Intercept)", TBTRUE = "TBTRUE")
  )

# --- y-axis limits (adjustable) ---
ylims <- data.frame(
  model = c("Distance Trend", "Mean depletion"),
  ymin = c(-1, -220),
  ymax = c(1, 200)
)

# --- Plotting function ---
make_plot <- function(model_name) {
  lims <- ylims %>% filter(model == model_name)
  
  ggplot(filter(plot_df, model == model_name),
         aes(x = dclust, y = est, color = term, fill = term)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.25, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("(Intercept)" = "darkorange", "TBTRUE" = "#1f78b4")) +
    scale_fill_manual(values = c("(Intercept)" = "darkorange", "TBTRUE" = "#1f78b4")) +
    coord_cartesian(ylim = c(lims$ymin, lims$ymax)) +
    labs(
      x = "declustering distance (km)",
      y = "Estimated effect",
      title = model_name
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.title = element_blank())
}

p1 <- make_plot("Mean depletion") + theme(legend.position = "none")
p2 <- make_plot("Distance Trend")

p1 + p2
