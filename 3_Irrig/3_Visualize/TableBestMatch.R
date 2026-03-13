getModel <- function(res, rank = 1) {
  # light objects: use top_models
  if ("top_models" %in% names(res)) {
    if (rank > length(res$top_models)) {
      stop("Requested rank > number of stored top_models.")
    }
    return(res$top_models[[rank]])
  }
  
  # fallback for old/full res objects (if you still have any)
  if ("results" %in% names(res) && "summary_df" %in% names(res)) {
    idx <- res$summary_df$idx[rank]
    return(res$results[[idx]]$fit)
  }
  
  stop("Object does not look like a pilot_match result.")
}

generate_res_table <- function(
    res_list,
    rank = 1,
    output = "latex",        # explicitly request latex
    stars = TRUE,
    labels = NULL,
    caption = NULL
) {
  library(modelsummary)
  library(purrr)
  
  # Extract model fits
  models <- map(res_list, getModel, rank = rank)
  
  # Column labels
  if (is.null(labels)) {
    names(models) <- paste0("Model_", seq_along(models))
  } else {
    if (length(labels) != length(models)) {
      stop("Length of 'labels' must match number of models in 'res_list'.")
    }
    names(models) <- labels
  }
  
  # --- LaTeX options: disable talltblr ---
  options(modelsummary_factory_latex = "kableExtra")  # forces tabular-based output
  
  # Build regression table
  modelsummary(
    models,
    stars    = stars,
    gof_omit = "AIC|BIC|Log.Lik|Deviance",
    title    = caption,
    output   = output,
    add_rows = NULL
  )
}


# source('Fct.R')
library(qs)
#load existing analyses
res_files <- list.files("resultOut", pattern = "\\.qs$", full.names = TRUE)
res_objects <- lapply(res_files, qs::qread)
names(res_objects) <- gsub("\\.qs$", "", basename(res_files))
list2env(res_objects, envir = .GlobalEnv)

#Main fig.
IrIntens$label="Irrigation"
Overdraft$label="Overdraft"
IrRivsInt$label="Irrig | Rivers"
IrGini$label='Irrigation'
IrRivsGini$label="Irrig | Rivers"
GWGini$label="GW Irrig."
GWRivsGini$label="GW Irrig. | Rivers"

x=list(IrIntens,Overdraft,IrRivsInt,IrGini,IrRivsGini,GWGini,GWRivsGini)
generate_res_table(
  x,
  rank = 1,
  output = "latex",        # explicitly request latex
  stars = TRUE,
  labels =  sapply(x, `[[`, "label"),
  caption = NULL
)
