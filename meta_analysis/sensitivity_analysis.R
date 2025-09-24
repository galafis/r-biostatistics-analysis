# ============================================================================
# SENSITIVITY ANALYSIS FOR META-ANALYSIS
# Biostatistical Analysis Platform for Medical Research
# ============================================================================

#' @title Sensitivity Analyses for Meta-Analyses
#' @description Implements leave-one-out, subgroup, meta-regression, and
#'   influence diagnostics using metafor/meta.
#' @author Gabriel Demetrios Lafis
#' @date September 2025

library(metafor)
library(meta)
library(dplyr)

# Leave-one-out ---------------------------------------------------------------

#' Leave-one-out sensitivity for metafor rma model
#' @export
leave1out_meta <- function(rma_model) {
  tryCatch(metafor::leave1out(rma_model), error = function(e) e)
}

# Subgroup analysis -----------------------------------------------------------

#' Subgroup analysis for binary outcomes (meta)
#' @export
subgroup_binary <- function(data,
                            events_exp, n_exp,
                            events_ctrl, n_ctrl,
                            study_col, subgroup_col,
                            measure = "OR",
                            method = "MH") {
  m <- meta::metabin(event.e = data[[events_exp]], n.e = data[[n_exp]],
                     event.c = data[[events_ctrl]], n.c = data[[n_ctrl]],
                     studlab = data[[study_col]], sm = measure,
                     byvar = data[[subgroup_col]], method = method,
                     comb.fixed = TRUE, comb.random = TRUE)
  return(m)
}

# Meta-regression -------------------------------------------------------------

#' Meta-regression using metafor rma
#' @param yi Effect sizes
#' @param vi Variances
#' @param moderators Data frame or matrix of moderators
#' @export
meta_regression <- function(yi, vi, moderators, method = "REML") {
  df <- as.data.frame(moderators)
  form <- as.formula(paste("~", paste(names(df), collapse = "+")))
  metafor::rma(yi = yi, vi = vi, mods = form, method = method)
}

# Influence diagnostics -------------------------------------------------------

#' Influence diagnostics for metafor rma
#' @export
influence_diagnostics <- function(rma_model) {
  list(influence = metafor::influence(rma_model),
       hatvalues = hatvalues(rma_model))
}

# Export utilities ------------------------------------------------------------

#' Export leave-one-out and influence plots
#' @export
export_sensitivity <- function(rma_model, filename_prefix) {
  loo <- tryCatch(leave1out_meta(rma_model), error = function(e) NULL)
  if (!is.null(loo)) {
    pdf(paste0(filename_prefix, "_leave1out.pdf"), width = 9, height = 6)
    plot(loo)
    dev.off()
  }
  
  pdf(paste0(filename_prefix, "_influence.pdf"), width = 9, height = 6)
  plot(influence(rma_model))
  dev.off()
}

# Example ---------------------------------------------------------------------

#' Example: Sensitivity analyses
example_sensitivity <- function() {
  set.seed(123)
  yi <- rnorm(20, 0.3, 0.25)
  vi <- runif(20, 0.02, 0.06)
  m <- rma(yi = yi, vi = vi, method = "REML")
  print(leave1out_meta(m))
  print(influence_diagnostics(m))
  export_sensitivity(m, "sensitivity_demo")
  invisible(m)
}

# ============================================================================
# END OF SENSITIVITY ANALYSIS MODULE
# ============================================================================
