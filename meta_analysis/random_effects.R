# ============================================================================
# RANDOM EFFECTS META-ANALYSIS
# Biostatistical Analysis Platform for Medical Research
# ============================================================================

#' @title Random Effects Meta-Analysis
#' @description Professional implementation of random effects meta-analysis
#'   using 'meta' and 'metafor' packages for clinical research
#' @author Gabriel Demetrios Lafis
#' @date September 2025

# Required Libraries ----------------------------------------------------------
library(meta)
library(metafor)
library(ggplot2)
library(dplyr)

# Main Random Effects Function ------------------------------------------------

#' Perform Random Effects Meta-Analysis (REML default)
#'
#' @param data Data frame with study data
#' @param effect_col Column name containing effect estimates
#' @param variance_col Column name containing variances  
#' @param study_col Column name containing study identifiers
#' @param tau2_method Method for tau^2 estimation ("REML", "DL", "SJ", "HE")
#' @param alpha Significance level (default 0.05)
#' @param digits Number of decimal places for output
#'
#' @return List containing meta-analysis results
#' @export
random_effects_meta <- function(data,
                                effect_col,
                                variance_col,
                                study_col,
                                tau2_method = "REML",
                                alpha = 0.05,
                                digits = 3) {
  
  # Validate input
  stopifnot(is.data.frame(data))
  required_cols <- c(effect_col, variance_col, study_col)
  if (!all(required_cols %in% names(data))) {
    stop("Required columns not found in data")
  }
  
  # Extract vectors
  effects <- data[[effect_col]]
  variances <- data[[variance_col]]
  studies <- data[[study_col]]
  
  # Fit random effects model using metafor
  re_model <- metafor::rma(yi = effects,
                           vi = variances,
                           slab = studies,
                           method = tau2_method,
                           level = (1 - alpha) * 100)
  
  # Heterogeneity and influence diagnostics
  inf <- metafor::influence(re_model)
  baujat <- metafor::baujat(re_model, plot = FALSE)
  
  # Prediction interval
  pred_int <- metafor::predict(re_model, transf = NULL)
  
  results <- list(
    method = paste0("Random Effects (", tau2_method, ")"),
    model = re_model,
    pooled_effect = round(re_model$beta, digits),
    pooled_se = round(re_model$se, digits),
    pooled_ci_lower = round(re_model$ci.lb, digits),
    pooled_ci_upper = round(re_model$ci.ub, digits),
    prediction_interval = c(lower = round(pred_int$cr.lb, digits),
                            upper = round(pred_int$cr.ub, digits)),
    z_value = round(re_model$zval, digits),
    p_value = round(re_model$pval, digits),
    heterogeneity = list(
      Q = round(re_model$QE, digits),
      df = re_model$k - 1,
      p_het = round(re_model$QEp, digits),
      I2 = round(re_model$I2, 1),
      H2 = round(re_model$H2, 3),
      tau2 = round(re_model$tau2, digits)
    ),
    influence = inf,
    baujat = baujat,
    studies = studies,
    n_studies = re_model$k
  )
  
  return(results)
}

# Forest Plot -----------------------------------------------------------------

#' Forest plot for Random Effects Meta-Analysis
#' @param meta_results Results from random_effects_meta
#' @param title Title
#' @param xlabel X label
#' @param show_pi Show prediction interval band
#' @return base plot created by metafor
#' @export
plot_random_effects_forest <- function(meta_results,
                                       title = "Random Effects Meta-Analysis",
                                       xlabel = "Effect Size",
                                       show_pi = TRUE) {
  
  p <- metafor::forest(meta_results$model,
                       slab = meta_results$studies,
                       header = TRUE,
                       xlab = xlabel,
                       main = title)
  
  if (show_pi && !is.null(meta_results$prediction_interval)) {
    abline(v = meta_results$prediction_interval, col = "steelblue", lty = 3)
  }
  
  return(p)
}

# Influence and Baujat Plots --------------------------------------------------

#' Influence plot for Random Effects model
#' @export
plot_influence_random <- function(meta_results) {
  metafor::influence(meta_results$model, plot = TRUE)
}

#' Baujat plot to identify studies contributing to heterogeneity
#' @export
plot_baujat_random <- function(meta_results) {
  metafor::baujat(meta_results$model)
}

# Binary Outcomes --------------------------------------------------------------

#' Random Effects Meta-Analysis for Binary Outcomes
#' @export
random_effects_binary <- function(data,
                                  events_exp,
                                  n_exp,
                                  events_ctrl,
                                  n_ctrl,
                                  study_col,
                                  measure = "OR",
                                  method.tau = "REML") {
  meta_bin <- meta::metabin(event.e = data[[events_exp]],
                            n.e = data[[n_exp]],
                            event.c = data[[events_ctrl]],
                            n.c = data[[n_ctrl]],
                            studlab = data[[study_col]],
                            sm = measure,
                            method.tau = method.tau,
                            comb.fixed = FALSE,
                            comb.random = TRUE)
  return(meta_bin)
}

# Continuous Outcomes ---------------------------------------------------------

#' Random Effects Meta-Analysis for Continuous Outcomes
#' @export
random_effects_continuous <- function(data,
                                      mean_exp, sd_exp, n_exp,
                                      mean_ctrl, sd_ctrl, n_ctrl,
                                      study_col,
                                      measure = "SMD",
                                      method.tau = "REML") {
  meta_cont <- meta::metacont(n.e = data[[n_exp]],
                              mean.e = data[[mean_exp]],
                              sd.e = data[[sd_exp]],
                              n.c = data[[n_ctrl]],
                              mean.c = data[[mean_ctrl]],
                              sd.c = data[[sd_ctrl]],
                              studlab = data[[study_col]],
                              sm = measure,
                              method.tau = method.tau,
                              comb.fixed = FALSE,
                              comb.random = TRUE)
  return(meta_cont)
}

# Export ----------------------------------------------------------------------

#' Export Random Effects Results
#' @export
export_random_effects <- function(meta_results,
                                  filename,
                                  format = "pdf",
                                  width = 10,
                                  height = 8) {
  if (format %in% c("pdf", "png")) {
    if (format == "pdf") pdf(paste0(filename, ".pdf"), width = width, height = height)
    if (format == "png") png(paste0(filename, ".png"), width = width*100, height = height*100, res = 300)
    plot_random_effects_forest(meta_results)
    dev.off()
  } else if (format == "csv") {
    df <- data.frame(
      pooled = meta_results$pooled_effect,
      se = meta_results$pooled_se,
      ci_lb = meta_results$pooled_ci_lower,
      ci_ub = meta_results$pooled_ci_upper,
      pi_lb = meta_results$prediction_interval[1],
      pi_ub = meta_results$prediction_interval[2],
      Q = meta_results$heterogeneity$Q,
      I2 = meta_results$heterogeneity$I2,
      tau2 = meta_results$heterogeneity$tau2
    )
    write.csv(df, paste0(filename, ".csv"), row.names = FALSE)
  }
  cat("Random effects results exported to:", paste0(filename, ".", format), "\n")
}

# Example ---------------------------------------------------------------------

#' Example: Random Effects Meta-Analysis
example_random_effects <- function() {
  set.seed(123)
  example_data <- data.frame(
    study = paste("Study", 1:12),
    effect = rnorm(12, 0.4, 0.25),
    variance = runif(12, 0.02, 0.08)
  )
  results <- random_effects_meta(
    data = example_data,
    effect_col = "effect",
    variance_col = "variance",
    study_col = "study",
    tau2_method = "REML"
  )
  summarize <- c(
    results$pooled_effect, results$pooled_ci_lower, results$pooled_ci_upper,
    results$heterogeneity$I2, results$heterogeneity$tau2
  )
  print(summarize)
  plot_random_effects_forest(results)
  return(results)
}

# ============================================================================
# END OF RANDOM EFFECTS META-ANALYSIS MODULE
# ============================================================================
