# ============================================================================
# FIXED EFFECTS META-ANALYSIS
# Biostatistical Analysis Platform for Medical Research
# ============================================================================

#' @title Fixed Effects Meta-Analysis
#' @description Professional implementation of fixed effects meta-analysis
#'   using 'meta' and 'metafor' packages for clinical research
#' @author Gabriel Demetrios Lafis
#' @date September 2025

# Required Libraries ----------------------------------------------------------
library(meta)
library(metafor)
library(ggplot2)
library(dplyr)
library(forestplot)

# Main Fixed Effects Function -------------------------------------------------

#' Perform Fixed Effects Meta-Analysis
#'
#' @param data Data frame with study data
#' @param effect_col Column name containing effect estimates
#' @param variance_col Column name containing variances  
#' @param study_col Column name containing study identifiers
#' @param method Method for effect size calculation ("SMD", "MD", "OR", "RR")
#' @param alpha Significance level (default 0.05)
#' @param digits Number of decimal places for output
#'
#' @return List containing meta-analysis results
#' @export
fixed_effects_meta <- function(data, 
                               effect_col, 
                               variance_col, 
                               study_col,
                               method = "SMD",
                               alpha = 0.05,
                               digits = 3) {
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("Data must be a data frame")
  }
  
  required_cols <- c(effect_col, variance_col, study_col)
  if (!all(required_cols %in% names(data))) {
    stop("Required columns not found in data")
  }
  
  # Extract vectors
  effects <- data[[effect_col]]
  variances <- data[[variance_col]]
  studies <- data[[study_col]]
  
  # Perform fixed effects meta-analysis using metafor
  fe_model <- metafor::rma(yi = effects,
                          vi = variances,
                          slab = studies,
                          method = "FE",
                          level = (1 - alpha) * 100)
  
  # Calculate weights
  weights <- 1 / variances
  weights_percent <- round((weights / sum(weights)) * 100, 2)
  
  # Create results summary
  results <- list(
    method = "Fixed Effects",
    model = fe_model,
    pooled_effect = round(fe_model$beta, digits),
    pooled_se = round(fe_model$se, digits),
    pooled_ci_lower = round(fe_model$ci.lb, digits),
    pooled_ci_upper = round(fe_model$ci.ub, digits),
    z_value = round(fe_model$zval, digits),
    p_value = round(fe_model$pval, digits),
    heterogeneity = list(
      Q = round(fe_model$QE, digits),
      df = fe_model$k - 1,
      p_het = round(fe_model$QEp, digits),
      I2 = ifelse(is.na(fe_model$I2), 0, round(fe_model$I2, 1))
    ),
    weights = weights_percent,
    studies = studies,
    n_studies = fe_model$k
  )
  
  return(results)
}

# Forest Plot Function --------------------------------------------------------

#' Create Professional Forest Plot for Fixed Effects Meta-Analysis
#'
#' @param meta_results Results from fixed_effects_meta function
#' @param title Plot title
#' @param xlabel X-axis label
#' @param show_weights Show study weights (default TRUE)
#' @param text_size Text size for plot elements
#'
#' @return ggplot object
#' @export
plot_fixed_effects_forest <- function(meta_results,
                                       title = "Fixed Effects Meta-Analysis",
                                       xlabel = "Effect Size",
                                       show_weights = TRUE,
                                       text_size = 12) {
  
  # Generate forest plot using metafor
  forest_plot <- metafor::forest(meta_results$model,
                                  slab = meta_results$studies,
                                  showweights = show_weights,
                                  header = TRUE,
                                  main = title,
                                  xlab = xlabel,
                                  cex = text_size/12,
                                  col = "blue")
  
  return(forest_plot)
}

# Summary Statistics Function -------------------------------------------------

#' Generate Summary Statistics Table
#'
#' @param meta_results Results from fixed_effects_meta function
#' @param export_format Export format ("console", "latex", "html")
#'
#' @return Formatted summary table
#' @export
summarize_fixed_effects <- function(meta_results, export_format = "console") {
  
  # Create summary table
  summary_table <- data.frame(
    Statistic = c("Pooled Effect", "Standard Error", "95% CI Lower", 
                  "95% CI Upper", "Z-value", "P-value", "Q-statistic",
                  "P-heterogeneity", "I-squared", "Number of Studies"),
    Value = c(
      meta_results$pooled_effect,
      meta_results$pooled_se,
      meta_results$pooled_ci_lower,
      meta_results$pooled_ci_upper,
      meta_results$z_value,
      meta_results$p_value,
      meta_results$heterogeneity$Q,
      meta_results$heterogeneity$p_het,
      paste0(meta_results$heterogeneity$I2, "%"),
      meta_results$n_studies
    )
  )
  
  # Format output based on export format
  if (export_format == "console") {
    print(summary_table)
  } else if (export_format == "latex") {
    return(knitr::kable(summary_table, format = "latex"))
  } else if (export_format == "html") {
    return(knitr::kable(summary_table, format = "html"))
  }
  
  return(summary_table)
}

# Binary Outcomes Function ---------------------------------------------------

#' Fixed Effects Meta-Analysis for Binary Outcomes
#'
#' @param data Data frame with study data
#' @param events_exp Events in experimental group
#' @param n_exp Sample size in experimental group
#' @param events_ctrl Events in control group
#' @param n_ctrl Sample size in control group
#' @param study_col Study identifiers
#' @param measure Effect measure ("OR", "RR", "RD")
#' @param method Calculation method ("MH", "Peto", "Inverse")
#'
#' @return Meta-analysis results for binary outcomes
#' @export
fixed_effects_binary <- function(data,
                                 events_exp,
                                 n_exp,
                                 events_ctrl,
                                 n_ctrl,
                                 study_col,
                                 measure = "OR",
                                 method = "MH") {
  
  # Perform meta-analysis using meta package
  meta_bin <- meta::metabin(event.e = data[[events_exp]],
                           n.e = data[[n_exp]],
                           event.c = data[[events_ctrl]],
                           n.c = data[[n_ctrl]],
                           studlab = data[[study_col]],
                           sm = measure,
                           method = method,
                           comb.fixed = TRUE,
                           comb.random = FALSE)
  
  return(meta_bin)
}

# Continuous Outcomes Function -----------------------------------------------

#' Fixed Effects Meta-Analysis for Continuous Outcomes
#'
#' @param data Data frame with study data
#' @param mean_exp Mean in experimental group
#' @param sd_exp Standard deviation in experimental group
#' @param n_exp Sample size in experimental group
#' @param mean_ctrl Mean in control group
#' @param sd_ctrl Standard deviation in control group
#' @param n_ctrl Sample size in control group
#' @param study_col Study identifiers
#' @param measure Effect measure ("SMD", "MD")
#'
#' @return Meta-analysis results for continuous outcomes
#' @export
fixed_effects_continuous <- function(data,
                                     mean_exp, sd_exp, n_exp,
                                     mean_ctrl, sd_ctrl, n_ctrl,
                                     study_col,
                                     measure = "SMD") {
  
  # Perform meta-analysis using meta package
  meta_cont <- meta::metacont(n.e = data[[n_exp]],
                             mean.e = data[[mean_exp]],
                             sd.e = data[[sd_exp]],
                             n.c = data[[n_ctrl]],
                             mean.c = data[[mean_ctrl]],
                             sd.c = data[[sd_ctrl]],
                             studlab = data[[study_col]],
                             sm = measure,
                             comb.fixed = TRUE,
                             comb.random = FALSE)
  
  return(meta_cont)
}

# Export Functions ------------------------------------------------------------

#' Export Meta-Analysis Results
#'
#' @param meta_results Results from meta-analysis
#' @param filename Output filename
#' @param format Export format ("pdf", "png", "csv")
#' @param width Plot width (inches)
#' @param height Plot height (inches)
#'
#' @export
export_fixed_effects <- function(meta_results, 
                                  filename,
                                  format = "pdf",
                                  width = 10,
                                  height = 8) {
  
  if (format %in% c("pdf", "png")) {
    # Export forest plot
    if (format == "pdf") {
      pdf(paste0(filename, ".pdf"), width = width, height = height)
    } else {
      png(paste0(filename, ".png"), width = width*100, height = height*100, res = 300)
    }
    
    plot_fixed_effects_forest(meta_results)
    dev.off()
    
  } else if (format == "csv") {
    # Export summary table
    summary_df <- summarize_fixed_effects(meta_results, "console")
    write.csv(summary_df, paste0(filename, ".csv"), row.names = FALSE)
  }
  
  cat("Results exported to:", paste0(filename, ".", format), "\n")
}

# Example Usage ---------------------------------------------------------------

#' Example: Fixed Effects Meta-Analysis
#'
#' This example demonstrates how to perform a fixed effects meta-analysis
#' using simulated clinical trial data
example_fixed_effects <- function() {
  
  # Simulate example data
  set.seed(123)
  example_data <- data.frame(
    study = paste("Study", 1:10),
    effect = rnorm(10, 0.5, 0.2),
    variance = runif(10, 0.01, 0.05),
    sample_size = sample(50:200, 10)
  )
  
  # Perform fixed effects meta-analysis
  results <- fixed_effects_meta(
    data = example_data,
    effect_col = "effect",
    variance_col = "variance",
    study_col = "study",
    method = "SMD"
  )
  
  # Display results
  cat("Fixed Effects Meta-Analysis Results\n")
  cat("====================================\n")
  summarize_fixed_effects(results)
  
  # Create forest plot
  plot_fixed_effects_forest(results, title = "Example Fixed Effects Meta-Analysis")
  
  return(results)
}

# ============================================================================
# END OF FIXED EFFECTS META-ANALYSIS MODULE
# ============================================================================
