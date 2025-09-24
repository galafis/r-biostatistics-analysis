# =============================================================================
# Randomized Clinical Trials Analysis
# Part of the r-biostatistics-analysis platform
# 
# This module provides comprehensive functions for analyzing randomized
# controlled trials (RCTs) including efficacy, safety, and non-inferiority
# analyses.
# =============================================================================

library(survival)
library(ggplot2)
library(dplyr)
library(tableone)
library(broom)

# -----------------------------------------------------------------------------
# Primary Efficacy Analysis
# -----------------------------------------------------------------------------

#' Analyze primary efficacy endpoint in RCT
#' 
#' @param data Data frame containing trial data
#' @param endpoint Primary endpoint variable name
#' @param treatment Treatment group variable name
#' @param covariates Optional vector of covariate names for adjusted analysis
#' @param alpha Significance level (default: 0.05)
#' @return List containing model results and treatment effect estimates
#' @export
efficacy_analysis <- function(data, endpoint, treatment, covariates = NULL, alpha = 0.05) {
  
  # Input validation
  if (!endpoint %in% names(data)) stop("Endpoint variable not found in data")
  if (!treatment %in% names(data)) stop("Treatment variable not found in data")
  
  if (is.null(covariates)) {
    # Unadjusted analysis
    formula_str <- paste(endpoint, "~", treatment)
    model <- lm(as.formula(formula_str), data = data)
  } else {
    # Adjusted analysis
    if (!all(covariates %in% names(data))) {
      stop("Some covariate variables not found in data")
    }
    formula_str <- paste(endpoint, "~", treatment, "+", 
                        paste(covariates, collapse = " + "))
    model <- lm(as.formula(formula_str), data = data)
  }
  
  # Extract results
  model_summary <- summary(model)
  treatment_coef <- model_summary$coefficients[2, ]
  ci <- confint(model, level = 1 - alpha)[2, ]
  
  # Effect size and Cohen's d
  pooled_sd <- sqrt(mean(residuals(model)^2))
  cohens_d <- treatment_coef[1] / pooled_sd
  
  return(list(
    model = model,
    treatment_effect = treatment_coef[1],
    standard_error = treatment_coef[2],
    t_statistic = treatment_coef[3],
    p_value = treatment_coef[4],
    ci_lower = ci[1],
    ci_upper = ci[2],
    cohens_d = cohens_d,
    formula = formula_str,
    sample_size = nrow(data),
    adjusted = !is.null(covariates)
  ))
}

# -----------------------------------------------------------------------------
# Non-inferiority Testing
# -----------------------------------------------------------------------------

#' Perform non-inferiority test
#' 
#' @param treatment_diff Treatment difference estimate
#' @param se_diff Standard error of treatment difference
#' @param margin Non-inferiority margin
#' @param df Degrees of freedom
#' @param alpha One-sided significance level (default: 0.025)
#' @return List containing test results
#' @export
non_inferiority_test <- function(treatment_diff, se_diff, margin, df, alpha = 0.025) {
  
  # Test statistic: H0: treatment_diff <= -margin vs H1: treatment_diff > -margin
  t_stat <- (treatment_diff + margin) / se_diff
  p_value <- 1 - pt(t_stat, df)
  
  # Critical value
  t_critical <- qt(1 - alpha, df)
  
  # Confidence interval
  ci_lower <- treatment_diff - qt(1 - alpha, df) * se_diff
  
  # Decision
  non_inferior <- p_value < alpha && ci_lower > -margin
  
  return(list(
    treatment_difference = treatment_diff,
    margin = margin,
    t_statistic = t_stat,
    p_value = p_value,
    t_critical = t_critical,
    ci_lower = ci_lower,
    non_inferior = non_inferior,
    conclusion = ifelse(non_inferior, "Non-inferior", "Inconclusive")
  ))
}

# -----------------------------------------------------------------------------
# Superiority Testing
# -----------------------------------------------------------------------------

#' Perform superiority test
#' 
#' @param treatment_diff Treatment difference estimate
#' @param se_diff Standard error of treatment difference
#' @param df Degrees of freedom
#' @param alpha Two-sided significance level (default: 0.05)
#' @return List containing test results
#' @export
superiority_test <- function(treatment_diff, se_diff, df, alpha = 0.05) {
  
  # Test statistic: H0: treatment_diff = 0 vs H1: treatment_diff ≠ 0
  t_stat <- treatment_diff / se_diff
  p_value <- 2 * (1 - pt(abs(t_stat), df))
  
  # Confidence interval
  t_critical <- qt(1 - alpha/2, df)
  ci_lower <- treatment_diff - t_critical * se_diff
  ci_upper <- treatment_diff + t_critical * se_diff
  
  # Decision
  superior <- p_value < alpha
  
  return(list(
    treatment_difference = treatment_diff,
    t_statistic = t_stat,
    p_value = p_value,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    superior = superior,
    conclusion = ifelse(superior, "Superior", "Not superior")
  ))
}

# -----------------------------------------------------------------------------
# Safety Analysis
# -----------------------------------------------------------------------------

#' Analyze safety endpoints in RCT
#' 
#' @param data Data frame containing safety data
#' @param treatment Treatment group variable
#' @param ae_vars Vector of adverse event variable names
#' @return List containing safety analysis results
#' @export
safety_analysis <- function(data, treatment, ae_vars) {
  
  safety_results <- list()
  
  for (ae in ae_vars) {
    if (!ae %in% names(data)) {
      warning(paste("AE variable", ae, "not found in data"))
      next
    }
    
    # Create contingency table
    tbl <- table(data[[treatment]], data[[ae]])
    
    # Fisher's exact test
    fisher_test <- fisher.test(tbl)
    
    # Risk difference and relative risk
    if (ncol(tbl) == 2 && nrow(tbl) == 2) {
      risk_ctrl <- tbl[1, 2] / sum(tbl[1, ])
      risk_trt <- tbl[2, 2] / sum(tbl[2, ])
      risk_diff <- risk_trt - risk_ctrl
      rel_risk <- ifelse(risk_ctrl > 0, risk_trt / risk_ctrl, NA)
    } else {
      risk_ctrl <- risk_trt <- risk_diff <- rel_risk <- NA
    }
    
    safety_results[[ae]] <- list(
      contingency_table = tbl,
      fisher_p_value = fisher_test$p.value,
      odds_ratio = fisher_test$estimate,
      or_ci_lower = fisher_test$conf.int[1],
      or_ci_upper = fisher_test$conf.int[2],
      risk_control = risk_ctrl,
      risk_treatment = risk_trt,
      risk_difference = risk_diff,
      relative_risk = rel_risk
    )
  }
  
  return(safety_results)
}

# -----------------------------------------------------------------------------
# Sample Size Calculation
# -----------------------------------------------------------------------------

#' Calculate sample size for RCT
#' 
#' @param effect_size Expected effect size (Cohen's d)
#' @param alpha Type I error rate (default: 0.05)
#' @param power Statistical power (default: 0.80)
#' @param ratio Allocation ratio (default: 1 for equal allocation)
#' @param two_sided Logical, two-sided test (default: TRUE)
#' @return List containing sample size calculations
#' @export
sample_size_rct <- function(effect_size, alpha = 0.05, power = 0.80, 
                           ratio = 1, two_sided = TRUE) {
  
  # Critical values
  if (two_sided) {
    z_alpha <- qnorm(1 - alpha/2)
  } else {
    z_alpha <- qnorm(1 - alpha)
  }
  z_beta <- qnorm(power)
  
  # Sample size calculation
  n1 <- ((z_alpha + z_beta)^2 * (1 + 1/ratio)) / effect_size^2
  n2 <- n1 * ratio
  n_total <- ceiling(n1 + n2)
  
  return(list(
    n_control = ceiling(n1),
    n_treatment = ceiling(n2),
    n_total = n_total,
    effect_size = effect_size,
    alpha = alpha,
    power = power,
    allocation_ratio = ratio,
    two_sided = two_sided
  ))
}

# -----------------------------------------------------------------------------
# Interim Analysis
# -----------------------------------------------------------------------------

#' Perform interim analysis with spending function
#' 
#' @param data Current data
#' @param endpoint Primary endpoint
#' @param treatment Treatment variable
#' @param target_n Target sample size
#' @param alpha_spend Alpha spending function value
#' @return List containing interim analysis results
#' @export
interim_analysis <- function(data, endpoint, treatment, target_n, alpha_spend) {
  
  current_n <- nrow(data)
  information_fraction <- current_n / target_n
  
  # Primary efficacy analysis
  efficacy_result <- efficacy_analysis(data, endpoint, treatment)
  
  # Adjusted critical value based on spending
  z_critical <- qnorm(1 - alpha_spend/2)
  z_observed <- abs(efficacy_result$t_statistic)
  
  # Futility assessment (conditional power)
  conditional_power <- calculate_conditional_power(
    current_effect = efficacy_result$treatment_effect,
    current_se = efficacy_result$standard_error,
    target_n = target_n,
    current_n = current_n
  )
  
  # Recommendations
  stop_efficacy <- z_observed > z_critical
  stop_futility <- conditional_power < 0.20
  
  return(list(
    information_fraction = information_fraction,
    current_n = current_n,
    target_n = target_n,
    efficacy_result = efficacy_result,
    z_critical = z_critical,
    z_observed = z_observed,
    conditional_power = conditional_power,
    stop_efficacy = stop_efficacy,
    stop_futility = stop_futility,
    recommendation = ifelse(stop_efficacy, "Stop for efficacy",
                           ifelse(stop_futility, "Stop for futility", "Continue"))
  ))
}

# Helper function for conditional power
calculate_conditional_power <- function(current_effect, current_se, target_n, current_n) {
  future_n <- target_n - current_n
  if (future_n <= 0) return(1)
  
  # Estimated final effect and SE
  final_se <- current_se * sqrt(current_n / target_n)
  z_final <- current_effect / final_se
  
  # Conditional power
  pnorm(z_final - 1.96)
}

# -----------------------------------------------------------------------------
# Utility Functions
# -----------------------------------------------------------------------------

#' Create summary table for RCT results
#' 
#' @param efficacy_result Result from efficacy_analysis()
#' @return Formatted data frame
#' @export
summarize_rct_results <- function(efficacy_result) {
  
  data.frame(
    Parameter = c("Treatment Effect", "Standard Error", "95% CI Lower", 
                 "95% CI Upper", "t-statistic", "p-value", "Cohen's d"),
    Value = c(
      round(efficacy_result$treatment_effect, 3),
      round(efficacy_result$standard_error, 3),
      round(efficacy_result$ci_lower, 3),
      round(efficacy_result$ci_upper, 3),
      round(efficacy_result$t_statistic, 3),
      format.pval(efficacy_result$p_value, digits = 3),
      round(efficacy_result$cohens_d, 3)
    ),
    stringsAsFactors = FALSE
  )
}

#' Generate RCT analysis report
#' 
#' @param data Trial data
#' @param endpoint Primary endpoint
#' @param treatment Treatment variable
#' @param covariates Optional covariates
#' @param safety_vars Safety variables
#' @return List containing comprehensive analysis
#' @export
generate_rct_report <- function(data, endpoint, treatment, 
                               covariates = NULL, safety_vars = NULL) {
  
  # Demographics
  demo_table <- CreateTableOne(
    vars = c(endpoint, covariates),
    strata = treatment,
    data = data
  )
  
  # Primary efficacy
  efficacy <- efficacy_analysis(data, endpoint, treatment, covariates)
  
  # Safety (if variables provided)
  safety <- if (!is.null(safety_vars)) {
    safety_analysis(data, treatment, safety_vars)
  } else NULL
  
  return(list(
    demographics = demo_table,
    efficacy = efficacy,
    safety = safety,
    summary_table = summarize_rct_results(efficacy)
  ))
}

# =============================================================================
# End of File
# =============================================================================
