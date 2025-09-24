# =============================================================================
# ANCOVA Analysis for RCT primary endpoint (compatible with rct_data.csv)
# =============================================================================
# Example in README: endpoint_score ~ treatment + baseline_score
# This script provides a robust, reusable ANCOVA workflow.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(broom)
  library(readr)
})

#' Read RCT efficacy dataset
#' @param path CSV path (default: data/clinical_trials/rct_data.csv)
#' @return tibble
#' @export
read_rct_data <- function(path = "data/clinical_trials/rct_data.csv") {
  readr::read_csv(path, show_col_types = FALSE)
}

#' Fit ANCOVA model with baseline adjustment
#' @param data data frame
#' @param endpoint character, endpoint variable (e.g., "endpoint_score")
#' @param treatment character, treatment factor (e.g., "treatment")
#' @param baseline character, baseline covariate (e.g., "baseline_score")
#' @param alpha numeric, two-sided alpha (0.05)
#' @return list with model, estimates, CI, adjusted means by group
#' @export
fit_ancova <- function(data, endpoint, treatment, baseline, alpha = 0.05) {
  stopifnot(endpoint %in% names(data), treatment %in% names(data), baseline %in% names(data))
  data[[treatment]] <- factor(data[[treatment]], levels = sort(unique(as.character(data[[treatment]]))))
  formula_str <- paste(endpoint, "~", treatment, "+", baseline)
  mod <- lm(as.formula(formula_str), data = data)

  smry <- summary(mod)
  cf   <- smry$coefficients
  # Assumes first treatment level is reference (e.g., Placebo)
  trt_row <- grep(paste0("^", treatment), rownames(cf))
  if (length(trt_row) == 0) {
    # Alternative: term is 'treatmentActive'
    trt_row <- which(grepl("^" %+% treatment, rownames(cf)))
  }
  est <- cf[trt_row[1], "Estimate"]
  se  <- cf[trt_row[1], "Std. Error"]
  df  <- smry$df[2]
  tcrit <- qt(1 - alpha/2, df)
  ci <- c(est - tcrit * se, est + tcrit * se)

  # Adjusted means at overall baseline mean
  mu_base <- mean(data[[baseline]], na.rm = TRUE)
  groups <- levels(data[[treatment]])
  newdata <- data.frame(
    baseline_tmp = mu_base
  )
  # Build adjusted means by group
  adj_means <- sapply(groups, function(g) {
    nd <- data.frame(baseline = mu_base)
    nd[[treatment]] <- factor(g, levels = groups)
    as.numeric(predict(mod, newdata = nd))
  })

  list(
    model = mod,
    treatment_contrast = names(adj_means)[2],
    treatment_effect = est,
    se = se,
    df = df,
    ci_lower = ci[1],
    ci_upper = ci[2],
    p_value = 2 * (1 - pt(abs(est/se), df)),
    adjusted_means = tibble::tibble(group = groups, adj_mean = as.numeric(adj_means)),
    formula = formula_str
  )
}

#' Diagnostic plots for ANCOVA
#' @export
ancova_diagnostics <- function(fit) {
  p1 <- ggplot2::ggplot(data.frame(resid = resid(fit$model), fitted = fitted(fit$model)),
                        aes(fitted, resid)) +
    geom_point(alpha = 0.6) + geom_hline(yintercept = 0, linetype = 2) +
    labs(title = "Residuals vs Fitted") + theme_minimal()
  p2 <- ggplot2::ggplot(data.frame(resid = resid(fit$model)), aes(sample = resid)) +
    stat_qq() + stat_qq_line() + theme_minimal() + labs(title = "Normal Q-Q")
  list(resid_fitted = p1, qq = p2)
}

#' End-to-end helper
#' @export
run_ancova <- function(path = "data/clinical_trials/rct_data.csv",
                       endpoint = "endpoint_score",
                       treatment = "treatment",
                       baseline = "baseline_score") {
  dat <- read_rct_data(path)
  fit <- fit_ancova(dat, endpoint, treatment, baseline)
  list(data = dat, fit = fit, diagnostics = ancova_diagnostics(fit))
}

# Example
# out <- run_ancova()
# out$fit$adjusted_means
# out$diagnostics$resid_fitted
