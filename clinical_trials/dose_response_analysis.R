# =============================================================================
# Dose-Response Analysis - Sigmoid Emax model (compatible with dose_response_data.csv)
# =============================================================================
# Data expected:
# - Columns: dose (mg), response (e.g., % inhibition), subject_id (optional)
# - Doses tested: 5, 10, 20, 40 mg (per README)
# - Goal: Fit Emax model, estimate ED50, Emax, Hill (h), generate plot and table
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(broom)
  library(readr)
})

#' Read dose-response dataset
#' @param path CSV path (default: data/clinical_trials/dose_response_data.csv)
#' @return tibble
#' @export
read_dr_data <- function(path = "data/clinical_trials/dose_response_data.csv") {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(dose = as.numeric(dose))
}

#' Fit Sigmoid Emax model: y = E0 + (Emax * dose^h)/(ED50^h + dose^h)
#' @param data data frame with dose, response
#' @param start optional named list of starting values: E0, Emax, ED50, h
#' @return nls model and tidy summaries
#' @export
fit_emax <- function(data, start = NULL) {
  stopifnot(all(c("dose", "response") %in% names(data)))
  # Reasonable starts
  if (is.null(start)) {
    E0   <- min(data$response, na.rm = TRUE)
    Emax <- max(data$response, na.rm = TRUE) - E0
    ED50 <- median(unique(data$dose), na.rm = TRUE)
    h    <- 1
    start <- list(E0 = E0, Emax = Emax, ED50 = ED50, h = h)
  }
  mod <- nls(response ~ E0 + (Emax * dose^h)/(ED50^h + dose^h),
             data = data,
             start = start,
             control = nls.control(maxiter = 200, warnOnly = TRUE))
  list(
    model = mod,
    params = broom::tidy(mod),
    glance = broom::glance(mod)
  )
}

#' Predict response over a dose grid
#' @export
predict_emax <- function(fit, doses = seq(0, 50, by = 0.5)) {
  coefs <- as.list(coef(fit$model))
  tibble::tibble(
    dose = doses,
    response = coefs$E0 + (coefs$Emax * doses^coefs$h)/(coefs$ED50^coefs$h + doses^coefs$h)
  )
}

#' Plot dose-response with fitted curve
#' @export
plot_dose_response <- function(data, fit) {
  pred <- predict_emax(fit)
  ggplot(data, aes(dose, response)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_line(data = pred, aes(dose, response), color = "#1f77b4", linewidth = 1) +
    labs(title = "Sigmoid Emax Dose-Response",
         x = "Dose (mg)", y = "Response (%)") +
    theme_minimal()
}

#' EDx computation (dose for x% of Emax over E0)
#' @param x target percent of maximal effect (e.g., 50 for ED50)
#' @export
compute_EDx <- function(fit, x = 50) {
  coefs <- as.list(coef(fit$model))
  # Solve for dose where effect reaches E0 + x% of Emax
  target <- coefs$E0 + (x/100) * coefs$Emax
  # Invert Emax model -> dose = ED50 * (x/(1-x))^(1/h)
  ratio <- (x/100)/(1 - x/100)
  dose <- coefs$ED50 * ratio^(1/coefs$h)
  tibble::tibble(x = x, dose = dose)
}

#' End-to-end analysis helper
#' @export
run_dose_response <- function(path = "data/clinical_trials/dose_response_data.csv") {
  dat <- read_dr_data(path)
  fit <- fit_emax(dat)
  list(
    data = dat,
    fit = fit,
    plot = plot_dose_response(dat, fit),
    ED50 = compute_EDx(fit, 50),
    ED80 = compute_EDx(fit, 80)
  )
}

# Example
# out <- run_dose_response()
# out$plot
# out$fit$params
