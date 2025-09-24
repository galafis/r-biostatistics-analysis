#' Causal Inference Utilities
#'
#' Provides doubly robust AIPW estimator and IV via 2SLS.
#' Requires: stats, AER (optional)
#' @export
#' @param data data.frame
#' @param treatment binary treatment variable name
#' @param outcome outcome variable name
#' @param covariates vector of covariate names
#' @return list with estimates and models
causal_inference <- function(data, treatment, outcome, covariates) {
  stopifnot(is.data.frame(data))
  if (!all(c(treatment, outcome, covariates) %in% names(data))) stop("Variables not in data")

  # Propensity model
  ps_fml <- as.formula(paste(treatment, "~", paste(covariates, collapse = " + ")))
  m_ps <- stats::glm(ps_fml, data = data, family = binomial())
  ps <- stats::predict(m_ps, type = "response")

  # Outcome models (treated/untreated)
  y1_fml <- as.formula(paste(outcome, "~", paste(covariates, collapse = " + ")))
  y0_fml <- y1_fml
  m_y1 <- stats::glm(y1_fml, data = subset(data, data[[treatment]] == 1))
  m_y0 <- stats::glm(y0_fml, data = subset(data, data[[treatment]] == 0))
  mu1 <- stats::predict(m_y1, newdata = data, type = "response")
  mu0 <- stats::predict(m_y0, newdata = data, type = "response")

  # AIPW estimator (ATE)
  A <- data[[treatment]]
  Y <- data[[outcome]]
  aipw <- mean(mu1 - mu0 + A * (Y - mu1) / pmax(ps, 1e-6) - (1 - A) * (Y - mu0) / pmax(1 - ps, 1e-6))

  list(ate_aipw = aipw, ps_model = m_ps, y1_model = m_y1, y0_model = m_y0)
}

#' Instrumental Variables via 2SLS (if AER available)
#' @export
iv_2sls <- function(data, treatment, outcome, instrument, covariates = NULL) {
  if (!requireNamespace("AER", quietly = TRUE)) stop("AER package required for IV")
  rhs <- paste(c(instrument, covariates), collapse = " + ")
  fml <- as.formula(paste(outcome, "~", treatment, "|", rhs))
  fit <- AER::ivreg(fml, data = data)
  summary(fit)
}
