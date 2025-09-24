#' Cohort Studies Analysis
#'
#' Provides cohort analysis utilities: incidence rates, risk ratios, adjusted Poisson/log-binomial models.
#' Requires: epiR, epitools, stats, utils
#' @export
cohort_analysis <- function(data, id = NULL, exposure, outcome, time = NULL,
                            conf.level = 0.95, adjusted_covariates = NULL,
                            model = c("poisson", "logbin")) {
  stopifnot(is.data.frame(data))
  model <- match.arg(model)
  if (!all(c(exposure, outcome) %in% names(data))) stop("Variables not in data")

  # Basic 2x2 table
  tab <- table(exposed = data[[exposure]], outcome = data[[outcome]])

  # epiR risk ratio with CI
  rr <- tryCatch({
    epiR::epi.2by2(tab, method = "cohort.count", conf.level = conf.level)
  }, error = function(e) NULL)

  # Incidence rate if time provided
  ir <- NULL
  if (!is.null(time) && time %in% names(data)) {
    # person-time by exposure
    pt_exp <- sum(data[[time]][data[[exposure]] == 1], na.rm = TRUE)
    pt_unexp <- sum(data[[time]][data[[exposure]] == 0], na.rm = TRUE)
    events_exp <- sum(data[[outcome]] == 1 & data[[exposure]] == 1, na.rm = TRUE)
    events_unexp <- sum(data[[outcome]] == 1 & data[[exposure]] == 0, na.rm = TRUE)
    rate_exp <- events_exp / pt_exp
    rate_unexp <- events_unexp / pt_unexp
    ir <- list(rate_exp = rate_exp, rate_unexp = rate_unexp, rate_ratio = rate_exp/rate_unexp)
  }

  # Adjusted model
  fit <- NULL
  if (!is.null(adjusted_covariates)) {
    covars <- paste(adjusted_covariates, collapse = " + ")
    fml <- stats::as.formula(paste(outcome, "~", exposure, if (nzchar(covars)) paste("+", covars) else ""))
    if (model == "poisson") {
      # Robust Poisson for RR
      fit <- stats::glm(fml, family = poisson(link = "log"), data = data)
      # robust SE via sandwich if available
      if (requireNamespace("sandwich", quietly = TRUE) && requireNamespace("lmtest", quietly = TRUE)) {
        cov <- sandwich::vcovCL(fit, cluster = if (!is.null(id) && id %in% names(data)) data[[id]] else NULL)
        fit$robust_se <- sqrt(diag(cov))
        fit$robust_ci <- cbind(
          est = coef(fit),
          lcl = coef(fit) - qnorm(1 - (1 - conf.level)/2) * fit$robust_se,
          ucl = coef(fit) + qnorm(1 - (1 - conf.level)/2) * fit$robust_se
        )
      }
    } else {
      # Log-binomial (may fail to converge)
      fit <- tryCatch(stats::glm(fml, family = binomial(link = "log"), data = data), error = function(e) NULL)
      if (is.null(fit)) {
        # fallback to Poisson with robust SE
        fit <- stats::glm(fml, family = poisson(link = "log"), data = data)
      }
    }
  }

  list(
    table = tab,
    risk_ratio = rr,
    incidence_rates = ir,
    adjusted_model = fit
  )
}

#' Kaplan-Meier-like cumulative incidence by exposure (for descriptive plots)
#' @export
cohort_km <- function(time, event, exposure, conf.int = TRUE) {
  stopifnot(length(time) == length(event), length(event) == length(exposure))
  if (!requireNamespace("survival", quietly = TRUE)) stop("survival package required")
  df <- data.frame(time = time, event = event, exposure = exposure)
  s <- survival::survfit(survival::Surv(time, event) ~ exposure, data = df, conf.int = conf.int)
  s
}
