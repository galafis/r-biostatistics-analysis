#' Cross-Sectional Study Analysis
#'
#' Provides prevalence estimation, prevalence ratio via Poisson robust, and adjusted models.
#' Requires: epiR, stats
#' @export
cross_sectional_analysis <- function(data, exposure, outcome, conf.level = 0.95,
                                     adjusted_covariates = NULL) {
  stopifnot(is.data.frame(data))
  if (!all(c(exposure, outcome) %in% names(data))) stop("Variables not in data")

  # Prevalence by exposure
  prev_exp <- mean(data[[outcome]][data[[exposure]] == 1], na.rm = TRUE)
  prev_unexp <- mean(data[[outcome]][data[[exposure]] == 0], na.rm = TRUE)
  pr <- prev_exp / prev_unexp

  # 2x2 with epiR
  tab <- table(exposed = data[[exposure]], outcome = data[[outcome]])
  pr_epi <- tryCatch({
    epiR::epi.2by2(tab, method = "cross.sectional", conf.level = conf.level)
  }, error = function(e) NULL)

  fit <- NULL
  if (!is.null(adjusted_covariates)) {
    covars <- paste(adjusted_covariates, collapse = " + ")
    fml <- stats::as.formula(paste(outcome, "~", exposure, if (nzchar(covars)) paste("+", covars) else ""))
    fit <- stats::glm(fml, family = poisson(link = "log"), data = data)
    if (requireNamespace("sandwich", quietly = TRUE)) {
      V <- sandwich::vcovHC(fit, type = "HC0")
      fit$robust_se <- sqrt(diag(V))
    }
  }

  list(
    table = tab,
    prevalence = list(exposed = prev_exp, unexposed = prev_unexp, ratio = pr),
    prevalence_ratio = pr_epi,
    adjusted_model = fit
  )
}
