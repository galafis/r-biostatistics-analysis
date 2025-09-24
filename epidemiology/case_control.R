#' Case-Control Study Analysis
#'
#' Provides odds ratio estimation with exact/fisher options, and adjusted logistic regression.
#' Requires: epiR, epitools, stats
#' @export
case_control_analysis <- function(data, exposure, outcome, conf.level = 0.95,
                                  exact = FALSE, adjusted_covariates = NULL) {
  stopifnot(is.data.frame(data))
  if (!all(c(exposure, outcome) %in% names(data))) stop("Variables not in data")

  # 2x2 table and OR
  tab <- table(exposed = data[[exposure]], case = data[[outcome]])
  or <- tryCatch({
    if (exact) {
      epitools::oddsratio(tab, method = "fisher")
    } else {
      epiR::epi.2by2(tab, method = "case.control", conf.level = conf.level)
    }
  }, error = function(e) NULL)

  fit <- NULL
  if (!is.null(adjusted_covariates)) {
    covars <- paste(adjusted_covariates, collapse = " + ")
    fml <- stats::as.formula(paste(outcome, "~", exposure, if (nzchar(covars)) paste("+", covars) else ""))
    fit <- stats::glm(fml, family = binomial(link = "logit"), data = data)
  }

  list(
    table = tab,
    odds_ratio = or,
    adjusted_model = fit
  )
}
