#' Propensity Score Methods
#'
#' Fits propensity scores and performs matching/weighting diagnostics.
#' Requires: MatchIt, stats, tableone
#' @export
ps_analysis <- function(data, treatment, covariates, method = c("matching","weighting"),
                        distance = c("logit","probit"), ratio = 1) {
  stopifnot(is.data.frame(data))
  method <- match.arg(method)
  distance <- match.arg(distance)
  if (!all(c(treatment, covariates) %in% names(data))) stop("Variables not in data")
  formula_str <- paste(treatment, "~", paste(covariates, collapse = " + "))
  fml <- as.formula(formula_str)
  if (method == "matching") {
    if (!requireNamespace("MatchIt", quietly = TRUE)) stop("MatchIt required")
    m <- MatchIt::matchit(fml, data = data, method = "nearest", distance = distance, ratio = ratio)
    matched <- MatchIt::match.data(m)
    bal <- NULL
    if (requireNamespace("tableone", quietly = TRUE)) {
      bal <- tableone::CreateTableOne(vars = covariates, strata = treatment, data = matched, test = FALSE)
    }
    list(matchit = m, matched_data = matched, balance = bal)
  } else {
    ps_model <- stats::glm(fml, data = data, family = binomial(link = distance))
    ps <- stats::predict(ps_model, type = "response")
    trt <- data[[treatment]]
    w <- ifelse(trt == 1, 1/ps, 1/(1-ps))
    list(model = ps_model, ps = ps, weights = w)
  }
}
