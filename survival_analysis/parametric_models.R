#' Parametric Survival Models: Flexible Distributions and Reporting
#'
#' Functions to fit parametric survival models (exponential, Weibull, log-normal,
#' log-logistic, Gompertz) using flexsurv, extract tidy summaries, compare models
#' via AIC/BIC, predict survival and hazard, and produce ready-to-report tables.
#'
#' @section Dependencies:
#' - flexsurv (flexsurvreg, summaries)
#' - survival (Surv)
#' - broom, dplyr, tidyr
#'
#' @examples
#' # fit suite of distributions
#' # mods <- ps_fit_set(df, time = "time", event = "status",
#' #                    covariates = c("treatment","age"),
#' #                    dists = c("weibull","lnorm","llogis"))
#' # cmp  <- ps_compare(mods)
#' # pred <- ps_predict_survival(mods$weibull$model, newdata = df[1:2,], times = c(6,12,24))
#'
#' @keywords survival parametric flexsurv weibull lnorm llogis gompertz reporting
#' @author Biostatistics Platform Contributors
#'
#' @importFrom survival Surv
#' @importFrom stats as.formula
#'
#' Fit a parametric survival model with flexsurv
#'
#' @param data data.frame
#' @param time,event characters
#' @param covariates character vector of covariates (can be NULL for intercept-only)
#' @param dist character distribution supported by flexsurvreg
#' @return list with model (flexsurvreg), formula, args
#' @export
ps_fit <- function(data, time, event, covariates = NULL, dist = "weibull") {
  requireNamespace("flexsurv", quietly = TRUE)
  stopifnot(time %in% names(data), event %in% names(data))
  if (!is.null(covariates)) stopifnot(all(covariates %in% names(data)))
  rhs <- if (is.null(covariates) || length(covariates) == 0) "1" else paste(covariates, collapse = " + ")
  fml <- as.formula(sprintf("survival::Surv(%s,%s) ~ %s", time, event, rhs))
  model <- flexsurv::flexsurvreg(formula = fml, data = data, dist = dist)
  structure(list(model = model, formula = fml, args = list(time=time,event=event,covariates=covariates,dist=dist)),
            class = c("ps_fit_result", class(model)))
}

#' Fit a set of parametric distributions
#' @param dists character vector (e.g., c("exp","weibull","lnorm","llogis","gompertz"))
#' @return named list of ps_fit_result
#' @export
ps_fit_set <- function(data, time, event, covariates = NULL,
                       dists = c("exp","weibull","lnorm","llogis","gompertz")) {
  res <- lapply(dists, function(d) {
    try(ps_fit(data, time, event, covariates, dist = d), silent = TRUE)
  })
  names(res) <- dists
  res
}

#' Compare parametric models by information criteria
#' @param models list as from ps_fit_set
#' @return data.frame with dist, AIC, BIC, logLik
#' @export
ps_compare <- function(models) {
  keep <- models[ sapply(models, function(x) inherits(x$model, "flexsurvreg")) ]
  df <- do.call(rbind, lapply(names(keep), function(nm) {
    m <- keep[[nm]]$model
    data.frame(dist = nm, AIC = AIC(m), BIC = BIC(m), logLik = as.numeric(logLik(m)), row.names = NULL)
  }))
  df[order(df$AIC), ]
}

#' Tidy table of parameters and hazard ratios (if applicable)
#' @param model flexsurvreg object
#' @return data.frame
#' @export
ps_tidy <- function(model) {
  stopifnot(inherits(model, "flexsurvreg"))
  s <- summary(model)
  out <- do.call(rbind, lapply(s, function(comp) {
    as.data.frame(comp)
  }))
  rownames(out) <- NULL
  out
}

#' Predict survival for newdata at times
#' @param model flexsurvreg
#' @param newdata data.frame
#' @param times numeric vector
#' @return matrix/data.frame of survival probabilities
#' @export
ps_predict_survival <- function(model, newdata, times) {
  stopifnot(inherits(model, "flexsurvreg"))
  pr <- flexsurv::summary.flexsurvreg(model, newdata = newdata, type = "survival", t = times)
  # pr is a list per row; coerce to matrix
  mat <- do.call(rbind, lapply(pr, function(lst) sapply(lst, function(x) x$est)))
  colnames(mat) <- paste0("t=", times)
  mat
}

#' Predict hazard for newdata at times
#' @export
ps_predict_hazard <- function(model, newdata, times) {
  stopifnot(inherits(model, "flexsurvreg"))
  pr <- flexsurv::summary.flexsurvreg(model, newdata = newdata, type = "hazard", t = times)
  mat <- do.call(rbind, lapply(pr, function(lst) sapply(lst, function(x) x$est)))
  colnames(mat) <- paste0("t=", times)
  mat
}

#' Convenience: full parametric workflow (fit set, compare, best model)
#' @export
ps_full_report <- function(data, time, event, covariates = NULL,
                           dists = c("exp","weibull","lnorm","llogis","gompertz")) {
  fits <- ps_fit_set(data, time, event, covariates, dists)
  cmp  <- tryCatch(ps_compare(fits), error = function(e) NULL)
  best <- NULL
  if (!is.null(cmp) && nrow(cmp) > 0) {
    best <- fits[[ cmp$dist[1] ]]
  }
  list(fits = fits, compare = cmp, best = best)
}
