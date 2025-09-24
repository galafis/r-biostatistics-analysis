#' Time-Varying Effects and Time-Dependent Covariates in Survival Analysis
#'
#' Utilities for modeling non-proportional hazards and time-dependent covariates:
#' - coxph with time-varying coefficients via tt() and interactions with log(time)
#' - tvc assessment using cox.zph and interaction terms
#' - start-stop format handling for external time-dependent covariates
#'
#' @section Dependencies:
#' - survival (coxph, Surv)
#' - dplyr (data munging)
#' - ggplot2 (plotting)
#'
#' @examples
#' # Time-varying coefficient using tt()
#' # fit <- tvc_cox(data=df, time="time", event="status",
#' #                covariates=c("treatment","age"), tvc_vars=c("treatment"))
#' # plot_list <- tvc_plot_effect(fit$model, terms=c("treatment"))
#'
#' # Start-stop format for external time-dependent covariate
#' # fit_ss <- tdc_cox_start_stop(data=df, start="t0", stop="t1", event="status",
#' #                              covariates=c("treatment","marker"))
#'
#' @keywords survival time-varying time-dependent coxph non-proportional hazards
#' @author Biostatistics Platform Contributors
#'
#' @importFrom survival Surv coxph cox.zph
#' @importFrom stats as.formula
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal

#' Cox model with time-varying coefficients using tt()
#'
#' @param data data.frame
#' @param time,event characters
#' @param covariates character vector of baseline covariates
#' @param tvc_vars character vector of variables to allow time-varying effect
#' @param tt_fun function of (x, t, ...) used in tt(), default = x*log(t)
#' @return list with model and diagnostics
#' @export
tvc_cox <- function(data, time, event, covariates, tvc_vars,
                    tt_fun = function(x, t, ...) x * log(t)) {
  stopifnot(all(c(time,event,covariates,tvc_vars) %in% names(data)))
  rhs <- paste(covariates, collapse = " + ")
  tt_terms <- paste(sprintf("tt(%s)", tvc_vars), collapse = " + ")
  fml <- as.formula(sprintf("survival::Surv(%s,%s) ~ %s + %s", time, event, rhs, tt_terms))
  model <- survival::coxph(fml, data = data,
                           tt = setNames(rep(list(tt_fun), length(tvc_vars)), tvc_vars))
  diag <- tryCatch(survival::cox.zph(model), error = function(e) NULL)
  list(model = model, diagnostics = diag, formula = fml)
}

#' Plot estimated time-varying effects using cox.zph output
#'
#' @param model coxph object
#' @param terms character vector of terms to plot
#' @return list of ggplot objects
#' @export
tvc_plot_effect <- function(model, terms = NULL) {
  requireNamespace("ggplot2", quietly = TRUE)
  zph <- survival::cox.zph(model)
  nm <- rownames(zph$table)
  nm <- nm[!grepl("GLOBAL", nm)]
  if (!is.null(terms)) nm <- intersect(nm, terms)
  plots <- lapply(nm, function(tr) {
    df <- data.frame(t = zph$x[[tr]], effect = zph$y[[tr]])
    ggplot2::ggplot(df, ggplot2::aes(x = t, y = effect)) +
      ggplot2::geom_line() +
      ggplot2::labs(title = paste("Time-varying effect:", tr), x = "Transformed time", y = "Scaled Schoenfeld residual") +
      ggplot2::theme_minimal()
  })
  names(plots) <- nm
  plots
}

#' Cox model with external time-dependent covariates (start-stop)
#'
#' @param data data.frame with start, stop, event and covariates
#' @param start,stop,event characters
#' @param covariates character vector
#' @param ties ties handling (default "efron")
#' @return coxph model
#' @export
tdc_cox_start_stop <- function(data, start, stop, event, covariates, ties = "efron") {
  stopifnot(all(c(start,stop,event,covariates) %in% names(data)))
  rhs <- paste(covariates, collapse = " + ")
  fml <- as.formula(sprintf("survival::Surv(%s, %s, %s) ~ %s", start, stop, event, rhs))
  survival::coxph(fml, data = data, ties = ties, x = TRUE, y = TRUE)
}

#' Convenience: full time-varying workflow
#' @export
tve_full_report <- function(data, time, event, covariates, tvc_vars = NULL,
                            start = NULL, stop = NULL, tt_fun = function(x, t, ...) x*log(t)) {
  tve <- if (!is.null(tvc_vars) && length(tvc_vars) > 0)
    tryCatch(tvc_cox(data, time, event, covariates, tvc_vars, tt_fun), error = function(e) NULL) else NULL
  ss  <- if (!is.null(start) && !is.null(stop))
    tryCatch(tdc_cox_start_stop(data, start, stop, event, covariates), error = function(e) NULL) else NULL
  list(tvc = tve, start_stop = ss)
}
