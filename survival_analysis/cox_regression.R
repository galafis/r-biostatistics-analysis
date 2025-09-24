#' Cox Proportional Hazards: Multivariable Modeling and Diagnostics
#'
#' Professional functions to fit Cox models, extract hazard ratios,
#' perform proportional hazards diagnostics, handle ties methods,
#' compute baseline survival, and export tidy results and plots.
#'
#' @section Dependencies:
#' - survival (coxph, cox.zph)
#' - broom (tidy)
#' - dplyr, tidyr (tidy outputs)
#' - ggplot2 (diagnostic plots)
#'
#' @examples
#' # Example (df with time, status, covariates including treatment):
#' # fit <- cox_fit(data = df, time = "time", event = "status",
#' #                covariates = c("treatment", "age", "sex"))
#' # diag <- cox_diagnostics(fit$model)
#' # hr   <- cox_hr_table(fit$model)
#' # base <- cox_baseline_survival(fit$model)
#' # rep  <- cox_full_report(data = df, time = "time", event = "status",
#' #                         covariates = c("treatment","age","sex"))
#'
#' @keywords survival cox proportional-hazards diagnostics hr biostatistics
#' @author Biostatistics Platform Contributors
#
#' @importFrom survival Surv coxph cox.zph basehaz
#' @importFrom stats as.formula model.matrix
#' @importFrom broom tidy glance
#' @importFrom dplyr mutate select arrange rename everything
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_minimal

#' Fit a multivariable Cox proportional hazards model
#'
#' @param data data.frame
#' @param time character, time-to-event variable
#' @param event character, event indicator (1=event, 0=censor)
#' @param covariates character vector of covariate names
#' @param strata optional character vector of variables to stratify baseline hazard
#' @param ties character, method for ties ("efron","breslow","exact")
#' @return list with model (coxph), formula, input args
#' @export
cox_fit <- function(data, time, event, covariates, strata = NULL, ties = "efron") {
  stopifnot(is.data.frame(data), time %in% names(data), event %in% names(data))
  stopifnot(all(covariates %in% names(data)))
  if (!is.null(strata)) stopifnot(all(strata %in% names(data)))
  rhs <- paste(covariates, collapse = " + ")
  if (!is.null(strata) && length(strata) > 0) {
    rhs <- paste(rhs, "+", paste(sprintf("strata(%s)", strata), collapse = " + "))
  }
  fml <- as.formula(sprintf("survival::Surv(%s,%s) ~ %s", time, event, rhs))
  model <- survival::coxph(fml, data = data, ties = ties, x = TRUE, y = TRUE, model = TRUE)
  structure(list(model = model, formula = fml,
                 args = list(time=time, event=event, covariates=covariates, strata=strata, ties=ties)),
            class = c("cox_fit_result", class(model)))
}

#' Extract tidy hazard ratio table with confidence intervals
#'
#' @param model coxph object
#' @param conf_level numeric, confidence level (default 0.95)
#' @return data.frame with term, HR, lower, upper, p_value
#' @export
cox_hr_table <- function(model, conf_level = 0.95) {
  stopifnot(inherits(model, "coxph"))
  est <- stats::coef(model)
  se  <- sqrt(diag(model$var))
  z   <- qnorm((1 + conf_level)/2)
  hr  <- exp(est)
  lower <- exp(est - z*se)
  upper <- exp(est + z*se)
  pval  <- 2*pnorm(abs(est/se), lower.tail = FALSE)
  out <- data.frame(term = names(est), HR = hr, lower = lower, upper = upper, p_value = pval, row.names = NULL)
  rownames(out) <- NULL
  out
}

#' Proportional hazards diagnostics (Schoenfeld residuals)
#'
#' @param model coxph object
#' @return list with cox.zph object and tidy table of global/term tests
#' @export
cox_diagnostics <- function(model) {
  stopifnot(inherits(model, "coxph"))
  zph <- survival::cox.zph(model)
  tbl <- data.frame(term = rownames(zph$table),
                    rho = zph$table[,"rho"],
                    chisq = zph$table[,"chisq"],
                    p_value = zph$table[,"p"], row.names = NULL)
  list(zph = zph, tidy = tbl)
}

#' Plot Schoenfeld residuals vs. transformed time for a covariate
#'
#' @param diag_list output of cox_diagnostics() or cox.zph object
#' @param term character, covariate term to plot; if NULL, returns list of plots
#' @return ggplot object or named list of plots
#' @export
cox_plot_schoenfeld <- function(diag_list, term = NULL) {
  requireNamespace("ggplot2", quietly = TRUE)
  zph <- if (inherits(diag_list, "cox.zph")) diag_list else diag_list$zph
  df_list <- lapply(names(zph$x), function(nm) {
    data.frame(term = nm, x = zph$x[[nm]], y = zph$y[[nm]])
  })
  dff <- do.call(rbind, df_list)
  make_plot <- function(df) {
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::geom_smooth(method = "loess", se = TRUE) +
      ggplot2::labs(x = "Transformed time", y = "Schoenfeld residual", title = unique(df$term)) +
      ggplot2::theme_minimal()
  }
  if (!is.null(term)) {
    stopifnot(term %in% unique(dff$term))
    return(make_plot(dff[dff$term == term, , drop = FALSE]))
  }
  split_df <- split(dff, dff$term)
  lapply(split_df, make_plot)
}

#' Baseline cumulative hazard and survival
#'
#' @param model coxph object
#' @return data.frame with time, cumhaz, surv
#' @export
cox_baseline_survival <- function(model) {
  stopifnot(inherits(model, "coxph"))
  bh <- survival::basehaz(model, centered = FALSE)
  # baseline survival S0(t) = exp(-H0(t))
  data.frame(time = bh$time, cumhaz = bh$hazard, surv = exp(-bh$hazard))
}

#' Predict individual survival curves for newdata
#'
#' @param model coxph object
#' @param newdata data.frame with covariates
#' @param times optional numeric vector of times
#' @return matrix/data.frame of survival probabilities (rows=newdata, cols=times)
#' @export
cox_predict_survival <- function(model, newdata, times = NULL) {
  stopifnot(inherits(model, "coxph"))
  sfit <- survfit(model, newdata = newdata)
  if (!is.null(times)) {
    s <- summary(sfit, times = times)
    # For multiple rows, surv is a list; coerce to matrix
    if (is.null(s$strata)) {
      return(matrix(s$surv, nrow = nrow(newdata), byrow = TRUE, dimnames = list(NULL, paste0("t=", s$time))))
    } else {
      return(s$surv)
    }
  }
  sfit
}

#' Convenience: full Cox analysis workflow
#'
#' @param data,time,event,covariates,strata,ties see above
#' @return list with model, HR table, PH diagnostics, baseline survival
#' @export
cox_full_report <- function(data, time, event, covariates, strata = NULL, ties = "efron") {
  fit <- cox_fit(data, time, event, covariates, strata, ties)
  hr  <- tryCatch(cox_hr_table(fit$model), error = function(e) NULL)
  diag <- tryCatch(cox_diagnostics(fit$model), error = function(e) NULL)
  base <- tryCatch(cox_baseline_survival(fit$model), error = function(e) NULL)
  list(model = fit$model, hr = hr, diagnostics = diag, baseline = base)
}
