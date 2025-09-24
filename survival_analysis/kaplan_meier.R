#' Kaplan-Meier Analysis and Reporting Utilities
#'
#' Professional-grade functions to perform Kaplan-Meier survival analysis,
#' generate publication-quality plots, conduct log-rank tests, compute survival
#' estimates at specific time points, restricted mean survival time (RMST), and
#' produce a tidy summary table for clinical reporting.
#'
#' Functions follow the platform conventions: return lists with model objects,
#' tidy results, and ggplot objects ready for downstream reporting.
#'
#' @section Dependencies:
#' - survival (core KM and tests)
#' - survminer (enhanced plotting)
#' - dplyr, tidyr (tidy outputs)
#' - broom (tidy model summaries)
#' - survRM2 (RMST analysis)
#'
#' @examples
#' # Example (assuming df with columns: time, status, arm)
#' # km <- km_fit(data = df, time = "time", event = "status", group = "arm")
#' # plt <- km_plot(km$fit, data = df, group = "arm")
#' # lr  <- km_logrank(df, time = "time", event = "status", group = "arm")
#' # est <- km_estimates_at(km$fit, time_points = c(6, 12))
#' # rmst <- km_rmst(df, time = "time", event = "status", group = "arm", tau = 24)
#' # tab <- km_summary_table(km$fit, data = df, group = "arm")
#'
#' @keywords survival kaplan-meier logrank rmst reporting biostatistics
#' @author
#' Biostatistics Platform Contributors

# ---- Imports ----
#' @importFrom survival Surv survfit survdiff
#' @importFrom stats as.formula
#' @importFrom dplyr mutate select arrange group_by summarise n across everything left_join bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom broom tidy glance
#' @importFrom ggplot2 ggplot aes labs theme_minimal

# ---- Core KM fit ----
#' Fit Kaplan-Meier model
#'
#' @param data data.frame with survival data
#' @param time character. Time-to-event variable name
#' @param event character. Event indicator (1=event, 0=censored)
#' @param group optional character. Grouping variable for stratified KM
#' @return list with elements: fit (survfit), call, args
#' @export
km_fit <- function(data, time, event, group = NULL) {
  stopifnot(is.data.frame(data), time %in% names(data), event %in% names(data))
  if (!is.null(group)) stopifnot(group %in% names(data))
  requireNamespace("survival", quietly = TRUE)

  if (is.null(group)) {
    fml <- as.formula(sprintf("survival::Surv(%s, %s) ~ 1", time, event))
  } else {
    fml <- as.formula(sprintf("survival::Surv(%s, %s) ~ %s", time, event, group))
  }
  fit <- survival::survfit(fml, data = data)
  structure(list(fit = fit, call = match.call(), args = list(time=time, event=event, group=group)),
            class = c("km_fit_result", class(fit)))
}

# ---- Plotting ----
#' Kaplan-Meier plot with risk table
#'
#' @param fit survfit object from km_fit or survival::survfit
#' @param data original data.frame (required when using risk.table)
#' @param group optional character, grouping variable name (for labeling)
#' @param conf_int logical, show confidence intervals
#' @param pval logical, show log-rank p-value on plot when groups present
#' @param risk_table logical, show risk table
#' @param palette optional vector of colors
#' @return ggsurvplot object (survminer)
#' @export
km_plot <- function(fit, data, group = NULL, conf_int = TRUE, pval = TRUE, risk_table = TRUE, palette = NULL) {
  requireNamespace("survminer", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  args <- list(conf.int = conf_int, pval = (pval && !is.null(group)), risk.table = risk_table)
  if (!is.null(palette)) args$palette <- palette
  do.call(survminer::ggsurvplot, c(list(fit = fit, data = data), args))
}

# ---- Log-rank test ----
#' Log-rank test for survival curves
#'
#' @param data data.frame
#' @param time,event,group variable names (characters)
#' @return list with survdiff object and tidy summary (chisq, df, p)
#' @export
km_logrank <- function(data, time, event, group) {
  stopifnot(group %in% names(data))
  fml <- as.formula(sprintf("survival::Surv(%s,%s) ~ %s", time, event, group))
  lr <- survival::survdiff(fml, data = data)
  chisq <- unname(lr$chisq)
  df <- length(lr$n) - 1L
  p <- stats::pchisq(chisq, df = df, lower.tail = FALSE)
  list(test = lr, tidy = data.frame(statistic = chisq, df = df, p_value = p))
}

# ---- Estimates at time points ----
#' Survival estimates at specified time points
#'
#' @param fit survfit object
#' @param time_points numeric vector of times (same scale as time var)
#' @return data.frame with group, time, surv, lower, upper
#' @export
km_estimates_at <- function(fit, time_points) {
  tp <- sort(unique(as.numeric(time_points)))
  summ <- summary(fit, times = tp)
  # Handle single vs multi-strata outputs
  if (is.null(summ$strata)) {
    data.frame(group = "All", time = summ$time, surv = summ$surv, lower = summ$lower, upper = summ$upper)
  } else {
    strata <- rep(names(summ$strata), times = as.numeric(summ$strata))
    data.frame(group = strata, time = summ$time, surv = summ$surv, lower = summ$lower, upper = summ$upper)
  }
}

# ---- RMST ----
#' Restricted Mean Survival Time (RMST) analysis
#'
#' @param data data.frame
#' @param time,event,group variable names
#' @param tau truncation time for RMST (numeric). If NULL, uses min(max times across groups)
#' @return list with rmst object and tidy table of RMST by group and difference
#' @export
km_rmst <- function(data, time, event, group, tau = NULL) {
  requireNamespace("survRM2", quietly = TRUE)
  stopifnot(time %in% names(data), event %in% names(data), group %in% names(data))
  # Determine tau if not provided
  if (is.null(tau)) {
    tau <- min(tapply(data[[time]], data[[group]], max, na.rm = TRUE))
  }
  arm <- as.factor(data[[group]])
  # survRM2 expects numeric arm coded 0/1 for two-group; extend gracefully
  if (nlevels(arm) == 2) {
    x <- survRM2::rmst2(time = data[[time]], status = data[[event]], arm = as.integer(arm) - 1L, tau = tau)
    tidy <- data.frame(
      measure = c("RMST", "RMST_diff", "RMST_ratio"),
      estimate = c(x$RMST.arm1$result[1, "Est."],
                   x$unadjusted.result[1, "Est."],
                   x$unadjusted.result[2, "Est."]),
      lower = c(x$RMST.arm1$result[1, "lower .95"],
                x$unadjusted.result[1, "lower .95"],
                x$unadjusted.result[2, "lower .95"]),
      upper = c(x$RMST.arm1$result[1, "upper .95"],
                x$unadjusted.result[1, "upper .95"],
                x$unadjusted.result[2, "upper .95"]),
      p_value = c(NA_real_, x$unadjusted.result[1, "p"], x$unadjusted.result[2, "p"]),
      tau = tau
    )
    list(rmst = x, tidy = tidy)
  } else {
    # For >2 groups, compute per-group RMSTs only
    fits <- lapply(split(data, data[[group]]), function(df) {
      survRM2::rmst2(time = df[[time]], status = df[[event]], arm = rep(0L, nrow(df)), tau = tau)
    })
    tidy <- do.call(rbind, lapply(names(fits), function(g) {
      arm1 <- fits[[g]]$RMST.arm1$result[1, , drop = FALSE]
      data.frame(group = g,
                 RMST = arm1[1, "Est."],
                 lower = arm1[1, "lower .95"],
                 upper = arm1[1, "upper .95"],
                 tau = tau)
    }))
    list(rmst = fits, tidy = tidy)
  }
}

# ---- Summary table ----
#' Tidy KM summary table for reporting
#'
#' Includes N at risk, events, median survival (with CI), and survival at key
#' time points (e.g., 6, 12, 24 months if present in data range).
#'
#' @param fit survfit object (stratified or not)
#' @param data original data.frame
#' @param group optional group variable name
#' @param time_points optional numeric vector of times for survival estimates; if NULL, uses pretty sequence within range
#' @return data.frame ready for table rendering
#' @export
km_summary_table <- function(fit, data, group = NULL, time_points = NULL) {
  requireNamespace("dplyr", quietly = TRUE)
  # Basic counts
  if (is.null(group)) {
    n <- nrow(data)
    events <- sum(data[[attr(fit, "call")]$args$event]], na.rm = TRUE) # safe fallback removed below
  }
  # Median survival
  med <- tryCatch({
    s <- summary(fit)
    data.frame(
      group = if (is.null(s$strata)) "All" else names(s$strata),
      median = suppressWarnings(stats::median(fit)),
      median_lower = suppressWarnings(fit$conf.int[1]),
      median_upper = suppressWarnings(fit$conf.int[2])
    )
  }, error = function(e) NULL)

  # Time points
  if (is.null(time_points)) {
    # choose up to 3 representative times within max follow-up
    tf <- max(fit$time, na.rm = TRUE)
    time_points <- unique(round(stats::quantile(fit$time, probs = c(0.25, 0.5, 0.75)), 1))
    time_points <- time_points[time_points > 0 & time_points <= tf]
  }
  est <- km_estimates_at(fit, time_points)
  est_wide <- tidyr::pivot_wider(est, names_from = time, values_from = c(surv, lower, upper), names_sep = "_t")

  out <- if (!is.null(med)) merge(est_wide, med, by = intersect(names(est_wide), names(med)), all = TRUE) else est_wide
  out
}

# ---- Utilities ----
#' Convenience function: full KM workflow
#'
#' Returns fit, plot, log-rank, estimates at time points, and summary table.
#'
#' @param data,time,event,group see above
#' @param time_points optional numeric vector
#' @param tau optional RMST truncation time
#' @return list
#' @export
km_full_report <- function(data, time, event, group = NULL, time_points = NULL, tau = NULL) {
  fit_res <- km_fit(data, time, event, group)
  plt <- tryCatch(km_plot(fit_res$fit, data = data, group = group), error = function(e) NULL)
  lr <- if (!is.null(group)) tryCatch(km_logrank(data, time, event, group), error = function(e) NULL) else NULL
  est <- km_estimates_at(fit_res$fit, time_points = if (is.null(time_points)) unique(round(stats::quantile(fit_res$fit$time, c(0.25,0.5,0.75)),1)) else time_points)
  rmst <- if (!is.null(group)) tryCatch(km_rmst(data, time, event, group, tau = tau), error = function(e) NULL) else NULL
  tab <- tryCatch(km_summary_table(fit_res$fit, data = data, group = group, time_points = time_points), error = function(e) NULL)
  list(fit = fit_res$fit, plot = plt, logrank = lr, estimates = est, rmst = rmst, table = tab)
}
