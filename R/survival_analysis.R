#' Survival Analysis Functions
#'
#' @description A collection of functions for survival analysis in clinical research.
#'
#' @importFrom survival Surv survfit coxph cox.zph survdiff
#' @importFrom survminer ggsurvplot ggcoxzph ggforest
#' @importFrom cmprsk crr cuminc
#' @importFrom dplyr %>% filter select mutate group_by summarise
#' @importFrom ggplot2 ggplot aes geom_line theme_minimal labs
#' @importFrom stats pchisq
#'
#' @name survival_analysis
NULL

#' Fit Kaplan-Meier survival curves
#'
#' @param data A data frame containing the survival data
#' @param time Name of the column containing time-to-event data
#' @param event Name of the column containing event indicator (1=event, 0=censored)
#' @param group Name of the column containing grouping variable (optional)
#' @param conf_int Whether to include confidence intervals
#' @param conf_level Confidence level for intervals (default: 0.95)
#'
#' @return A survfit object
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit Kaplan-Meier curve
#' km_fit <- fit_kaplan_meier(lung, "time", "status")
#'
#' # Fit Kaplan-Meier curves by sex
#' km_fit_by_sex <- fit_kaplan_meier(lung, "time", "status", "sex")
#' }
fit_kaplan_meier <- function(data, time, event, group = NULL, conf_int = TRUE, conf_level = 0.95) {
  # Check if required columns exist
  if (!time %in% names(data)) {
    stop("Time column not found in data: ", time)
  }
  if (!event %in% names(data)) {
    stop("Event column not found in data: ", event)
  }
  if (!is.null(group) && !group %in% names(data)) {
    stop("Group column not found in data: ", group)
  }
  
  # Create survival formula
  if (is.null(group)) {
    surv_formula <- stats::as.formula(paste("survival::Surv(", time, ",", event, ") ~ 1"))
  } else {
    surv_formula <- stats::as.formula(paste("survival::Surv(", time, ",", event, ") ~", group))
  }
  
  # Fit Kaplan-Meier model
  km_fit <- survival::survfit(surv_formula, data = data, conf.type = "log-log", conf.int = conf_level)
  
  return(km_fit)
}

#' Plot Kaplan-Meier survival curves
#'
#' @param fit A survfit object from fit_kaplan_meier
#' @param data The original data used to fit the model
#' @param time Name of the column containing time-to-event data
#' @param event Name of the column containing event indicator
#' @param group Name of the column containing grouping variable (optional)
#' @param risk_table Whether to include risk table
#' @param conf_int Whether to include confidence intervals
#' @param palette Color palette for groups
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param legend_title Legend title
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit Kaplan-Meier curves by sex
#' km_fit <- fit_kaplan_meier(lung, "time", "status", "sex")
#'
#' # Plot survival curves
#' plot_kaplan_meier(km_fit, lung, "time", "status", "sex", risk_table = TRUE)
#' }
plot_kaplan_meier <- function(fit, data, time, event, group = NULL, risk_table = FALSE,
                             conf_int = TRUE, palette = "jco", title = "Kaplan-Meier Survival Curve",
                             xlab = "Time", ylab = "Survival Probability", legend_title = "Group") {
  # Create survival plot
  p <- survminer::ggsurvplot(
    fit = fit,
    data = data,
    risk.table = risk_table,
    pval = !is.null(group),
    conf.int = conf_int,
    palette = palette,
    xlab = xlab,
    ylab = ylab,
    title = title,
    legend.title = legend_title,
    risk.table.height = 0.25,
    ggtheme = ggplot2::theme_minimal()
  )
  
  return(p)
}

#' Fit Cox proportional hazards model
#'
#' @param data A data frame containing the survival data
#' @param time Name of the column containing time-to-event data
#' @param event Name of the column containing event indicator (1=event, 0=censored)
#' @param covariates Character vector of covariate column names
#' @param strata Name of the column for stratification (optional)
#' @param robust Whether to use robust standard errors
#'
#' @return A coxph object
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit Cox model with age and sex as covariates
#' cox_model <- fit_cox_model(lung, "time", "status", c("age", "sex"))
#'
#' # Fit stratified Cox model
#' cox_model_strat <- fit_cox_model(lung, "time", "status", c("age"), strata = "sex")
#' }
fit_cox_model <- function(data, time, event, covariates, strata = NULL, robust = FALSE) {
  # Check if required columns exist
  if (!time %in% names(data)) {
    stop("Time column not found in data: ", time)
  }
  if (!event %in% names(data)) {
    stop("Event column not found in data: ", event)
  }
  for (cov in covariates) {
    if (!cov %in% names(data)) {
      stop("Covariate not found in data: ", cov)
    }
  }
  if (!is.null(strata) && !strata %in% names(data)) {
    stop("Stratification variable not found in data: ", strata)
  }
  
  # Create formula for Cox model
  if (is.null(strata)) {
    cox_formula <- stats::as.formula(paste("survival::Surv(", time, ",", event, ") ~", 
                                          paste(covariates, collapse = " + ")))
  } else {
    cox_formula <- stats::as.formula(paste("survival::Surv(", time, ",", event, ") ~", 
                                          paste(covariates, collapse = " + "), 
                                          " + strata(", strata, ")"))
  }
  
  # Fit Cox model
  cox_fit <- survival::coxph(cox_formula, data = data, robust = robust)
  
  return(cox_fit)
}

#' Test proportional hazards assumption
#'
#' @param cox_model A coxph object from fit_cox_model
#' @param plot Whether to create diagnostic plots
#'
#' @return A cox.zph object with test results
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit Cox model
#' cox_model <- fit_cox_model(lung, "time", "status", c("age", "sex"))
#'
#' # Test proportional hazards assumption
#' ph_test <- test_proportional_hazards(cox_model, plot = TRUE)
#' }
test_proportional_hazards <- function(cox_model, plot = FALSE) {
  # Test proportional hazards assumption
  ph_test <- survival::cox.zph(cox_model)
  
  # Create diagnostic plots if requested
  if (plot) {
    p <- survminer::ggcoxzph(ph_test)
    print(p)
  }
  
  return(ph_test)
}

#' Plot forest plot of hazard ratios
#'
#' @param cox_model A coxph object from fit_cox_model
#' @param data The original data used to fit the model
#' @param covariates Character vector of covariate names to include (optional)
#' @param conf_level Confidence level for intervals (default: 0.95)
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit Cox model
#' cox_model <- fit_cox_model(lung, "time", "status", c("age", "sex", "ph.ecog"))
#'
#' # Create forest plot
#' forest_plot <- plot_forest(cox_model, lung)
#' }
plot_forest <- function(cox_model, data, covariates = NULL, conf_level = 0.95) {
  # Create forest plot
  p <- survminer::ggforest(
    model = cox_model,
    data = data,
    main = "Hazard Ratios",
    cpositions = c(0.02, 0.22, 0.4),
    fontsize = 0.8,
    refLabel = "reference",
    noDigits = 2
  )
  
  return(p)
}

#' Perform log-rank test for survival curves
#'
#' @param data A data frame containing the survival data
#' @param time Name of the column containing time-to-event data
#' @param event Name of the column containing event indicator (1=event, 0=censored)
#' @param group Name of the column containing grouping variable
#'
#' @return A survdiff object with test results
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform log-rank test by sex
#' lr_test <- logrank_test(lung, "time", "status", "sex")
#' }
logrank_test <- function(data, time, event, group) {
  # Check if required columns exist
  if (!time %in% names(data)) {
    stop("Time column not found in data: ", time)
  }
  if (!event %in% names(data)) {
    stop("Event column not found in data: ", event)
  }
  if (!group %in% names(data)) {
    stop("Group column not found in data: ", group)
  }
  
  # Create formula for log-rank test
  lr_formula <- stats::as.formula(paste("survival::Surv(", time, ",", event, ") ~", group))
  
  # Perform log-rank test
  lr_test <- survival::survdiff(lr_formula, data = data)
  
  # Calculate p-value
  p_value <- 1 - stats::pchisq(lr_test$chisq, length(lr_test$n) - 1)
  
  # Add p-value to result
  lr_test$p_value <- p_value
  
  return(lr_test)
}

#' Fit competing risks model
#'
#' @param data A data frame containing the survival data
#' @param time Name of the column containing time-to-event data
#' @param event Name of the column containing event indicator (0=censored, 1,2,...=event types)
#' @param event_type Event type of interest
#' @param covariates Character vector of covariate column names
#'
#' @return A crr object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create sample competing risks data
#' cr_data <- data.frame(
#'   time = c(5, 10, 15, 20, 25),
#'   event = c(1, 0, 2, 1, 0),
#'   age = c(65, 70, 55, 60, 75),
#'   sex = c(1, 0, 1, 0, 1)
#' )
#'
#' # Fit competing risks model for event type 1
#' cr_model <- fit_competing_risks(cr_data, "time", "event", 1, c("age", "sex"))
#' }
fit_competing_risks <- function(data, time, event, event_type, covariates) {
  # Check if required columns exist
  if (!time %in% names(data)) {
    stop("Time column not found in data: ", time)
  }
  if (!event %in% names(data)) {
    stop("Event column not found in data: ", event)
  }
  for (cov in covariates) {
    if (!cov %in% names(data)) {
      stop("Covariate not found in data: ", cov)
    }
  }
  
  # Extract time, event, and covariates
  time_values <- data[[time]]
  event_values <- data[[event]]
  cov_matrix <- as.matrix(data[, covariates, drop = FALSE])
  
  # Fit competing risks model
  cr_model <- cmprsk::crr(time_values, event_values, cov_matrix, failcode = event_type)
  
  return(cr_model)
}

#' Calculate cumulative incidence for competing risks
#'
#' @param data A data frame containing the survival data
#' @param time Name of the column containing time-to-event data
#' @param event Name of the column containing event indicator (0=censored, 1,2,...=event types)
#' @param group Name of the column containing grouping variable (optional)
#'
#' @return A cuminc object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create sample competing risks data
#' cr_data <- data.frame(
#'   time = c(5, 10, 15, 20, 25),
#'   event = c(1, 0, 2, 1, 0),
#'   group = c("A", "B", "A", "B", "A")
#' )
#'
#' # Calculate cumulative incidence by group
#' ci <- cumulative_incidence(cr_data, "time", "event", "group")
#' }
cumulative_incidence <- function(data, time, event, group = NULL) {
  # Check if required columns exist
  if (!time %in% names(data)) {
    stop("Time column not found in data: ", time)
  }
  if (!event %in% names(data)) {
    stop("Event column not found in data: ", event)
  }
  if (!is.null(group) && !group %in% names(data)) {
    stop("Group column not found in data: ", group)
  }
  
  # Extract time and event
  time_values <- data[[time]]
  event_values <- data[[event]]
  
  # Calculate cumulative incidence
  if (is.null(group)) {
    ci <- cmprsk::cuminc(time_values, event_values)
  } else {
    group_values <- data[[group]]
    ci <- cmprsk::cuminc(time_values, event_values, group_values)
  }
  
  return(ci)
}

#' Plot cumulative incidence curves
#'
#' @param ci A cuminc object from cumulative_incidence
#' @param event_types Event types to include in the plot
#' @param palette Color palette for event types
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Y-axis label
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate cumulative incidence
#' ci <- cumulative_incidence(cr_data, "time", "event", "group")
#'
#' # Plot cumulative incidence curves
#' plot_cumulative_incidence(ci, c(1, 2))
#' }
plot_cumulative_incidence <- function(ci, event_types, palette = "jco",
                                     title = "Cumulative Incidence", 
                                     xlab = "Time", ylab = "Cumulative Incidence") {
  # Extract data from cuminc object
  ci_data <- list()
  for (i in seq_along(ci$est)) {
    name <- names(ci$est)[i]
    event_type <- as.numeric(sub(".*event=([0-9]+).*", "\\1", name))
    group <- sub(".*group=([^,]+).*", "\\1", name)
    
    if (event_type %in% event_types) {
      ci_data[[name]] <- data.frame(
        time = ci$time[[i]],
        est = ci$est[[i]],
        var = ci$var[[i]],
        event_type = event_type,
        group = group
      )
    }
  }
  
  # Combine data
  ci_df <- do.call(rbind, ci_data)
  
  # Create plot
  p <- ggplot2::ggplot(ci_df, ggplot2::aes(x = time, y = est, color = group, linetype = factor(event_type))) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab,
      color = "Group",
      linetype = "Event Type"
    )
  
  return(p)
}

#' Calculate restricted mean survival time (RMST)
#'
#' @param fit A survfit object from fit_kaplan_meier
#' @param tau Time horizon for RMST calculation
#' @param conf_level Confidence level for intervals (default: 0.95)
#'
#' @return A data frame with RMST estimates and confidence intervals
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit Kaplan-Meier curves by sex
#' km_fit <- fit_kaplan_meier(lung, "time", "status", "sex")
#'
#' # Calculate RMST with 365-day horizon
#' rmst <- calculate_rmst(km_fit, 365)
#' }
calculate_rmst <- function(fit, tau, conf_level = 0.95) {
  # Check if tau is within the range of observed times
  max_time <- max(fit$time)
  if (tau > max_time) {
    warning("Time horizon (tau) is beyond the maximum observed time. Results may be unreliable.")
  }
  
  # Calculate RMST
  rmst_result <- survminer::surv_rmst(fit, tau = tau, conf.level = conf_level)
  
  return(rmst_result)
}

#' Calculate concordance index (C-index)
#'
#' @param cox_model A coxph object from fit_cox_model
#' @param data The original data used to fit the model
#' @param time Name of the column containing time-to-event data
#' @param event Name of the column containing event indicator
#'
#' @return A numeric value representing the C-index
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit Cox model
#' cox_model <- fit_cox_model(lung, "time", "status", c("age", "sex"))
#'
#' # Calculate C-index
#' c_index <- calculate_c_index(cox_model, lung, "time", "status")
#' }
calculate_c_index <- function(cox_model, data, time, event) {
  # Check if required columns exist
  if (!time %in% names(data)) {
    stop("Time column not found in data: ", time)
  }
  if (!event %in% names(data)) {
    stop("Event column not found in data: ", event)
  }
  
  # Extract C-index from model
  c_index <- cox_model$concordance["concordance"]
  
  return(c_index)
}

#' Perform time-dependent ROC analysis
#'
#' @param data A data frame containing the survival data
#' @param time Name of the column containing time-to-event data
#' @param event Name of the column containing event indicator
#' @param marker Name of the column containing the marker or risk score
#' @param pred_times Vector of time points for prediction
#' @param span Smoothing parameter for ROC curve (default: 0.25)
#'
#' @return A list with time-dependent ROC results
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit Cox model
#' cox_model <- fit_cox_model(lung, "time", "status", c("age", "sex"))
#'
#' # Calculate risk scores
#' lung$risk_score <- predict(cox_model, type = "risk")
#'
#' # Perform time-dependent ROC analysis
#' td_roc <- time_dependent_roc(lung, "time", "status", "risk_score", c(365, 730))
#' }
time_dependent_roc <- function(data, time, event, marker, pred_times, span = 0.25) {
  # Check if required columns exist
  if (!time %in% names(data)) {
    stop("Time column not found in data: ", time)
  }
  if (!event %in% names(data)) {
    stop("Event column not found in data: ", event)
  }
  if (!marker %in% names(data)) {
    stop("Marker column not found in data: ", marker)
  }
  
  # Load required package
  if (!requireNamespace("timeROC", quietly = TRUE)) {
    stop("Package 'timeROC' is required for time-dependent ROC analysis")
  }
  
  # Perform time-dependent ROC analysis
  td_roc <- timeROC::timeROC(
    T = data[[time]],
    delta = data[[event]],
    marker = data[[marker]],
    cause = 1,
    times = pred_times,
    iid = TRUE
  )
  
  return(td_roc)
}

#' Plot time-dependent ROC curves
#'
#' @param td_roc A time-dependent ROC object from time_dependent_roc
#' @param title Plot title
#' @param legend_pos Legend position
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform time-dependent ROC analysis
#' td_roc <- time_dependent_roc(lung, "time", "status", "risk_score", c(365, 730))
#'
#' # Plot ROC curves
#' plot_time_dependent_roc(td_roc)
#' }
plot_time_dependent_roc <- function(td_roc, title = "Time-Dependent ROC Curves", legend_pos = "bottomright") {
  # Load required package
  if (!requireNamespace("timeROC", quietly = TRUE)) {
    stop("Package 'timeROC' is required for plotting time-dependent ROC curves")
  }
  
  # Plot ROC curves
  plot(td_roc, time = td_roc$times, title = title, legend = legend_pos)
  
  # Return invisibly
  invisible(td_roc)
}

