#' Clinical Trials Analysis Functions
#'
#' @description A collection of functions for clinical trials analysis.
#'
#' @importFrom stats power.t.test power.prop.test t.test prop.test fisher.test
#' @importFrom stats lm glm binomial quasibinomial poisson quasipoisson
#' @importFrom stats anova aov TukeyHSD
#' @importFrom dplyr %>% filter select mutate group_by summarise
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_minimal labs
#' @importFrom tableone CreateTableOne
#'
#' @name clinical_trials
NULL

#' Calculate sample size for comparing two means
#'
#' @param delta Expected difference between means
#' @param sd Standard deviation (assumed equal in both groups)
#' @param power Desired power (default: 0.8)
#' @param sig_level Significance level (default: 0.05)
#' @param alternative Alternative hypothesis ("two.sided", "less", or "greater")
#' @param ratio Allocation ratio (n2/n1) (default: 1)
#'
#' @return A list with sample size calculation results
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate sample size for detecting a difference of 10 units with SD of 20
#' ss_means <- sample_size_means(delta = 10, sd = 20)
#' }
sample_size_means <- function(delta, sd, power = 0.8, sig_level = 0.05, 
                             alternative = "two.sided", ratio = 1) {
  # Calculate sample size
  result <- stats::power.t.test(
    delta = delta,
    sd = sd,
    power = power,
    sig.level = sig_level,
    type = "two.sample",
    alternative = alternative,
    strict = TRUE
  )
  
  # Adjust for unequal allocation
  if (ratio != 1) {
    n1 <- result$n * (1 + ratio) / (2 * ratio)
    n2 <- n1 * ratio
    result$n <- NULL
    result$n1 <- ceiling(n1)
    result$n2 <- ceiling(n2)
    result$total_n <- ceiling(n1 + n2)
  } else {
    result$n <- ceiling(result$n)
    result$total_n <- ceiling(2 * result$n)
  }
  
  return(result)
}

#' Calculate sample size for comparing two proportions
#'
#' @param p1 Expected proportion in group 1
#' @param p2 Expected proportion in group 2
#' @param power Desired power (default: 0.8)
#' @param sig_level Significance level (default: 0.05)
#' @param alternative Alternative hypothesis ("two.sided", "less", or "greater")
#' @param ratio Allocation ratio (n2/n1) (default: 1)
#'
#' @return A list with sample size calculation results
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate sample size for detecting a difference between 0.3 and 0.5
#' ss_props <- sample_size_proportions(p1 = 0.3, p2 = 0.5)
#' }
sample_size_proportions <- function(p1, p2, power = 0.8, sig_level = 0.05, 
                                   alternative = "two.sided", ratio = 1) {
  # Calculate sample size
  result <- stats::power.prop.test(
    p1 = p1,
    p2 = p2,
    power = power,
    sig.level = sig_level,
    alternative = alternative,
    strict = TRUE
  )
  
  # Adjust for unequal allocation
  if (ratio != 1) {
    n1 <- result$n * (1 + ratio) / (2 * ratio)
    n2 <- n1 * ratio
    result$n <- NULL
    result$n1 <- ceiling(n1)
    result$n2 <- ceiling(n2)
    result$total_n <- ceiling(n1 + n2)
  } else {
    result$n <- ceiling(result$n)
    result$total_n <- ceiling(2 * result$n)
  }
  
  return(result)
}

#' Create baseline characteristics table
#'
#' @param data A data frame containing the data
#' @param vars Character vector of variable names to include
#' @param group Name of the grouping variable (optional)
#' @param cat_vars Character vector of categorical variable names
#' @param test Whether to include p-values for group comparisons
#' @param digits Number of digits for rounding
#'
#' @return A tableone object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create baseline characteristics table by sex
#' baseline_table <- create_baseline_table(lung, 
#'                                        vars = c("age", "ph.ecog", "ph.karno"),
#'                                        group = "sex",
#'                                        cat_vars = "ph.ecog")
#' }
create_baseline_table <- function(data, vars, group = NULL, cat_vars = NULL, test = TRUE, digits = 2) {
  # Check if required columns exist
  for (var in vars) {
    if (!var %in% names(data)) {
      stop("Variable not found in data: ", var)
    }
  }
  if (!is.null(group) && !group %in% names(data)) {
    stop("Group variable not found in data: ", group)
  }
  
  # Create TableOne object
  table_one <- tableone::CreateTableOne(
    vars = vars,
    data = data,
    factorVars = cat_vars,
    strata = group,
    test = test
  )
  
  # Print formatted table
  formatted_table <- print(table_one, printToggle = FALSE, nonnormal = vars, digits = digits)
  
  # Return both the TableOne object and the formatted table
  return(list(table_one = table_one, formatted_table = formatted_table))
}

#' Perform t-test for continuous outcomes
#'
#' @param data A data frame containing the data
#' @param outcome Name of the outcome variable
#' @param group Name of the grouping variable
#' @param paired Whether the test is paired
#' @param var_equal Whether to assume equal variances
#' @param conf_level Confidence level for intervals (default: 0.95)
#'
#' @return A t.test object
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform t-test for age by sex
#' t_test_result <- perform_t_test(lung, "age", "sex")
#' }
perform_t_test <- function(data, outcome, group, paired = FALSE, var_equal = FALSE, conf_level = 0.95) {
  # Check if required columns exist
  if (!outcome %in% names(data)) {
    stop("Outcome variable not found in data: ", outcome)
  }
  if (!group %in% names(data)) {
    stop("Group variable not found in data: ", group)
  }
  
  # Extract outcome and group variables
  outcome_values <- data[[outcome]]
  group_values <- data[[group]]
  
  # Check if group variable has exactly two levels for unpaired test
  if (!paired) {
    group_levels <- unique(group_values)
    if (length(group_levels) != 2) {
      stop("Group variable must have exactly two levels for unpaired t-test")
    }
    
    # Split outcome by group
    group1_values <- outcome_values[group_values == group_levels[1]]
    group2_values <- outcome_values[group_values == group_levels[2]]
    
    # Perform t-test
    result <- stats::t.test(
      x = group1_values,
      y = group2_values,
      paired = paired,
      var.equal = var_equal,
      conf.level = conf_level
    )
  } else {
    # For paired test, create formula
    t_formula <- stats::as.formula(paste(outcome, "~", group))
    
    # Perform paired t-test
    result <- stats::t.test(
      formula = t_formula,
      data = data,
      paired = paired,
      conf.level = conf_level
    )
  }
  
  return(result)
}

#' Perform chi-square or Fisher's exact test for categorical outcomes
#'
#' @param data A data frame containing the data
#' @param outcome Name of the outcome variable
#' @param group Name of the grouping variable
#' @param fisher Whether to use Fisher's exact test
#' @param conf_level Confidence level for intervals (default: 0.95)
#'
#' @return A prop.test or fisher.test object
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform chi-square test for status by sex
#' chi_sq_result <- perform_categorical_test(lung, "status", "sex")
#'
#' # Perform Fisher's exact test
#' fisher_result <- perform_categorical_test(lung, "status", "sex", fisher = TRUE)
#' }
perform_categorical_test <- function(data, outcome, group, fisher = FALSE, conf_level = 0.95) {
  # Check if required columns exist
  if (!outcome %in% names(data)) {
    stop("Outcome variable not found in data: ", outcome)
  }
  if (!group %in% names(data)) {
    stop("Group variable not found in data: ", group)
  }
  
  # Create contingency table
  cont_table <- table(data[[group]], data[[outcome]])
  
  # Perform test
  if (fisher) {
    result <- stats::fisher.test(cont_table, conf.level = conf_level)
  } else {
    result <- stats::prop.test(cont_table, conf.level = conf_level)
  }
  
  return(result)
}

#' Perform ANOVA for continuous outcomes with multiple groups
#'
#' @param data A data frame containing the data
#' @param outcome Name of the outcome variable
#' @param group Name of the grouping variable
#' @param posthoc Whether to perform post-hoc tests
#'
#' @return A list with ANOVA results and post-hoc tests (if requested)
#' @export
#'
#' @examples
#' \dontrun{
#' # Create sample data with multiple groups
#' data <- data.frame(
#'   outcome = rnorm(100),
#'   group = sample(c("A", "B", "C"), 100, replace = TRUE)
#' )
#'
#' # Perform ANOVA with post-hoc tests
#' anova_result <- perform_anova(data, "outcome", "group", posthoc = TRUE)
#' }
perform_anova <- function(data, outcome, group, posthoc = FALSE) {
  # Check if required columns exist
  if (!outcome %in% names(data)) {
    stop("Outcome variable not found in data: ", outcome)
  }
  if (!group %in% names(data)) {
    stop("Group variable not found in data: ", group)
  }
  
  # Create formula for ANOVA
  anova_formula <- stats::as.formula(paste(outcome, "~", group))
  
  # Perform ANOVA
  anova_model <- stats::aov(anova_formula, data = data)
  anova_result <- stats::anova(anova_model)
  
  # Perform post-hoc tests if requested
  if (posthoc) {
    posthoc_result <- stats::TukeyHSD(anova_model)
    return(list(anova = anova_result, posthoc = posthoc_result))
  } else {
    return(list(anova = anova_result))
  }
}

#' Fit linear regression model
#'
#' @param data A data frame containing the data
#' @param outcome Name of the outcome variable
#' @param covariates Character vector of covariate column names
#' @param robust Whether to use robust standard errors
#'
#' @return A lm object
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit linear regression model for age with sex and ph.ecog as covariates
#' lm_model <- fit_linear_model(lung, "age", c("sex", "ph.ecog"))
#' }
fit_linear_model <- function(data, outcome, covariates, robust = FALSE) {
  # Check if required columns exist
  if (!outcome %in% names(data)) {
    stop("Outcome variable not found in data: ", outcome)
  }
  for (cov in covariates) {
    if (!cov %in% names(data)) {
      stop("Covariate not found in data: ", cov)
    }
  }
  
  # Create formula for linear model
  lm_formula <- stats::as.formula(paste(outcome, "~", paste(covariates, collapse = " + ")))
  
  # Fit linear model
  lm_model <- stats::lm(lm_formula, data = data)
  
  # Apply robust standard errors if requested
  if (robust) {
    if (!requireNamespace("sandwich", quietly = TRUE)) {
      stop("Package 'sandwich' is required for robust standard errors")
    }
    
    # Calculate robust standard errors
    robust_se <- sandwich::vcovHC(lm_model, type = "HC1")
    
    # Update model with robust standard errors
    lm_model$robust_se <- sqrt(diag(robust_se))
    lm_model$is_robust <- TRUE
  }
  
  return(lm_model)
}

#' Fit generalized linear model
#'
#' @param data A data frame containing the data
#' @param outcome Name of the outcome variable
#' @param covariates Character vector of covariate column names
#' @param family Model family (e.g., "binomial", "poisson")
#' @param link Link function (default depends on family)
#' @param robust Whether to use robust standard errors
#'
#' @return A glm object
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit logistic regression model for status with age and sex as covariates
#' glm_model <- fit_glm(lung, "status", c("age", "sex"), family = "binomial")
#' }
fit_glm <- function(data, outcome, covariates, family, link = NULL, robust = FALSE) {
  # Check if required columns exist
  if (!outcome %in% names(data)) {
    stop("Outcome variable not found in data: ", outcome)
  }
  for (cov in covariates) {
    if (!cov %in% names(data)) {
      stop("Covariate not found in data: ", cov)
    }
  }
  
  # Create formula for GLM
  glm_formula <- stats::as.formula(paste(outcome, "~", paste(covariates, collapse = " + ")))
  
  # Determine family and link
  if (family == "binomial") {
    if (is.null(link)) link <- "logit"
    fam <- stats::binomial(link = link)
  } else if (family == "poisson") {
    if (is.null(link)) link <- "log"
    fam <- stats::poisson(link = link)
  } else if (family == "quasibinomial") {
    if (is.null(link)) link <- "logit"
    fam <- stats::quasibinomial(link = link)
  } else if (family == "quasipoisson") {
    if (is.null(link)) link <- "log"
    fam <- stats::quasipoisson(link = link)
  } else if (family == "gaussian") {
    if (is.null(link)) link <- "identity"
    fam <- stats::gaussian(link = link)
  } else {
    stop("Unsupported family: ", family)
  }
  
  # Fit GLM
  glm_model <- stats::glm(glm_formula, family = fam, data = data)
  
  # Apply robust standard errors if requested
  if (robust) {
    if (!requireNamespace("sandwich", quietly = TRUE)) {
      stop("Package 'sandwich' is required for robust standard errors")
    }
    
    # Calculate robust standard errors
    robust_se <- sandwich::vcovHC(glm_model, type = "HC1")
    
    # Update model with robust standard errors
    glm_model$robust_se <- sqrt(diag(robust_se))
    glm_model$is_robust <- TRUE
  }
  
  return(glm_model)
}

#' Calculate odds ratio from logistic regression
#'
#' @param glm_model A glm object from fit_glm with family = "binomial"
#' @param conf_level Confidence level for intervals (default: 0.95)
#'
#' @return A data frame with odds ratios and confidence intervals
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit logistic regression model
#' glm_model <- fit_glm(lung, "status", c("age", "sex"), family = "binomial")
#'
#' # Calculate odds ratios
#' odds_ratios <- calculate_odds_ratio(glm_model)
#' }
calculate_odds_ratio <- function(glm_model, conf_level = 0.95) {
  # Check if model is a logistic regression
  if (!inherits(glm_model, "glm") || !identical(glm_model$family$family, "binomial")) {
    stop("Model must be a logistic regression (glm with binomial family)")
  }
  
  # Calculate odds ratios and confidence intervals
  coefs <- stats::coef(glm_model)
  
  # Use robust standard errors if available
  if (!is.null(glm_model$robust_se)) {
    se <- glm_model$robust_se
  } else {
    se <- sqrt(diag(stats::vcov(glm_model)))
  }
  
  # Calculate z-value for confidence interval
  z_value <- stats::qnorm(1 - (1 - conf_level) / 2)
  
  # Calculate odds ratios and confidence intervals
  or <- exp(coefs)
  lower_ci <- exp(coefs - z_value * se)
  upper_ci <- exp(coefs + z_value * se)
  
  # Create data frame with results
  result <- data.frame(
    term = names(coefs),
    odds_ratio = or,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    p_value = 2 * (1 - stats::pnorm(abs(coefs / se)))
  )
  
  return(result)
}

#' Calculate relative risk from log-binomial model
#'
#' @param glm_model A glm object from fit_glm with family = "binomial" and link = "log"
#' @param conf_level Confidence level for intervals (default: 0.95)
#'
#' @return A data frame with relative risks and confidence intervals
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit log-binomial model
#' glm_model <- fit_glm(lung, "status", c("age", "sex"), family = "binomial", link = "log")
#'
#' # Calculate relative risks
#' relative_risks <- calculate_relative_risk(glm_model)
#' }
calculate_relative_risk <- function(glm_model, conf_level = 0.95) {
  # Check if model is a log-binomial regression
  if (!inherits(glm_model, "glm") || 
      !identical(glm_model$family$family, "binomial") || 
      !identical(glm_model$family$link, "log")) {
    stop("Model must be a log-binomial regression (glm with binomial family and log link)")
  }
  
  # Calculate relative risks and confidence intervals
  coefs <- stats::coef(glm_model)
  
  # Use robust standard errors if available
  if (!is.null(glm_model$robust_se)) {
    se <- glm_model$robust_se
  } else {
    se <- sqrt(diag(stats::vcov(glm_model)))
  }
  
  # Calculate z-value for confidence interval
  z_value <- stats::qnorm(1 - (1 - conf_level) / 2)
  
  # Calculate relative risks and confidence intervals
  rr <- exp(coefs)
  lower_ci <- exp(coefs - z_value * se)
  upper_ci <- exp(coefs + z_value * se)
  
  # Create data frame with results
  result <- data.frame(
    term = names(coefs),
    relative_risk = rr,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    p_value = 2 * (1 - stats::pnorm(abs(coefs / se)))
  )
  
  return(result)
}

#' Perform interim analysis with alpha spending
#'
#' @param observed_z Z-statistic from current analysis
#' @param information_fraction Fraction of information accrued (0 to 1)
#' @param alpha Overall significance level (default: 0.05)
#' @param spending_function Alpha spending function ("obrien_fleming", "pocock", "power")
#' @param power_param Power parameter for power spending function (default: 3)
#'
#' @return A list with interim analysis results
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform interim analysis with O'Brien-Fleming spending function
#' interim_result <- interim_analysis(observed_z = 2.5, information_fraction = 0.5)
#' }
interim_analysis <- function(observed_z, information_fraction, alpha = 0.05, 
                            spending_function = "obrien_fleming", power_param = 3) {
  # Check if information fraction is valid
  if (information_fraction <= 0 || information_fraction > 1) {
    stop("Information fraction must be between 0 and 1")
  }
  
  # Calculate critical value based on spending function
  if (spending_function == "obrien_fleming") {
    # O'Brien-Fleming spending function
    critical_z <- stats::qnorm(1 - alpha / 2) / sqrt(information_fraction)
  } else if (spending_function == "pocock") {
    # Pocock spending function
    critical_z <- stats::qnorm(1 - alpha / (2 * information_fraction))
  } else if (spending_function == "power") {
    # Power spending function
    spent_alpha <- alpha * information_fraction^power_param
    critical_z <- stats::qnorm(1 - spent_alpha / 2)
  } else {
    stop("Unsupported spending function: ", spending_function)
  }
  
  # Calculate p-value
  p_value <- 2 * (1 - stats::pnorm(abs(observed_z)))
  
  # Determine decision
  decision <- ifelse(abs(observed_z) >= critical_z, "Reject H0", "Continue")
  
  # Create result
  result <- list(
    observed_z = observed_z,
    critical_z = critical_z,
    p_value = p_value,
    information_fraction = information_fraction,
    spending_function = spending_function,
    decision = decision
  )
  
  return(result)
}

#' Calculate number needed to treat (NNT)
#'
#' @param p_treatment Event rate in treatment group
#' @param p_control Event rate in control group
#' @param conf_level Confidence level for intervals (default: 0.95)
#' @param n_treatment Sample size in treatment group
#' @param n_control Sample size in control group
#'
#' @return A list with NNT and confidence interval
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate NNT with 95% confidence interval
#' nnt_result <- calculate_nnt(p_treatment = 0.15, p_control = 0.25, 
#'                            n_treatment = 100, n_control = 100)
#' }
calculate_nnt <- function(p_treatment, p_control, conf_level = 0.95, 
                         n_treatment = NULL, n_control = NULL) {
  # Calculate absolute risk reduction
  arr <- p_control - p_treatment
  
  # Calculate NNT
  nnt <- 1 / arr
  
  # Calculate confidence interval if sample sizes are provided
  if (!is.null(n_treatment) && !is.null(n_control)) {
    # Standard error of risk difference
    se_arr <- sqrt(p_treatment * (1 - p_treatment) / n_treatment + 
                  p_control * (1 - p_control) / n_control)
    
    # Z-value for confidence interval
    z_value <- stats::qnorm(1 - (1 - conf_level) / 2)
    
    # Confidence interval for ARR
    arr_lower <- arr - z_value * se_arr
    arr_upper <- arr + z_value * se_arr
    
    # Confidence interval for NNT
    # Note: NNT CI is the reciprocal of ARR CI, with bounds swapped
    if (arr_lower * arr_upper > 0) {
      # Both bounds have the same sign
      nnt_lower <- 1 / arr_upper
      nnt_upper <- 1 / arr_lower
    } else {
      # Bounds have different signs, indicating ARR includes zero
      nnt_lower <- -Inf
      nnt_upper <- Inf
    }
    
    result <- list(
      nnt = nnt,
      nnt_lower = nnt_lower,
      nnt_upper = nnt_upper,
      arr = arr,
      arr_lower = arr_lower,
      arr_upper = arr_upper
    )
  } else {
    # Return only point estimate if sample sizes are not provided
    result <- list(
      nnt = nnt,
      arr = arr
    )
  }
  
  return(result)
}

