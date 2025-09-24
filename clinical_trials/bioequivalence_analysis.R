# =============================================================================
# Bioequivalence Analysis (2x2 Crossover) - Educational Script
# Compatible with data/clinical_trials/bioequivalence_data.csv
# =============================================================================
# Purpose:
# - Read crossover bioequivalence dataset
# - Compute log-transformed PK parameters (Cmax, AUC)
# - Fit linear model with sequence, subject(sequence), period, treatment
# - Obtain geometric LSMeans ratio (Test/Reference)
# - Construct 90% CI and decide BE within 80-125%
# - Provide clean functions for reproducible use
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(broom)
})

# -----------------------------------------------------------------------------
# Expected columns:
# subject_id, sequence (TR/RT), period (1/2), treatment (Reference/Test),
# Cmax, AUC, Tmax
# -----------------------------------------------------------------------------

#' Read bioequivalence dataset
#' @param path CSV path (default: data/clinical_trials/bioequivalence_data.csv)
#' @return tibble with factor variables
#' @export
read_be_data <- function(path = "data/clinical_trials/bioequivalence_data.csv") {
  dat <- readr::read_csv(path, show_col_types = FALSE)
  dat <- dat %>% mutate(
    subject_id = factor(subject_id),
    sequence   = factor(sequence),
    period     = factor(period),
    treatment  = factor(treatment, levels = c("Reference", "Test"))
  )
  dat
}

#' BE ANOVA on log scale for a PK variable
#' @param data data.frame
#' @param pk_var character: one of Cmax, AUC
#' @return list with estimates and CI
#' @export
be_anova <- function(data, pk_var = "Cmax") {
  stopifnot(pk_var %in% names(data))
  data$log_pk <- log(data[[pk_var]])
  data <- data %>% mutate(subj_seq = interaction(sequence, subject_id, drop = TRUE))
  fit <- lm(log_pk ~ sequence + subj_seq + period + treatment, data = data)

  cf <- summary(fit)$coefficients
  if (!"treatmentTest" %in% rownames(cf)) stop("Treatment contrast not found")
  est <- unname(cf["treatmentTest", "Estimate"])  # log(Test/Ref)
  se  <- unname(cf["treatmentTest", "Std. Error"])
  df  <- unname(summary(fit)$df[2])

  alpha <- 0.10
  tcrit <- qt(1 - alpha/2, df)
  ci_log <- c(est - tcrit * se, est + tcrit * se)

  list(
    model = fit,
    pk_var = pk_var,
    estimate_log = est,
    se_log = se,
    df = df,
    gmr = exp(est),
    ci90_lower = exp(ci_log[1]),
    ci90_upper = exp(ci_log[2]),
    be_within_80_125 = (exp(ci_log[1]) >= 0.80) && (exp(ci_log[2]) <= 1.25)
  )
}

#' Two One-Sided Tests (TOST)
#' @export
be_tost <- function(est_log, se_log, df, lower = log(0.8), upper = log(1.25)) {
  t1 <- (est_log - lower) / se_log
  p1 <- 1 - pt(t1, df)
  t2 <- (upper - est_log) / se_log
  p2 <- 1 - pt(t2, df)
  list(t_lower = t1, p_lower = p1, t_upper = t2, p_upper = p2,
       decision = (p1 < 0.05) && (p2 < 0.05))
}

#' Run BE for both Cmax and AUC (if available)
#' @export
run_be_all <- function(data) {
  res_cmax <- be_anova(data, "Cmax")
  res_auc  <- if ("AUC" %in% names(data)) be_anova(data, "AUC") else NULL
  list(Cmax = res_cmax, AUC = res_auc)
}

#' Summary table
#' @export
be_summary_table <- function(res) {
  tibble::tibble(
    PK = res$pk_var,
    GMR = round(res$gmr, 3),
    CI90_L = round(res$ci90_lower, 3),
    CI90_U = round(res$ci90_upper, 3),
    BE_80_125 = ifelse(res$be_within_80_125, "Yes", "No")
  )
}

# Example
# dat <- read_be_data()
# out <- run_be_all(dat)
# be_summary_table(out$Cmax)
# if (!is.null(out$AUC)) be_summary_table(out$AUC)

# =============================================================================
# End of file
# =============================================================================
