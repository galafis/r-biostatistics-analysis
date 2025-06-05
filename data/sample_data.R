#' Sample datasets for biostatistics analysis
#'
#' This file contains sample datasets for demonstrating and testing
#' the biostatistics analysis functions.

#' @name lung_cancer_survival
#' @title Lung Cancer Survival Dataset
#' @description A dataset containing survival information for lung cancer patients
#' @format A data frame with 228 rows and 10 variables:
#' \describe{
#'   \item{id}{Patient identifier}
#'   \item{time}{Survival time in months}
#'   \item{status}{Censoring status (1=event occurred, 0=censored)}
#'   \item{age}{Patient age at diagnosis}
#'   \item{sex}{Patient sex (M=male, F=female)}
#'   \item{ph.ecog}{ECOG performance score (0=good, 5=dead)}
#'   \item{ph.karno}{Karnofsky performance score (0=bad, 100=good)}
#'   \item{meal.cal}{Calories consumed at meals}
#'   \item{wt.loss}{Weight loss in last six months}
#'   \item{stage}{Cancer stage (I, II, III, IV)}
#' }
#' @source Simulated data based on the lung dataset in the survival package
lung_cancer_survival <- data.frame(
  id = 1:228,
  time = round(rexp(228, 1/12) + 1),
  status = sample(0:1, 228, replace = TRUE, prob = c(0.2, 0.8)),
  age = round(rnorm(228, 65, 10)),
  sex = sample(c("M", "F"), 228, replace = TRUE, prob = c(0.6, 0.4)),
  ph.ecog = sample(0:4, 228, replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.15, 0.05)),
  ph.karno = round(rnorm(228, 70, 20)),
  meal.cal = round(rnorm(228, 1000, 300)),
  wt.loss = round(rnorm(228, 10, 8)),
  stage = sample(c("I", "II", "III", "IV"), 228, replace = TRUE, prob = c(0.15, 0.25, 0.35, 0.25))
)

# Ensure values are within reasonable ranges
lung_cancer_survival$ph.karno <- pmin(pmax(lung_cancer_survival$ph.karno, 0), 100)
lung_cancer_survival$meal.cal <- pmax(lung_cancer_survival$meal.cal, 0)
lung_cancer_survival$wt.loss <- pmax(lung_cancer_survival$wt.loss, 0)
lung_cancer_survival$age <- pmin(pmax(lung_cancer_survival$age, 30), 90)

#' @name clinical_trial_data
#' @title Clinical Trial Dataset
#' @description A dataset containing results from a randomized clinical trial
#' @format A data frame with 200 rows and 8 variables:
#' \describe{
#'   \item{patient_id}{Patient identifier}
#'   \item{treatment}{Treatment group (Treatment or Placebo)}
#'   \item{baseline}{Baseline measurement}
#'   \item{week4}{Measurement at week 4}
#'   \item{week8}{Measurement at week 8}
#'   \item{week12}{Measurement at week 12 (primary endpoint)}
#'   \item{age}{Patient age}
#'   \item{sex}{Patient sex (M=male, F=female)}
#'   \item{site}{Clinical site identifier}
#' }
#' @source Simulated data
clinical_trial_data <- data.frame(
  patient_id = 1:200,
  treatment = rep(c("Treatment", "Placebo"), each = 100),
  baseline = round(rnorm(200, 100, 15)),
  age = round(rnorm(200, 55, 12)),
  sex = sample(c("M", "F"), 200, replace = TRUE),
  site = sample(paste0("Site", 1:10), 200, replace = TRUE)
)

# Generate outcome measurements with treatment effect
set.seed(123)
effect_size <- 10
clinical_trial_data$week4 <- clinical_trial_data$baseline - 
  (ifelse(clinical_trial_data$treatment == "Treatment", effect_size * 0.4, effect_size * 0.2) + rnorm(200, 0, 5))

clinical_trial_data$week8 <- clinical_trial_data$week4 - 
  (ifelse(clinical_trial_data$treatment == "Treatment", effect_size * 0.3, effect_size * 0.1) + rnorm(200, 0, 5))

clinical_trial_data$week12 <- clinical_trial_data$week8 - 
  (ifelse(clinical_trial_data$treatment == "Treatment", effect_size * 0.3, effect_size * 0.1) + rnorm(200, 0, 5))

# Ensure values are within reasonable ranges
clinical_trial_data$baseline <- pmax(clinical_trial_data$baseline, 50)
clinical_trial_data$week4 <- pmax(clinical_trial_data$week4, 30)
clinical_trial_data$week8 <- pmax(clinical_trial_data$week8, 20)
clinical_trial_data$week12 <- pmax(clinical_trial_data$week12, 10)
clinical_trial_data$age <- pmin(pmax(clinical_trial_data$age, 18), 85)

