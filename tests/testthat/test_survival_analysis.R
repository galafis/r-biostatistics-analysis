# Test suite for survival analysis functions
library(testthat)
library(rbiostatistics)
library(survival)

context("Survival Analysis Functions")

# Sample data for testing
test_data <- data.frame(
  time = c(6, 7, 10, 15, 19, 25),
  status = c(1, 0, 1, 1, 0, 1),
  group = c("A", "A", "A", "B", "B", "B"),
  age = c(65, 72, 58, 61, 74, 69),
  sex = c("M", "F", "M", "F", "M", "F")
)

test_that("create_survival_object works correctly", {
  surv_obj <- create_survival_object(test_data, time_col = "time", status_col = "status")
  expect_s3_class(surv_obj, "Surv")
  expect_equal(length(surv_obj), nrow(test_data))
})

test_that("fit_kaplan_meier works correctly", {
  km_fit <- fit_kaplan_meier(test_data, time_col = "time", status_col = "status", group_col = "group")
  expect_s3_class(km_fit, "survfit")
  expect_equal(length(km_fit$strata), 2)  # Two groups: A and B
})

test_that("calculate_survival_metrics works correctly", {
  km_fit <- fit_kaplan_meier(test_data, time_col = "time", status_col = "status")
  metrics <- calculate_survival_metrics(km_fit)
  expect_type(metrics, "list")
  expect_true("median_survival" %in% names(metrics))
  expect_true("survival_probabilities" %in% names(metrics))
})

test_that("fit_cox_model works correctly", {
  cox_model <- fit_cox_model(test_data, time_col = "time", status_col = "status", 
                            covariates = c("group", "age", "sex"))
  expect_s3_class(cox_model, "coxph")
  expect_equal(length(cox_model$coefficients), 3)  # Three covariates
})

test_that("calculate_hazard_ratio works correctly", {
  cox_model <- fit_cox_model(test_data, time_col = "time", status_col = "status", 
                            covariates = c("group", "age"))
  hr <- calculate_hazard_ratio(cox_model)
  expect_type(hr, "list")
  expect_true("hazard_ratios" %in% names(hr))
  expect_true("confidence_intervals" %in% names(hr))
})

test_that("log_rank_test works correctly", {
  result <- log_rank_test(test_data, time_col = "time", status_col = "status", group_col = "group")
  expect_type(result, "list")
  expect_true("p_value" %in% names(result))
  expect_true("test_statistic" %in% names(result))
})

