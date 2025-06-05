# Test suite for clinical trials analysis functions
library(testthat)
library(rbiostatistics)

context("Clinical Trials Analysis Functions")

# Sample data for testing
test_data <- data.frame(
  patient_id = 1:20,
  treatment = rep(c("Treatment", "Placebo"), each = 10),
  baseline = c(rnorm(10, mean = 100, sd = 15), rnorm(10, mean = 100, sd = 15)),
  endpoint = c(rnorm(10, mean = 85, sd = 15), rnorm(10, mean = 95, sd = 15)),
  age = sample(40:80, 20, replace = TRUE),
  sex = sample(c("M", "F"), 20, replace = TRUE),
  site = sample(c("Site1", "Site2", "Site3"), 20, replace = TRUE)
)

test_that("calculate_sample_size works correctly", {
  n <- calculate_sample_size(effect_size = 0.5, power = 0.8, alpha = 0.05)
  expect_type(n, "double")
  expect_gt(n, 0)
  
  # Check that increasing power increases sample size
  n_higher_power <- calculate_sample_size(effect_size = 0.5, power = 0.9, alpha = 0.05)
  expect_gt(n_higher_power, n)
})

test_that("randomize_subjects works correctly", {
  subjects <- data.frame(id = 1:100)
  randomized <- randomize_subjects(subjects, groups = c("Treatment", "Placebo"), ratio = c(1, 1))
  
  expect_equal(nrow(randomized), 100)
  expect_true("group" %in% names(randomized))
  
  # Check balanced allocation
  treatment_count <- sum(randomized$group == "Treatment")
  placebo_count <- sum(randomized$group == "Placebo")
  expect_equal(treatment_count, placebo_count)
})

test_that("analyze_treatment_effect works correctly", {
  result <- analyze_treatment_effect(test_data, outcome = "endpoint", 
                                    treatment = "treatment", 
                                    covariates = c("baseline", "age", "sex"))
  
  expect_type(result, "list")
  expect_true("p_value" %in% names(result))
  expect_true("effect_size" %in% names(result))
  expect_true("confidence_interval" %in% names(result))
})

test_that("perform_interim_analysis works correctly", {
  result <- perform_interim_analysis(test_data, outcome = "endpoint", 
                                    treatment = "treatment", 
                                    alpha_spending = 0.025)
  
  expect_type(result, "list")
  expect_true("continue_trial" %in% names(result))
  expect_true("adjusted_p_value" %in% names(result))
})

test_that("calculate_number_needed_to_treat works correctly", {
  # Create binary outcome data
  binary_data <- data.frame(
    treatment = rep(c("Treatment", "Placebo"), each = 50),
    outcome = c(rbinom(50, 1, 0.3), rbinom(50, 1, 0.5))  # 30% vs 50% event rate
  )
  
  nnt <- calculate_number_needed_to_treat(binary_data, 
                                         treatment_col = "treatment", 
                                         outcome_col = "outcome")
  
  expect_type(nnt, "double")
  expect_gt(nnt, 0)
})

test_that("analyze_subgroups works correctly", {
  result <- analyze_subgroups(test_data, 
                             outcome = "endpoint", 
                             treatment = "treatment", 
                             subgroups = c("sex", "site"))
  
  expect_type(result, "list")
  expect_equal(length(result), 2)  # Two subgroup variables
  expect_true(all(c("sex", "site") %in% names(result)))
})

