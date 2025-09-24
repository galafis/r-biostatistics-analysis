# roc_analysis.R
# Diagnostic ROC analysis utilities using pROC
# Author: Project r-biostatistics-analysis
# Description: Functions to compute ROC curves, AUC (with CI), optimal thresholds, and plotting.
# Dependencies: pROC
# Example:
#   set.seed(1)
#   score <- rnorm(100)
#   status <- rbinom(100, 1, plogis(score))
#   res <- roc_analyze(score, status)
#   print(res$auc)
#   res$plot()

#' Compute ROC analysis with AUC, CI, and optimal cutoff
#' @param predictor numeric vector of test scores/probabilities (higher = more likely positive)
#' @param outcome binary vector/factor (1/TRUE/"positive" = disease present)
#' @param direction passed to pROC::roc (default auto); use ">" if higher predictor means more positive
#' @param ci logical; compute 95% CI for AUC (default TRUE)
#' @param ci_method method for CI ("delong" default), see pROC::ci.auc
#' @return list with roc (pROC object), auc, auc_ci (if requested), coords_best (Youden), sensitivity, specificity, threshold, and a plot() closure
#' @export
roc_analyze <- function(predictor, outcome, direction = NULL, ci = TRUE, ci_method = c("delong","bootstrap","venkatraman")) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required. Please install it: install.packages('pROC')")
  }
  ci_method <- match.arg(ci_method)
  # Coerce outcome to binary 0/1
  if (is.factor(outcome)) {
    if (nlevels(outcome) != 2) stop("Outcome must have 2 levels")
    outcome_bin <- as.integer(outcome == levels(outcome)[2])
  } else if (is.logical(outcome)) {
    outcome_bin <- as.integer(outcome)
  } else {
    outcome_bin <- as.integer(outcome)
  }
  roc_obj <- pROC::roc(response = outcome_bin, predictor = predictor, direction = direction)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  auc_ci <- if (isTRUE(ci)) as.numeric(pROC::ci.auc(roc_obj, method = ci_method)) else NULL
  best <- pROC::coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity","ppv","npv"))
  plot_fun <- function(...) {
    pROC::plot.roc(roc_obj, print.auc = TRUE, legacy.axes = TRUE, ...)
    invisible(roc_obj)
  }
  list(
    roc = roc_obj,
    auc = auc_val,
    auc_ci = auc_ci,
    threshold = unname(best["threshold"]),
    sensitivity = unname(best["sensitivity"]),
    specificity = unname(best["specificity"]),
    ppv = unname(best["ppv"]),
    npv = unname(best["npv"]),
    plot = plot_fun
  )
}

#' Compare two ROC curves (paired or unpaired) using DeLong test
#' @param predictor1 numeric scores for model/test 1
#' @param predictor2 numeric scores for model/test 2
#' @param outcome common binary outcome (if paired). If unpaired, provide outcome1/outcome2 instead.
#' @param paired logical; if TRUE uses paired DeLong, else unpaired
#' @param outcome1,outcome2 optional outcomes for unpaired comparison
#' @return list with test statistic and p-value from pROC::roc.test
#' @export
roc_compare <- function(predictor1, predictor2, outcome = NULL, paired = TRUE, outcome1 = NULL, outcome2 = NULL) {
  if (!requireNamespace("pROC", quietly = TRUE)) stop("Package 'pROC' is required")
  if (paired) {
    if (is.null(outcome)) stop("Provide common outcome for paired comparison")
    roc1 <- pROC::roc(outcome, predictor1)
    roc2 <- pROC::roc(outcome, predictor2)
  } else {
    if (is.null(outcome1) || is.null(outcome2)) stop("Provide outcome1 and outcome2 for unpaired comparison")
    roc1 <- pROC::roc(outcome1, predictor1)
    roc2 <- pROC::roc(outcome2, predictor2)
  }
  tst <- pROC::roc.test(roc1, roc2, paired = paired, method = "delong")
  list(statistic = unname(tst$statistic), p_value = tst$p.value, method = tst$method)
}
