# ============================================================================
# PUBLICATION BIAS ASSESSMENT
# Biostatistical Analysis Platform for Medical Research
# ============================================================================

#' @title Publication Bias and Small-Study Effects
#' @description Implements Egger's test, Begg's test, trim-and-fill, and
#'   contour-enhanced funnel plots using metafor/meta.
#' @author Gabriel Demetrios Lafis
#' @date September 2025

library(metafor)
library(meta)
library(ggplot2)

# Core tests ------------------------------------------------------------------

#' Assess publication bias given a metafor rma model
#' @param rma_model An object from metafor::rma
#' @export
publication_bias <- function(rma_model) {
  egger <- tryCatch(metafor::regtest(rma_model, model = "lm"), error = function(e) e)
  begg  <- tryCatch(metafor::ranktest(rma_model), error = function(e) e)
  tf    <- tryCatch(metafor::trimfill(rma_model), error = function(e) e)
  list(egger = egger, begg = begg, trimfill = tf)
}

#' Contour-enhanced funnel plot
#' @export
funnel_contour <- function(rma_model, main = "Contour-enhanced funnel plot") {
  metafor::funnel(rma_model, main = main, shade = c(.9, .75, .5))
}

# Using meta objects -----------------------------------------------------------

#' Publication bias for meta::meta objects (metabin/metacont/metagen)
#' @export
publication_bias_meta <- function(meta_obj) {
  metabias_res <- tryCatch(meta::metabias(meta_obj, method.bias = "linreg"), error = function(e) e)
  funnel(meta_obj)
  list(metabias = metabias_res)
}

# Export utilities -------------------------------------------------------------

#' Export funnel plot and trim-and-fill comparison
#' @export
export_publication_bias <- function(rma_model, filename_prefix) {
  pdf(paste0(filename_prefix, "_funnel.pdf"), width = 7, height = 7)
  funnel_contour(rma_model)
  dev.off()
  
  tf <- tryCatch(trimfill(rma_model), error = function(e) NULL)
  if (!is.null(tf)) {
    pdf(paste0(filename_prefix, "_trimfill_forest.pdf"), width = 8, height = 6)
    forest(tf)
    dev.off()
  }
}

# Example ---------------------------------------------------------------------

#' Example: Publication bias workflow
example_publication_bias <- function() {
  set.seed(123)
  yi <- rnorm(15, 0.2, 0.3)
  vi <- runif(15, 0.02, 0.08)
  mod <- rma(yi = yi, vi = vi, method = "REML")
  res <- publication_bias(mod)
  funnel_contour(mod)
  return(res)
}

# ============================================================================
# END OF PUBLICATION BIAS MODULE
# ============================================================================
