# ============================================================================
# NETWORK META-ANALYSIS (Frequentist)
# Biostatistical Analysis Platform for Medical Research
# ============================================================================

#' @title Network Meta-Analysis
#' @description Frequentist network meta-analysis using 'netmeta'. Supports
#'   pairwise data, random/fixed effects, ranking (P-scores), and plots.
#' @author Gabriel Demetrios Lafis
#' @date September 2025

# Required Libraries ----------------------------------------------------------
library(netmeta)
library(meta)
library(ggplot2)
library(dplyr)

# Pairwise Preparation --------------------------------------------------------

#' Prepare pairwise comparisons for network meta-analysis
#' @param data Data frame with trial-level arm data
#' @param studlab Trial identifier
#' @param treat Treatment variable
#' @param event Events (for binary outcomes)
#' @param n Sample size (for binary outcomes)
#' @param mean,sd Numeric outcomes (for continuous outcomes)
#' @param type Outcome type: "binary" or "continuous"
#' @param sm Summary measure: e.g., OR, RR, RD, MD, SMD
#' @export
prepare_pairwise <- function(data,
                             studlab, treat,
                             event = NULL, n = NULL,
                             mean = NULL, sd = NULL,
                             type = "binary",
                             sm = "OR") {
  if (type == "binary") {
    pw <- pairwise(treat = data[[treat]],
                   event = data[[event]],
                   n = data[[n]],
                   studlab = data[[studlab]],
                   data = data,
                   sm = sm)
  } else {
    pw <- pairwise(treat = data[[treat]],
                   mean = data[[mean]],
                   sd = data[[sd]],
                   n = data[[n]],
                   studlab = data[[studlab]],
                   data = data,
                   sm = sm)
  }
  return(pw)
}

# Main Network Meta-Analysis --------------------------------------------------

#' Run Network Meta-Analysis
#' @param pw Pairwise object from prepare_pairwise
#' @param random Logical; random effects (default TRUE)
#' @param reference Reference treatment (baseline)
#' @param method.tau Between-study variance method (e.g., "REML", "DL")
#' @param small.values Are smaller values better (for ranking)
#' @export
network_meta <- function(pw,
                         random = TRUE,
                         reference = NULL,
                         method.tau = "REML",
                         small.values = "good") {
  
  nm <- netmeta(TE = pw$TE, seTE = pw$seTE, treat1 = pw$treat1, treat2 = pw$treat2,
                studlab = pw$studlab,
                sm = attr(pw, "sm"),
                comb.fixed = !random,
                comb.random = random,
                reference.group = reference,
                method.tau = method.tau,
                small.values = small.values)
  
  # Ranking (P-scores)
  ranking <- netrank(nm, small.values = small.values)
  
  list(model = nm, ranking = ranking)
}

# Consistency and Inconsistency ----------------------------------------------

#' Assess global and local inconsistency
#' @export
network_inconsistency <- function(nm) {
  decomp <- decomp.design(nm$model)
  netheat <- netheat(nm$model, random = nm$model$random)
  netsplit <- netsplit(nm$model)
  list(decomp = decomp, netheat = netheat, netsplit = netsplit)
}

# Plots -----------------------------------------------------------------------

#' Network geometry plot
#' @export
plot_network_geometry <- function(nm) {
  netgraph(nm$model, number.of.studies = TRUE, plastic = FALSE)
}

#' Forest and league table
#' @export
plot_network_forest <- function(nm) {
  forest(nm$model)
}

#' League table (relative effects between all treatments)
#' @export
league_table <- function(nm, digits = 2) {
  return(netleague(nm$model, bracket = "(" , digits = digits))
}

# Export ----------------------------------------------------------------------

#' Export network meta-analysis results
#' @export
export_network_meta <- function(nm, filename_prefix) {
  # Save forest plot
  pdf(paste0(filename_prefix, "_forest.pdf"), width = 10, height = 8)
  forest(nm$model)
  dev.off()
  
  # Save network plot
  pdf(paste0(filename_prefix, "_network.pdf"), width = 8, height = 8)
  netgraph(nm$model)
  dev.off()
  
  # Save league table
  lt <- capture.output(print(netleague(nm$model)))
  writeLines(lt, paste0(filename_prefix, "_league.txt"))
}

# Example ---------------------------------------------------------------------

#' Example: Network Meta-Analysis
example_network_meta <- function() {
  # Example synthetic binary data
  set.seed(123)
  df <- data.frame(
    studlab = rep(paste0("Trial", 1:6), each = 2),
    treat = c("A", "B",  "A", "C",  "B", "C",  "A", "D",  "B", "D",  "C", "D"),
    event = rbinom(12, size = 100, prob = runif(12, 0.2, 0.5)),
    n = 100
  )
  
  pw <- prepare_pairwise(df, studlab = "studlab", treat = "treat",
                         event = "event", n = "n", type = "binary", sm = "OR")
  nm <- network_meta(pw, random = TRUE, reference = "A")
  plot_network_geometry(nm)
  plot_network_forest(nm)
  print(league_table(nm))
  return(nm)
}

# ============================================================================
# END OF NETWORK META-ANALYSIS MODULE
# ============================================================================
