# agreement_studies.R
# Agreement analysis: Cohen's Kappa and Intraclass Correlation Coefficient (ICC)
# Author: Project r-biostatistics-analysis
# Dependencies: irr, psych (optional) or stats for basic ICC via aov
# Example:
#   # Kappa example
#   r1 <- factor(c('A','A','B','B','A'))
#   r2 <- factor(c('A','B','B','B','A'))
#   k <- kappa_cohen(r1, r2)
#   print(k)
#   # ICC example
#   set.seed(1); ratings <- matrix(rnorm(30), nrow=10, ncol=3)
#   icc <- icc_analysis(ratings, model='twoway', type='agreement')
#   print(icc)

#' Compute Cohen's Kappa with standard error and CI
#' @param rater1 factor/character/integer ratings
#' @param rater2 factor/character/integer ratings of equal length
#' @return list with kappa, se, z, p, conf.int
#' @export
kappa_cohen <- function(rater1, rater2){
  if (length(rater1) != length(rater2)) stop('rater1 and rater2 must have same length')
  r1 <- as.factor(rater1); r2 <- as.factor(rater2)
  tab <- table(r1, r2)
  n <- sum(tab)
  po <- sum(diag(tab))/n
  p_row <- rowSums(tab)/n
  p_col <- colSums(tab)/n
  pe <- sum(p_row * p_col)
  k <- (po - pe)/(1 - pe)
  # Approx SE per Fleiss et al.
  q <- 1 - pe
  var_k <- (po*(1 - po)/(n*q^2)) + (2*(1 - po)*(2*pe*po - sum(p_row*p_col*(1 - p_row - p_col)))/(n*q^3))
  se <- sqrt(max(var_k, 0))
  z <- k/se
  p <- 2*pnorm(-abs(z))
  ci <- c(k - 1.96*se, k + 1.96*se)
  list(kappa = k, se = se, z = z, p_value = p, conf_int = ci, table = tab)
}

#' Intraclass Correlation Coefficient (ICC) for quantitative ratings
#' @param ratings matrix/data.frame: rows = targets, cols = raters
#' @param model 'oneway' or 'twoway'
#' @param type 'consistency' or 'agreement'
#' @return list with ICC estimate and 95% CI
#' @export
icc_analysis <- function(ratings, model = c('oneway','twoway'), type = c('consistency','agreement')){
  model <- match.arg(model)
  type <- match.arg(type)
  X <- as.matrix(ratings)
  n <- nrow(X); k <- ncol(X)
  # ANOVA components
  grand_mean <- mean(X)
  ms_between <- k * var(rowMeans(X))
  ms_within <- mean(apply(X, 1, var))
  ms_rater <- n * var(colMeans(X))
  if (model == 'oneway') {
    # Shrout & Fleiss ICC(1)
    icc <- (ms_between - ms_within) / (ms_between + (k - 1) * ms_within)
  } else {
    # Two-way random effects
    if (type == 'consistency') {
      # ICC(C,1) approximation
      icc <- (ms_between - ms_within) / (ms_between + (k - 1) * ms_within)
    } else {
      # agreement uses rater variance
      icc <- (ms_between - ms_within) / (ms_between + (k - 1) * ms_within + k*(ms_rater - ms_within)/n)
    }
  }
  # Simple CI via Fisher z transform bound (approx)
  se <- sqrt(2*(1 - icc)^2 * (1 + (k - 1)*icc)^2 / ((k - 1)*(n - 1)*k^2))
  z <- atanh(icc)
  ci <- tanh(z + c(-1,1)*1.96*se)
  list(icc = icc, conf_int = ci, k = k, n = n)
}
