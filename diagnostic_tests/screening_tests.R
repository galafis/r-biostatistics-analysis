# screening_tests.R
# Screening test evaluation utilities
# Author: Project r-biostatistics-analysis
# Functions: number_needed_to_screen, net_benefit_decision_curve (simple), prevalence_standardization
# Example:
#   nns <- number_needed_to_screen(sensitivity=0.9, specificity=0.95, prevalence=0.02, harm=0.001)
#   nns

#' Number needed to screen (NNS) accounting for harms
#' @param sensitivity test sensitivity (0-1)
#' @param specificity test specificity (0-1)
#' @param prevalence disease prevalence (0-1)
#' @param harm per-person harm or disutility of screening (0-1), e.g., 0.001
#' @return estimated NNS to achieve one additional true case detected net of harms
#' @export
number_needed_to_screen <- function(sensitivity, specificity, prevalence, harm = 0){
  stopifnot(sensitivity>=0, sensitivity<=1, specificity>=0, specificity<=1,
            prevalence>=0, prevalence<=1, harm>=0)
  # Expected net benefit per screened person: TP gain - FP harm
  tp_gain <- sensitivity * prevalence
  fp_rate <- (1 - specificity) * (1 - prevalence)
  net_benefit <- tp_gain - harm - fp_rate * harm
  if (net_benefit <= 0) return(Inf)
  1 / net_benefit
}

#' Simple net benefit across threshold probabilities (decision curve)
#' @param risk_scores vector of predicted risks (0-1)
#' @param outcome binary outcome (0/1)
#' @param thresholds numeric vector of decision thresholds (0-1)
#' @return data.frame with threshold and net benefit
#' @export
net_benefit_decision_curve <- function(risk_scores, outcome, thresholds = seq(0.05,0.5,by=0.05)){
  y <- as.integer(outcome)
  s <- as.numeric(risk_scores)
  out <- lapply(thresholds, function(t){
    treat <- s >= t
    tp <- mean(treat & y==1)
    fp <- mean(treat & y==0)
    nb <- tp - fp * (t/(1-t))
    c(threshold=t, net_benefit=nb)
  })
  as.data.frame(do.call(rbind, out))
}

#' Adjust PPV and NPV for a new prevalence via Bayes theorem
#' @param sensitivity specificity values (0-1)
#' @param specificity specificity (0-1)
#' @param prevalence_new target prevalence (0-1)
#' @return list with PPV_new and NPV_new
#' @export
prevalence_standardization <- function(sensitivity, specificity, prevalence_new){
  stopifnot(sensitivity>=0, sensitivity<=1, specificity>=0, specificity<=1,
            prevalence_new>=0, prevalence_new<=1)
  ppv_num <- sensitivity * prevalence_new
  ppv_den <- ppv_num + (1-specificity) * (1-prevalence_new)
  ppv_new <- if (ppv_den>0) ppv_num/ppv_den else NA_real_
  npv_num <- specificity * (1-prevalence_new)
  npv_den <- npv_num + (1-sensitivity) * prevalence_new
  npv_new <- if (npv_den>0) npv_num/npv_den else NA_real_
  list(PPV_new = ppv_new, NPV_new = npv_new)
}
