# diagnostic_accuracy.R
# Diagnostic accuracy measures using epiR and base R
# Author: Project r-biostatistics-analysis
# Functions: compute_accuracy (Se, Sp, PPV, NPV, LR+, LR-, accuracy, F1), epi_confusion_epir
# Example:
#   tab <- matrix(c(80, 20, 10, 90), nrow=2, byrow=TRUE)
#   dimnames(tab) <- list(Test=c('Positive','Negative'), Ref=c('Disease','NoDisease'))
#   res <- compute_accuracy(tab)
#   print(res$measures)

#' Compute diagnostic accuracy measures from 2x2 table
#' @param table2x2 2x2 matrix with rows = Test (Positive, Negative), cols = Reference (Disease, NoDisease)
#' @param ci logical; if TRUE compute Wald 95% CIs
#' @return list with measures data.frame and confusion matrix
#' @export
compute_accuracy <- function(table2x2, ci = TRUE) {
  stopifnot(is.matrix(table2x2), all(dim(table2x2) == c(2,2)))
  TP <- table2x2[1,1]; FP <- table2x2[1,2]; FN <- table2x2[2,1]; TN <- table2x2[2,2]
  Se <- TP/(TP+FN)
  Sp <- TN/(TN+FP)
  PPV <- TP/(TP+FP)
  NPV <- TN/(TN+FN)
  Acc <- (TP+TN)/sum(table2x2)
  LRpos <- Se/(1-Sp)
  LRneg <- (1-Se)/Sp
  F1 <- 2*TP/(2*TP + FP + FN)
  measures <- data.frame(Measure=c('Sensitivity','Specificity','PPV','NPV','Accuracy','LR+','LR-','F1'),
                         Value=c(Se,Sp,PPV,NPV,Acc,LRpos,LRneg,F1))
  if (ci) {
    wald_ci <- function(p,n){ se <- sqrt(p*(1-p)/n); c(p-1.96*se, p+1.96*se) }
    ci_df <- rbind(
      wald_ci(Se, TP+FN),
      wald_ci(Sp, TN+FP),
      wald_ci(PPV, TP+FP),
      wald_ci(NPV, TN+FN),
      wald_ci(Acc, sum(table2x2))
    )
    measures$Lower <- c(ci_df[,1], NA, NA, NA)
    measures$Upper <- c(ci_df[,2], NA, NA, NA)
  }
  list(confusion=table2x2, measures=measures)
}

#' Use epiR to compute a rich set of diagnostic parameters
#' @param table2x2 2x2 matrix (same orientation as compute_accuracy)
#' @return epiR::epi.tests result if epiR is available, else NULL with message
#' @export
epi_confusion_epir <- function(table2x2){
  if (!requireNamespace('epiR', quietly = TRUE)) {
    message("Package 'epiR' not installed; returning NULL")
    return(NULL)
  }
  # epi.tests expects matrix with outcome in rows; adjust accordingly
  # Reorder to: rows = outcome (Disease, NoDisease), cols = test (Positive, Negative)
  m <- rbind(c(table2x2[1,1], table2x2[2,1]), c(table2x2[1,2], table2x2[2,2]))
  dimnames(m) <- list(outcome=c('Disease','NoDisease'), test=c('Positive','Negative'))
  epiR::epi.tests(m)
}
