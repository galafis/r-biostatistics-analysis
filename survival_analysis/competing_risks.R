#' Competing Risks Analysis: Cumulative Incidence and Fine-Gray Models
#'
#' Functions to perform competing risks analysis: estimate cumulative incidence
#' functions (CIF) with Gray's test, and fit subdistribution hazard models
#' (Fine-Gray) for covariate effects. Tidy outputs and plotting helpers.
#'
#' @section Dependencies:
#' - cmprsk (cuminc, crr)
#' - survival (for data prep)
#' - ggplot2 (plotting)
#'
#' @examples
#' # Assume df with time, status (0=censor, 1=event of interest, 2=competing), group
#' # cif <- cr_cif(df, time="time", status="status", group="group", event=1)
#' # fg  <- cr_finegray(df, time="time", status="status", covariates=c("treatment","age"), event=1)
#'
#' @keywords survival competing-risks CIF fine-gray Gray-test subdistribution
#' @author Biostatistics Platform Contributors
#'
#' @importFrom stats as.formula

#' Cumulative incidence estimation and Gray's test
#'
#' @param data data.frame
#' @param time character
#' @param status character (0=censor, 1=event of interest, >1 competing)
#' @param group optional character for groups/strata
#' @param event integer code of event of interest (default 1)
#' @return list with cmprsk::cuminc object and (optional) Gray test p-value
#' @export
cr_cif <- function(data, time, status, group = NULL, event = 1) {
  requireNamespace("cmprsk", quietly = TRUE)
  stopifnot(time %in% names(data), status %in% names(data))
  ftime  <- data[[time]]
  fstatus <- data[[status]]
  if (!is.null(group)) {
    g <- factor(data[[group]])
    ci <- cmprsk::cuminc(ftime = ftime, fstatus = fstatus, group = g, cencode = 0)
  } else {
    ci <- cmprsk::cuminc(ftime = ftime, fstatus = fstatus, cencode = 0)
  }
  # Gray's test is embedded in cuminc when group provided
  pval <- NULL
  if (!is.null(group)) {
    # Extract p-value for event of interest if available
    if (!is.null(ci$Tests)) pval <- ci$Tests[grep(paste0("^", event, "$"), rownames(ci$Tests), perl = TRUE), "pv"]
  }
  list(cif = ci, gray_p = pval)
}

#' Tidy CIF into data.frame for plotting/reporting
#' @param ci cmprsk cuminc object
#' @param event integer event code to extract
#' @return data.frame with time, est, var, group
#' @export
cr_cif_tidy <- function(ci, event = 1) {
  comps <- names(ci)
  keep <- grep(paste0("^", event, " "), comps)
  out <- do.call(rbind, lapply(keep, function(idx) {
    comp <- ci[[idx]]
    grp  <- sub(paste0("^", event, " "), "", names(ci)[idx])
    data.frame(time = comp$time, est = comp$est, var = comp$var, group = grp)
  }))
  rownames(out) <- NULL
  out
}

#' Fine-Gray subdistribution hazards regression
#'
#' @param data data.frame
#' @param time,status characters as above
#' @param covariates character vector
#' @param event integer event code of interest (default 1)
#' @return list with model (cmprsk::crr) and tidy HR-like table (sHR)
#' @export
cr_finegray <- function(data, time, status, covariates, event = 1) {
  requireNamespace("cmprsk", quietly = TRUE)
  stopifnot(all(c(time,status,covariates) %in% names(data)))
  y <- data[[time]]; d <- data[[status]]
  X <- model.matrix(as.formula(paste("~", paste(covariates, collapse = "+"))), data)[, -1, drop = FALSE]
  fit <- cmprsk::crr(ftime = y, fstatus = d, cov1 = X, failcode = event, cencode = 0)
  # Tidy coefficients as subdistribution HR (sHR)
  beta <- fit$coef
  se   <- sqrt(diag(fit$var))
  z    <- beta / se
  p    <- 2*pnorm(abs(z), lower.tail = FALSE)
  sHR  <- exp(beta)
  out <- data.frame(term = names(beta), sHR = sHR,
                    lower = exp(beta - 1.96*se), upper = exp(beta + 1.96*se), p_value = p,
                    row.names = NULL)
  list(model = fit, tidy = out)
}

#' Convenience: full competing risks workflow
#' @export
cr_full_report <- function(data, time, status, group = NULL, covariates = NULL, event = 1) {
  ci <- cr_cif(data, time, status, group, event)
  tidy <- tryCatch(cr_cif_tidy(ci$cif, event = event), error = function(e) NULL)
  fg <- if (!is.null(covariates)) tryCatch(cr_finegray(data, time, status, covariates, event), error = function(e) NULL) else NULL
  list(cif = ci$cif, cif_tidy = tidy, finegray = fg)
}
