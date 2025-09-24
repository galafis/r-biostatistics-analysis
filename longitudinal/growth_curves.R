# growth_curves.R
# Curvas de crescimento e modelos polinomiais/mistos para trajetórias longitudinais
# Pacotes: lme4, nlme, stats, splines

#' Ajusta curvas de crescimento com polinômios no tempo e efeitos aleatórios
#' @param data data.frame em formato longo
#' @param outcome desfecho contínuo
#' @param time variável de tempo (numérica)
#' @param subject ID do sujeito
#' @param degree grau do polinômio (default 2)
#' @param group fator opcional (ex.: tratamento)
#' @return lista com modelo e previsão
#' @examples
#' # fit <- growth_poly_mixed(df, outcome="y", time="time", subject="id", degree=2, group="trt")
#' # summary(fit$model)
#' # fit$predict_fn(data.frame(id=1, time=seq(0,3,0.5), trt="Ativo"))
growth_poly_mixed <- function(data, outcome, time, subject, degree=2, group=NULL){
  if (!requireNamespace("lme4", quietly=TRUE)) stop("Package 'lme4' is required")
  stopifnot(is.data.frame(data), all(c(outcome,time,subject) %in% names(data)))
  deg <- as.integer(degree)
  if (deg < 1) stop("degree must be >= 1")
  poly_terms <- paste(sprintf("I(%s^%d)", time, 1:deg), collapse = " + ")
  rhs <- poly_terms
  if (!is.null(group)) rhs <- paste(rhs, "+", group, "+", paste0("(", group, ")*", poly_terms))
  rand <- sprintf("(1 + %s | %s)", time, subject)
  fml <- stats::as.formula(sprintf("%s ~ %s + %s", outcome, rhs, rand))
  model <- lme4::lmer(fml, data=data, REML=TRUE)
  predict_fn <- function(newdata, re.form=NA){
    stats::predict(model, newdata=newdata, re.form=re.form, allow.new.levels=TRUE)
  }
  list(model=model, formula=fml, predict_fn=predict_fn)
}

#' Ajusta curvas de crescimento não lineares com splines cúbicos naturais
#' @param df_knots número de graus de liberdade do spline (default 3)
#' @examples
#' # fit <- growth_spline_mixed(df, outcome="y", time="time", subject="id", df_knots=4, group="trt")
#' # summary(fit$model)
#' growth_spline_mixed <- function(data, outcome, time, subject, df_knots=3, group=NULL){
  if (!requireNamespace("lme4", quietly=TRUE)) stop("Package 'lme4' is required")
  if (!requireNamespace("splines", quietly=TRUE)) stop("Package 'splines' is required")
  stopifnot(is.data.frame(data), all(c(outcome,time,subject) %in% names(data)))
  spline <- sprintf("splines::ns(%s, df=%d)", time, as.integer(df_knots))
  rhs <- spline
  if (!is.null(group)) rhs <- paste(rhs, "+", group, "+", paste0("(", group, ")*", spline))
  rand <- sprintf("(1 + %s | %s)", time, subject)
  fml <- stats::as.formula(sprintf("%s ~ %s + %s", outcome, rhs, rand))
  model <- lme4::lmer(fml, data=data, REML=TRUE)
  list(model=model, formula=fml)
}
