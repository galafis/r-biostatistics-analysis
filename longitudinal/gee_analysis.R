# gee_analysis.R
# Análise de dados longitudinais com Equações de Estimação Generalizadas (GEE)
# Pacote: geepack; famílias gaussian/binomial/poisson; correlações AR-1/EX/UN

#' Ajusta um modelo GEE para desfechos contínuos, binários ou contagens
#' @param data data.frame em formato longo
#' @param formula fórmula do modelo (ex.: y ~ time * treatment + age)
#' @param id variável de cluster/sujeito (string)
#' @param family família: gaussian(), binomial(), poisson(), etc.
#' @param corstr estrutura de correlação: "ar1", "exchangeable", "independence", "unstructured"
#' @param weights pesos opcionais
#' @return lista com: modelo (geeglm), QIC, coeficientes com IC robustos
#' @examples
#' # library(geepack)
#' # fit <- gee_fit(df, formula = y ~ time * trt, id = "id", family = gaussian(), corstr = "ar1")
#' # summary(fit$model)
#' # fit$coef_ci
gee_fit <- function(data, formula, id, family = gaussian(),
                    corstr = c("ar1","exchangeable","independence","unstructured"),
                    weights = NULL){
  if (!requireNamespace("geepack", quietly=TRUE)) stop("Package 'geepack' is required")
  corstr <- match.arg(corstr)
  stopifnot(is.data.frame(data), id %in% names(data))
  id_vec <- data[[id]]
  model <- geepack::geeglm(formula=formula, id=id_vec, data=data,
                           family=family, corstr=corstr, weights=weights)
  # QIC para comparação de modelos
  q <- tryCatch(geepack::QIC(model), error=function(e) NA)
  # Intervalos de confiança robustos
  sm <- summary(model)
  est <- sm$coefficients
  se_rob <- est[,2]
  z <- qnorm(0.975)
  ci <- cbind(Estimate=est[,1], LCL=est[,1]-z*se_rob, UCL=est[,1]+z*se_rob, `Pr(>|z|)`=est[,4])
  list(model=model, QIC=q, coef_ci=ci)
}

#' Predição marginal a partir de um modelo GEE
#' @param model objeto retornado por geepack::geeglm
#' @param newdata data.frame com combinações de interesse
#' @param type tipo de predição: "link" ou "response"
gee_predict <- function(model, newdata, type=c("response","link")){
  type <- match.arg(type)
  stats::predict(model, newdata=newdata, type=type)
}
