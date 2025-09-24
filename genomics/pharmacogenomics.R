# pharmacogenomics.R
# Farmacogenômica: análise de associação droga-genótipo, metabolismo, predição de dose
# Pacotes: snpStats, survival (tempo até evento adverso), VGAM/gnm para modelos

#' Associação entre genótipo e resposta ao fármaco (contínua/binária)
#' @param G matriz genotípica (indiv x SNP, 0/1/2)
#' @param response vetor de resposta (numérico contínuo ou binário 0/1)
#' @param covars data.frame de covariáveis clínicas (idade, sexo, etc.)
#' @param family "gaussian" ou "binomial"
#' @return data.frame com beta, se, p por SNP
pharmaco_assoc <- function(G, response, covars=NULL, family=c("gaussian","binomial")){
  family <- match.arg(family)
  n <- nrow(G)
  df <- data.frame(y=response)
  if (!is.null(covars)) df <- cbind(df, covars)
  test <- function(g){
    if (family=="binomial"){
      fit <- try(stats::glm(y ~ g + ., data=df, family=stats::binomial()), silent=TRUE)
    } else {
      fit <- try(stats::lm(y ~ g + ., data=df), silent=TRUE)
    }
    if (inherits(fit, "try-error")) return(c(NA,NA,NA))
    co <- summary(fit)$coefficients
    if (!("g" %in% rownames(co))) return(c(NA,NA,NA))
    c(beta=co["g",1], se=co["g",2], p=co["g",4])
  }
  res <- t(apply(as.matrix(G), 2, test))
  res <- as.data.frame(res)
  res$snp <- colnames(G)
  res$padj_bh <- stats::p.adjust(res$p, method="BH")
  res
}

#' Predição de dose terapêutica baseada em genótipo (ex.: warfarina)
#' @param data data.frame com dose e variáveis clínicas/genéticas
#' @param dose_var nome da coluna de dose
#' @param predictors variáveis preditoras (incluindo SNPs codificados 0/1/2)
#' @return modelo ajustado e função de predição
dose_prediction <- function(data, dose_var, predictors){
  stopifnot(dose_var %in% names(data), all(predictors %in% names(data)))
  form <- stats::as.formula(paste(dose_var, "~", paste(predictors, collapse=" + ")))
  fit <- stats::lm(form, data=data)
  predict_fn <- function(newdata){ stats::predict(fit, newdata=newdata) }
  list(model=fit, predict_fn=predict_fn)
}
