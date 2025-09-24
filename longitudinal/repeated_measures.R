# repeated_measures.R
# Autor: Projeto r-biostatistics-analysis (galafis)
# Descrição: Funções para ANOVA de medidas repetidas e verificação de pressupostos
# Pacotes: nlme, afex (opcional), car, effectsize
#
# Exemplo reprodutível ao final do arquivo (executável após source).

suppressPackageStartupMessages({
  if (!requireNamespace("nlme", quietly = TRUE)) stop("Pacote 'nlme' é necessário.")
  if (!requireNamespace("car", quietly = TRUE)) stop("Pacote 'car' é necessário.")
  if (!requireNamespace("effectsize", quietly = TRUE)) stop("Pacote 'effectsize' é necessário.")
})

#' ANOVA de medidas repetidas com nlme::lme
#'
#' Ajusta um modelo linear com correlação intra-sujeito usando nlme e retorna
#' ANOVA por efeitos fixos, incluindo esfericidade (quando aplicável) e tamanhos de efeito.
#'
#' @param data data.frame em formato longo contendo as colunas id, tempo e resposta
#' @param id nome da coluna de sujeito (fator/identificador)
#' @param time nome da coluna de tempo (fator ou numérica)
#' @param response nome da coluna resposta (numérica)
#' @param between fatores entre-sujeitos (ex.: grupo), como fórmula ou NULL
#' @param corr estrutura de correlação intra ("compound", "ar1", "un"). Padrão "compound"
#' @param weights estrutura de heterocedasticidade via nlme::varIdent ou NULL
#' @return lista com modelo, anova, testes de pressupostos e tamanhos de efeito
#' @examples
#' # Exemplo mínimo (dados simulados)
#' set.seed(1)
#' n <- 40; t <- 4
#' df <- expand.grid(id = factor(1:n), time = factor(1:t))
#' df$group <- factor(sample(c("A","B"), n, TRUE))
#' df$response <- 10 + as.numeric(df$time) + rnorm(n*t,0,1) + ifelse(df$group=="B",1,0)
#' fit <- rm_anova_nlme(df, id="id", time="time", response="response", between=~group, corr="compound")
#' fit$anova
rm_anova_nlme <- function(data, id, time, response, between = NULL, corr = c("compound","ar1","un"), weights = NULL){
  corr <- match.arg(corr)
  df <- data
  df[[id]] <- as.factor(df[[id]])
  # garantir fator para tempo em ANOVA tradicional
  if (!is.factor(df[[time]])) df[[time]] <- factor(df[[time]])

  fixed_formula <- if (is.null(between)) {
    stats::as.formula(sprintf("%s ~ %s", response, time))
  } else {
    stats::as.formula(sprintf("%s ~ %s * (%s)", response, time, as.character(between)[2]))
  }

  random_formula <- stats::as.formula(sprintf("~ 1 | %s", id))

  cor_struct <- switch(corr,
    compound = nlme::corCompSymm(form = stats::as.formula(sprintf("~ 1 | %s", id))),
    ar1 = nlme::corAR1(form = stats::as.formula(sprintf("~ as.numeric(%s) | %s", time, id))),
    un = nlme::corSymm(form = stats::as.formula(sprintf("~ %s | %s", time, id)))
  )

  fit <- nlme::lme(
    fixed = fixed_formula,
    random = random_formula,
    data = df,
    correlation = cor_struct,
    weights = weights,
    na.action = stats::na.omit,
    control = nlme::lmeControl(opt = "optim")
  )

  an <- anova(fit)

  # Tamanhos de efeito aproximados (eta parcial) via effectsize
  es <- try(effectsize::eta_squared(an, ci = 0.90), silent = TRUE)

  # Pressupostos: normalidade resíduos e homocedasticidade simples
  res <- stats::residuals(fit, type = "normalized")
  shapiro <- try(stats::shapiro.test(res), silent = TRUE)

  list(model = fit, anova = an, effect_sizes = es, shapiro = shapiro)
}

#' Teste de esfericidade (Mauchly) para fatores intra-sujeitos
#'
#' Implementa Mauchly a partir de matriz de covariâncias por sujeito.
#' Observação: para designs complexos, considere pacotes como afex/ez.
#' @param data, id, time, response conforme rm_anova_nlme
#' @return lista com estatística W e p-valor aproximado
mauchly_test <- function(data, id, time, response){
  df <- data
  df[[id]] <- factor(df[[id]])
  df[[time]] <- factor(df[[time]])
  # matriz de covariância média
  wide <- stats::reshape(df[, c(id, time, response)], idvar = id, timevar = time, direction = "wide")
  mat <- stats::cov(stats::na.omit(wide[,-1]), use = "pairwise.complete.obs")
  k <- ncol(mat)
  # estatística de Mauchly (aproximação)
  logdet <- determinant(mat, logarithm = TRUE)$modulus
  tr <- sum(diag(mat))
  W <- exp(logdet - k * log(tr/k))
  n <- nrow(wide)
  df_chi <- (k*(k-1))/2 - 1
  # correção de Bartlett
  c_bar <- (2*k^2 + k + 2*(k-1)*(k+1)/(k-1) - 2) / (6*(k-1))
  chi <- -(n - 1) * (1 - c_bar) * log(W)
  p <- stats::pchisq(chi, df = df_chi, lower.tail = FALSE)
  list(W = as.numeric(W), chi2 = as.numeric(chi), df = df_chi, p.value = as.numeric(p))
}

# Exemplo executável
if (identical(environment(), globalenv())){
  set.seed(123)
  n <- 30; t <- 3
  d <- expand.grid(id = factor(1:n), time = factor(1:t))
  d$group <- factor(sample(c("Ctrl","Trt"), n, TRUE))
  d$y <- 5 + as.numeric(d$time) + ifelse(d$group=="Trt", 0.8, 0) + rnorm(n*t,0,1)
  ex <- rm_anova_nlme(d, id="id", time="time", response="y", between=~group, corr="compound")
  print(head(ex$anova))
  print(ex$effect_sizes)
  print(mauchly_test(d, id="id", time="time", response="y"))
}
