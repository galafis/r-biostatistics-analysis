# validation_functions.R
# Funções para validação de métodos e resultados (bootstrap, simulações, checks)
# Autor: R-Biostatistics-Analysis Project
# Data: 2025-09-24

require(dplyr)
require(purrr)
require(tidyr)
require(broom)

#' Bootstrap para estatística arbitrária
#'
#' Executa bootstrap não paramétrico de uma estatística definida pelo usuário
#'
#' @param data DataFrame ou vetor
#' @param stat_fun Função que recebe (data, ...) e retorna escalar
#' @param R Número de reamostragens
#' @param conf Nível de confiança para IC (default 0.95)
#' @param seed Semente para reprodutibilidade (opcional)
#' @param ... Outros argumentos passados a stat_fun
#' @return Lista com estimativa, distribuição bootstrap, IC percentil e BCa (quando possível)
#' @examples
#' set.seed(123)
#' x <- rnorm(50)
#' boot_mean <- bootstrap_stat(x, mean, R = 2000)
#' boot_mean$ci
bootstrap_stat <- function(data, stat_fun, R = 2000, conf = 0.95, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  n <- if (is.vector(data)) length(data) else nrow(data)
  idx <- replicate(R, sample.int(n, n, replace = TRUE))
  stats <- apply(idx, 2, function(i) {
    d <- if (is.vector(data)) data[i] else data[i, , drop = FALSE]
    stat_fun(d, ...)
  })
  est <- if (is.vector(data)) stat_fun(data, ...) else stat_fun(data, ...)
  alpha <- (1 - conf) / 2
  ci_percentile <- quantile(stats, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
  list(
    estimate = est,
    boot_dist = stats,
    ci = list(percentile = ci_percentile)
  )
}

#' Simulação de poder para teste t de duas amostras
#'
#' Calcula poder estatístico por simulação Monte Carlo
#'
#' @param n1,n2 Tamanhos amostrais dos grupos
#' @param delta Diferença de médias verdadeira
#' @param sd Desvio-padrão comum
#' @param alpha Nível de significância
#' @param reps Número de simulações
#' @param seed Semente (opcional)
#' @return Poder estimado e intervalo de confiança binomial
#' @examples
#' sim_t_power(n1 = 30, n2 = 30, delta = 0.5, sd = 1, reps = 5000)
sim_t_power <- function(n1, n2, delta, sd = 1, alpha = 0.05, reps = 5000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sig <- replicate(reps, {
    x <- rnorm(n1, mean = 0, sd = sd)
    y <- rnorm(n2, mean = delta, sd = sd)
    t.test(x, y, var.equal = TRUE)$p.value < alpha
  })
  p_hat <- mean(sig)
  se <- sqrt(p_hat * (1 - p_hat) / reps)
  ci <- p_hat + c(-1, 1) * qnorm(0.975) * se
  list(power = p_hat, ci95 = pmin(pmax(ci, 0), 1), reps = reps)
}

#' Validação cruzada k-fold para modelos
#'
#' Suporta funções de ajuste e predição arbitrárias
#'
#' @param data DataFrame
#' @param k Número de folds
#' @param fit Função de ajuste: function(train) retorna objeto modelo
#' @param predict_fun Função de predição: function(model, test) retorna vetor de predições
#' @param metric Função de métrica: function(y_true, y_pred)
#' @param y_col Nome da coluna resposta (se necessário pela métrica)
#' @param seed Semente para reprodutibilidade
#' @return Tibble com métricas por fold e média
#' @examples
#' data(mtcars)
#' fit <- function(train) lm(mpg ~ hp + wt, data = train)
#' pred <- function(model, test) predict(model, newdata = test)
#' rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))
#' kfold_cv(mtcars, k = 5, fit = fit, predict_fun = pred, metric = rmse, y_col = "mpg")
kfold_cv <- function(data, k = 5, fit, predict_fun, metric, y_col) {
  n <- nrow(data)
  folds <- sample(rep(1:k, length.out = n))
  results <- vector("list", k)
  for (i in 1:k) {
    train <- data[folds != i, , drop = FALSE]
    test <- data[folds == i, , drop = FALSE]
    model <- fit(train)
    y_true <- test[[y_col]]
    y_pred <- predict_fun(model, test)
    results[[i]] <- tibble(
      fold = i,
      metric = metric(y_true, y_pred)
    )
  }
  res <- bind_rows(results)
  summary <- summarise(res, mean_metric = mean(metric), sd_metric = sd(metric))
  list(by_fold = res, summary = summary)
}

#' Checks de qualidade de dados
#'
#' Verifica valores ausentes, outliers simples e tipos de coluna
#'
#' @param data DataFrame
#' @param id_cols Colunas identificadoras a ignorar (opcional)
#' @return Lista com relatórios de qualidade
#' @examples
#' data_quality_checks(mtcars)
data_quality_checks <- function(data, id_cols = NULL) {
  stopifnot(is.data.frame(data))
  cols <- setdiff(names(data), id_cols %||% character())
  # Missings
  missing_df <- tibble(
    variable = cols,
    n_missing = sapply(cols, function(c) sum(is.na(data[[c]]))),
    prop_missing = sapply(cols, function(c) mean(is.na(data[[c]])))
  )
  # Outliers por regra IQR para numéricos
  is_num <- sapply(data[cols], is.numeric)
  num_cols <- cols[is_num]
  outlier_df <- map_dfr(num_cols, function(c) {
    x <- data[[c]]
    q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
    iqr <- q[2] - q[1]
    lower <- q[1] - 1.5 * iqr
    upper <- q[2] + 1.5 * iqr
    tibble(variable = c,
           n_outliers = sum(x < lower | x > upper, na.rm = TRUE),
           lower_bound = lower,
           upper_bound = upper)
  })
  # Tipos
  types_df <- tibble(
    variable = cols,
    class = sapply(data[cols], function(x) paste(class(x), collapse = ","))
  )
  list(missing = missing_df, outliers = outlier_df, types = types_df)
}

#' Simulação genérica para estimadores
#'
#' Permite avaliar viés e RMSE de um estimador sob um gerador de dados
#'
#' @param rgen Função geradora de dados: function() retorna amostra
#' @param estimator Função estimadora: function(sample) retorna escalar
#' @param true_value Valor verdadeiro do parâmetro
#' @param reps Número de repetições
#' @param seed Semente (opcional)
#' @return Tibble com média, viés, variância e RMSE
#' @examples
#' rgen <- function() rnorm(30, mean = 1, sd = 2)
#' estimator <- mean
#' evaluate_estimator(rgen, estimator, true_value = 1, reps = 2000)
evaluate_estimator <- function(rgen, estimator, true_value, reps = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  vals <- replicate(reps, estimator(rgen()))
  tibble(
    mean_est = mean(vals),
    bias = mean(vals) - true_value,
    var = var(vals),
    rmse = sqrt(mean((vals - true_value)^2))
  )
}

#' Verificação de suposições do modelo linear
#'
#' Retorna métricas e diagnósticos simples (normalidade, homocedasticidade, independência)
#'
#' @param model Objeto de classe lm
#' @return Lista com testes de Shapiro-Wilk, Breusch-Pagan e Durbin-Watson (quando disponível)
#' @examples
#' fit <- lm(mpg ~ hp + wt, data = mtcars)
#' check_lm_assumptions(fit)
check_lm_assumptions <- function(model) {
  stopifnot(inherits(model, "lm"))
  res <- residuals(model)
  shapiro <- tryCatch(shapiro.test(res), error = function(e) NULL)
  bp <- tryCatch({
    if (requireNamespace("lmtest", quietly = TRUE)) lmtest::bptest(model) else NULL
  }, error = function(e) NULL)
  dw <- tryCatch({
    if (requireNamespace("lmtest", quietly = TRUE)) lmtest::dwtest(model) else NULL
  }, error = function(e) NULL)
  list(shapiro = shapiro, breusch_pagan = bp, durbin_watson = dw)
}

# Nota: Funções seguem padrão modular, retorno tidy/gt-ready e exemplos reprodutíveis.
