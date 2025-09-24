# ============================================================================
# METHOD VALIDATION
# Validação cruzada de métodos estatísticos e reprodutibilidade
# Projeto: R Biostatistics Analysis
# Autor: [Nome do Autor]
# Data: 2024-2025
# ============================================================================

# Dependências
suppressPackageStartupMessages({
  require(dplyr)
  require(purrr)
  require(tidyr)
  require(broom)
  require(rsample)
  require(yardstick)
  require(ggplot2)
})

#' k-fold cross-validation genérica
#'
#' @param data data.frame com dados
#' @param fit_fun função que recebe (train_data) e retorna um modelo/objeto
#' @param predict_fun função que recebe (fit, new_data) e retorna vetor de predições
#' @param outcome nome da coluna outcome (para métricas supervisionadas) ou NULL
#' @param folds número de dobras k (default: 5)
#' @param repeats repetições de CV (default: 1)
#' @param metrics vetor de métricas yardstick, ex: c("rmse","mae","rsq")
#' @param seed semente para reprodutibilidade
#' @return tibble com métricas por fold e resumo
cv_kfold <- function(data, fit_fun, predict_fun, outcome = NULL,
                     folds = 5, repeats = 1, metrics = c("rmse","mae","rsq"), seed = 123) {
  stopifnot(is.data.frame(data))
  set.seed(seed)

  # Preparar partições
  resamples <- map(1:repeats, ~ vfold_cv(data, v = folds, repeats = 1))
  resamples <- bind_rows(resamples, .id = "repeat")

  # Loop por folds
  evals <- resamples %>% mutate(
    results = pmap(list(splits, repeat, id), function(s, rep_id, fold_id){
      train <- rsample::analysis(s)
      test  <- rsample::assessment(s)

      fit <- fit_fun(train)
      preds <- predict_fun(fit, test)

      if (!is.null(outcome)) {
        dfm <- tibble(truth = test[[outcome]], estimate = preds)
        # Calcular métricas solicitadas
        mets <- list()
        if ("rmse" %in% metrics) mets$rmse <- yardstick::rmse(dfm, truth, estimate)$.estimate
        if ("mae"  %in% metrics) mets$mae  <- yardstick::mae(dfm, truth, estimate)$.estimate
        if ("rsq"  %in% metrics) mets$rsq  <- yardstick::rsq(dfm, truth, estimate)$.estimate
        if ("roc_auc" %in% metrics) {
          # Para classificadores binários com probabilidade
          if (is.numeric(dfm$estimate)) {
            if (is.factor(dfm$truth) || is.character(dfm$truth)) {
              dfm <- dfm %>% mutate(truth = factor(truth))
              mets$roc_auc <- yardstick::roc_auc(dfm, truth, estimate, event_level = "second")$.estimate
            }
          }
        }
        tibble(metric = names(mets), value = unlist(mets)) %>% mutate(repeat = rep_id, fold = fold_id)
      } else {
        tibble(metric = c("n_train","n_test"), value = c(nrow(train), nrow(test)),
               repeat = rep_id, fold = fold_id)
      }
    })
  ) %>% select(results) %>% unnest(results)

  summary <- evals %>% group_by(metric) %>% summarise(mean = mean(value), sd = sd(value), .groups = "drop")

  list(per_fold = evals, summary = summary)
}

#' Bootstrap para avaliação de incerteza das métricas
bootstrap_validation <- function(data, fit_fun, predict_fun, outcome,
                                B = 1000, metrics = c("rmse","mae","rsq"), seed = 123) {
  set.seed(seed)
  boots <- rsample::bootstraps(data, times = B)

  res <- map_dfr(seq_len(B), function(i){
    s <- boots$splits[[i]]
    train <- rsample::analysis(s)
    test  <- rsample::assessment(s)

    fit <- fit_fun(train)
    preds <- predict_fun(fit, test)

    dfm <- tibble(truth = test[[outcome]], estimate = preds)
    out <- list()
    if ("rmse" %in% metrics) out$rmse <- yardstick::rmse(dfm, truth, estimate)$.estimate
    if ("mae"  %in% metrics) out$mae  <- yardstick::mae(dfm, truth, estimate)$.estimate
    if ("rsq"  %in% metrics) out$rsq  <- yardstick::rsq(dfm, truth, estimate)$.estimate
    tibble(.draw = i, !!!out)
  })

  ci <- res %>% pivot_longer(-.draw, names_to = "metric", values_to = "value") %>%
    group_by(metric) %>% summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      p2.5 = quantile(value, 0.025, na.rm = TRUE),
      p97.5 = quantile(value, 0.975, na.rm = TRUE), .groups = "drop"
    )

  list(draws = res, summary = ci)
}

#' Teste de reprodutibilidade entre execuções
reproducibility_test <- function(run_fun, seeds = c(1,2,3,4,5), compare_fun = NULL) {
  if (is.null(compare_fun)) {
    compare_fun <- function(a,b){
      inter <- intersect(names(a), names(b))
      scores <- purrr::map_dbl(inter, function(nm){
        xa <- a[[nm]]; xb <- b[[nm]]
        if (is.numeric(xa) && is.numeric(xb)) {
          1 - (mean(abs(xa - xb), na.rm = TRUE) / (sd(xa, na.rm = TRUE) + sd(xb, na.rm = TRUE) + 1e-8))
        } else { NA_real_ }
      })
      mean(scores, na.rm = TRUE)
    }
  }
  runs <- purrr::map(seeds, ~ run_fun(seed = .x))
  combs <- t(combn(seq_along(runs), 2))
  res <- purrr::map_dfr(seq_len(nrow(combs)), function(i){
    i1 <- combs[i,1]; i2 <- combs[i,2]
    tibble(run1 = i1, run2 = i2, agreement = compare_fun(runs[[i1]], runs[[i2]]))
  })
  res
}

# Exemplos
# set.seed(42)
# n <- 200
# df <- tibble(x = rnorm(n), y = 1 + 2*x + rnorm(n, sd = 0.5))
# fit_fun <- function(d) lm(y ~ x, data = d)
# predict_fun <- function(fit, newd) as.numeric(predict(fit, newd))
# cv <- cv_kfold(df, fit_fun, predict_fun, outcome = "y", folds = 5, repeats = 2)
# boot <- bootstrap_validation(df, fit_fun, predict_fun, outcome = "y", B = 200)
# repro <- reproducibility_test(function(seed){set.seed(seed); list(pred = predict_fun(fit_fun(df), df)[1:10])})

cat("method_validation.R carregado com sucesso\n")
