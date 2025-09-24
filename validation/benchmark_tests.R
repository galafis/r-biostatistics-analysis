# ============================================================================
# BENCHMARK TESTS
# Testes de benchmark e comparações de desempenho entre métodos
# Projeto: R Biostatistics Analysis
# Autor: [Nome do Autor]
# Data: 2024-2025
# ============================================================================

suppressPackageStartupMessages({
  require(dplyr)
  require(purrr)
  require(tidyr)
  require(broom)
  require(bench)
  require(ggplot2)
})

#' Benchmark de funções de ajuste e predição
#'
#' @param data dataset base
#' @param methods lista nomeada de listas com fit_fun e predict_fun
#'   ex: list(lm = list(fit = function(d) lm(y~x,d), pred = function(m,newd) predict(m,newd)))
#' @param outcome nome do outcome para métricas
#' @param nrep repetições do benchmark
#' @return lista com tempos (bench::mark) e métricas de acurácia
benchmark_methods <- function(data, methods, outcome, nrep = 10) {
  stopifnot(is.list(methods), outcome %in% names(data))

  # Métricas de acurácia via reamostragem simples holdout repetido
  set.seed(123)
  acc <- map_dfr(seq_len(nrep), function(r){
    idx <- sample(seq_len(nrow(data)), floor(0.7*nrow(data)))
    train <- data[idx, , drop=FALSE]
    test  <- data[-idx, , drop=FALSE]
    map_dfr(names(methods), function(mn){
      fit <- methods[[mn]]$fit(train)
      preds <- methods[[mn]]$pred(fit, test)
      tibble(method = mn,
             rmse = yardstick::rmse_vec(truth = test[[outcome]], estimate = preds),
             mae  = yardstick::mae_vec(truth = test[[outcome]], estimate = preds))
    }) %>% mutate(rep = r)
  })

  # Benchmark de tempo apenas no ajuste (fit)
  bmarks <- bench::mark(
    !!!setNames(lapply(names(methods), function(mn){
      function(){ methods[[mn]]$fit(data) }
    }), names(methods)),
    iterations = nrep,
    check = FALSE
  )

  list(time = bmarks, accuracy = acc,
       accuracy_summary = acc %>% group_by(method) %>% summarise(across(c(rmse, mae), list(mean=mean, sd=sd), .names = "{.col}_{.fn}")))
}

#' Visualização dos resultados de benchmark
plot_benchmark <- function(bench_res){
  p1 <- autoplot(bench_res$time) + ggtitle("Tempo de execução (fit)")
  acc_long <- bench_res$accuracy %>% pivot_longer(c(rmse, mae), names_to = "metric", values_to = "value")
  p2 <- ggplot(acc_long, aes(x = method, y = value, fill = metric)) + geom_boxplot(alpha = 0.6) + theme_minimal() + coord_flip()
  list(time_plot = p1, accuracy_plot = p2)
}

#' Exemplo:
#' df <- tibble(x = rnorm(1000), y = 1 + 2*x + rnorm(1000))
#' methods <- list(
#'   lm = list(fit = function(d) lm(y~x, data=d), pred = function(m,nd) as.numeric(predict(m, nd))),
#'   rlm = list(fit = function(d) MASS::rlm(y~x, data=d), pred = function(m,nd) as.numeric(predict(m, nd)))
#' )
#' res <- benchmark_methods(df, methods, outcome = "y", nrep = 30)
#' plots <- plot_benchmark(res)

cat("benchmark_tests.R carregado com sucesso\n")
