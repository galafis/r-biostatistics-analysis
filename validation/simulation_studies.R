# ============================================================================
# SIMULATION STUDIES
# Estudos de simulação para avaliação de métodos e cenários
# Projeto: R Biostatistics Analysis
# Autor: [Nome do Autor]
# Data: 2024-2025
# ============================================================================

suppressPackageStartupMessages({
  require(dplyr)
  require(purrr)
  require(tidyr)
  require(broom)
  require(ggplot2)
  require(future.apply)
})

#' Função genérica para estudos de simulação
#'
#' @param nsims número de simulações
#' @param dgp função geradora de dados: function(seed, params) -> data.frame
#' @param analyze função de análise: function(data, params) -> lista com estimativas
#' @param params lista de parâmetros (passada a dgp e analyze)
#' @param parallel se TRUE, usa future.apply para paralelizar
#' @param seed semente base
#' @return lista com resultados por simulação e sumários
run_simulation <- function(nsims = 1000, dgp, analyze, params = list(), parallel = TRUE, seed = 123) {
  stopifnot(is.function(dgp), is.function(analyze))
  set.seed(seed)

  idx <- seq_len(nsims)
  do_fun <- function(i){
    di <- dgp(seed + i, params)
    ai <- analyze(di, params)
    tibble(.sim = i, !!!ai)
  }

  if (parallel) {
    oplan <- future::plan()
    on.exit(future::plan(oplan), add = TRUE)
    future::plan(future::multisession, workers = max(1, future::availableCores() - 1))
    res <- future_lapply(idx, do_fun)
    res <- bind_rows(res)
  } else {
    res <- map_dfr(idx, do_fun)
  }

  # Sumários (média, viés, EP empírico, RMSE)
  long <- res %>% pivot_longer(-.sim, names_to = "parameter", values_to = "estimate")
  summary <- long %>% group_by(parameter) %>% summarise(
    mean = mean(estimate, na.rm = TRUE),
    sd = sd(estimate, na.rm = TRUE),
    rmse = sqrt(mean((estimate - mean(estimate, na.rm = TRUE))^2, na.rm = TRUE)),
    .groups = "drop"
  )

  list(results = res, summary = summary)
}

#' Grid de cenários para simulação
#'
#' @param grid tibble com combinações de parâmetros
#' @param dgp função geradora de dados
#' @param analyze função de análise
#' @param nsims número de repetições por cenário
#' @param parallel paralelização
#' @return tibble com resultados agregados por cenário
run_simulation_grid <- function(grid, dgp, analyze, nsims = 1000, parallel = TRUE, seed = 123) {
  stopifnot(is.data.frame(grid))
  grid <- grid %>% mutate(.scenario = row_number())

  res <- purrr::map_dfr(seq_len(nrow(grid)), function(i){
    params <- as.list(grid[i, , drop = FALSE])
    sim <- run_simulation(nsims = nsims, dgp = dgp, analyze = analyze, params = params,
                          parallel = parallel, seed = seed + i)
    sim$summary %>% mutate(.scenario = i)
  })

  res %>% left_join(grid %>% select(.scenario, everything()), by = ".scenario")
}

#' Exemplo: regressão linear sob diferentes tamanhos amostrais e SNR
#'
#' dgp <- function(seed, params){
#'   set.seed(seed)
#'   n <- params$n; beta <- params$beta; sigma <- params$sigma
#'   x <- rnorm(n); y <- beta[1] + beta[2]*x + rnorm(n, sd = sigma)
#'   tibble(x = x, y = y)
#' }
#' analyze <- function(data, params){
#'   fit <- lm(y ~ x, data = data)
#'   c(beta0 = coef(fit)[1], beta1 = coef(fit)[2])
#' }
#' grid <- tidyr::expand_grid(n = c(50, 100, 500), sigma = c(0.5, 1, 2)) %>%
#'   mutate(beta = list(c(1,2)))
#' sim <- run_simulation_grid(grid, dgp, analyze, nsims = 500)

cat("simulation_studies.R carregado com sucesso\n")
