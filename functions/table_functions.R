# table_functions.R
# Funções para criação de tabelas profissionais em análises biostatísticas
# Autor: R-Biostatistics-Analysis Project
# Data: 2025-09-24

# Carregar bibliotecas necessárias
require(gt)
require(dplyr)
require(tibble)

#' Criar tabela descritiva básica
#'
#' Esta função cria uma tabela descritiva com estatísticas básicas
#' para variáveis numéricas e categóricas
#'
#' @param data DataFrame com os dados
#' @param variables Vetor com nomes das variáveis a incluir
#' @param group_var Variável para estratificação (opcional)
#' @param title Título da tabela
#' @return Objeto gt com tabela formatada
#' @examples
#' # Exemplo básico
#' data(mtcars)
#' create_descriptive_table(mtcars, c("mpg", "hp", "cyl"), title = "Características dos Veículos")
#'
#' # Com estratificação
#' create_descriptive_table(mtcars, c("mpg", "hp"), group_var = "cyl")
create_descriptive_table <- function(data, variables, group_var = NULL, title = "Tabela Descritiva") {
  
  if (!is.null(group_var)) {
    # Tabela estratificada
    summary_data <- data %>%
      group_by(!!sym(group_var)) %>%
      summarise(
        across(all_of(variables), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    q25 = ~quantile(.x, 0.25, na.rm = TRUE),
                    q75 = ~quantile(.x, 0.75, na.rm = TRUE)),
               .names = "{.col}_{.fn}"),
        .groups = "drop"
      )
  } else {
    # Tabela simples
    summary_data <- data %>%
      summarise(
        across(all_of(variables), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    q25 = ~quantile(.x, 0.25, na.rm = TRUE),
                    q75 = ~quantile(.x, 0.75, na.rm = TRUE)),
               .names = "{.col}_{.fn}")
      )
  }
  
  # Criar tabela gt
  gt_table <- summary_data %>%
    gt() %>%
    tab_header(title = title) %>%
    fmt_number(decimals = 2) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) %>%
    tab_options(
      table.font.size = px(12),
      heading.title.font.size = px(14),
      heading.title.font.weight = "bold"
    )
  
  return(gt_table)
}

#' Criar tabela de resultados de teste estatístico
#'
#' Formata resultados de testes estatísticos em tabela profissional
#'
#' @param test_results Lista com resultados de testes
#' @param test_names Nomes dos testes realizados
#' @param title Título da tabela
#' @return Objeto gt formatado
#' @examples
#' # Exemplo com teste t
#' t_result <- t.test(mtcars$mpg ~ mtcars$am)
#' create_test_results_table(list(t_result), "Teste t para MPG")
create_test_results_table <- function(test_results, test_names = NULL, title = "Resultados dos Testes Estatísticos") {
  
  if (is.null(test_names)) {
    test_names <- paste("Teste", seq_along(test_results))
  }
  
  results_df <- data.frame(
    Teste = test_names,
    Estatistica = sapply(test_results, function(x) round(x$statistic, 4)),
    p_valor = sapply(test_results, function(x) 
      ifelse(x$p.value < 0.001, "< 0.001", sprintf("%.3f", x$p.value))),
    IC_inferior = sapply(test_results, function(x) 
      ifelse(!is.null(x$conf.int), round(x$conf.int[1], 4), NA)),
    IC_superior = sapply(test_results, function(x) 
      ifelse(!is.null(x$conf.int), round(x$conf.int[2], 4), NA)),
    stringsAsFactors = FALSE
  )
  
  gt_table <- results_df %>%
    gt() %>%
    tab_header(title = title) %>%
    cols_label(
      Teste = "Teste",
      Estatistica = "Estatística",
      p_valor = "p-valor",
      IC_inferior = "IC 95% (inf)",
      IC_superior = "IC 95% (sup)"
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(columns = p_valor, rows = p_valor == "< 0.001")
    ) %>%
    tab_footnote(
      footnote = "IC = Intervalo de Confiança",
      locations = cells_column_labels(columns = c(IC_inferior, IC_superior))
    )
  
  return(gt_table)
}

#' Criar tabela de correlações
#'
#' Gera matriz de correlação formatada profissionalmente
#'
#' @param data DataFrame com variáveis numéricas
#' @param method Método de correlação ("pearson", "spearman", "kendall")
#' @param show_pvalues Mostrar p-valores (lógico)
#' @param title Título da tabela
#' @return Objeto gt com matriz de correlação
#' @examples
#' create_correlation_table(mtcars[, c("mpg", "hp", "wt", "qsec")])
create_correlation_table <- function(data, method = "pearson", show_pvalues = TRUE, title = "Matriz de Correlação") {
  
  # Calcular correlações
  cor_matrix <- cor(data, use = "complete.obs", method = method)
  
  if (show_pvalues) {
    # Calcular p-valores
    cor_test <- function(x, y) {
      test <- cor.test(x, y, method = method)
      return(test$p.value)
    }
    
    p_matrix <- outer(1:ncol(data), 1:ncol(data), 
                      Vectorize(function(i, j) {
                        if (i == j) return(NA)
                        cor_test(data[,i], data[,j])
                      }))
    
    # Combinar correlação e p-valor
    combined_matrix <- matrix(
      sprintf("%.3f%s", 
              cor_matrix,
              ifelse(is.na(p_matrix) | p_matrix >= 0.05, "", 
                     ifelse(p_matrix < 0.001, "***",
                            ifelse(p_matrix < 0.01, "**", "*")))),
      nrow = nrow(cor_matrix),
      dimnames = list(rownames(cor_matrix), colnames(cor_matrix))
    )
  } else {
    combined_matrix <- round(cor_matrix, 3)
  }
  
  # Converter para data frame
  cor_df <- as.data.frame(combined_matrix)
  cor_df$Variable <- rownames(cor_df)
  cor_df <- cor_df[, c("Variable", setdiff(names(cor_df), "Variable"))]
  
  # Criar tabela gt
  gt_table <- cor_df %>%
    gt() %>%
    tab_header(title = title) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    )
  
  if (show_pvalues) {
    gt_table <- gt_table %>%
      tab_footnote(
        footnote = "* p < 0.05, ** p < 0.01, *** p < 0.001",
        locations = cells_title()
      )
  }
  
  return(gt_table)
}

#' Criar tabela de regressão
#'
#' Formata resultados de modelos de regressão
#'
#' @param model Modelo de regressão (lm, glm, etc.)
#' @param title Título da tabela
#' @param show_ci Mostrar intervalos de confiança
#' @return Objeto gt formatado
#' @examples
#' model <- lm(mpg ~ hp + wt, data = mtcars)
#' create_regression_table(model, "Regressão Linear: MPG")
create_regression_table <- function(model, title = "Resultados da Regressão", show_ci = TRUE) {
  
  # Extrair coeficientes
  coef_summary <- summary(model)$coefficients
  
  if (show_ci) {
    ci <- confint(model)
    results_df <- data.frame(
      Variavel = rownames(coef_summary),
      Coeficiente = round(coef_summary[, "Estimate"], 4),
      Erro_Padrao = round(coef_summary[, "Std. Error"], 4),
      IC_inf = round(ci[, 1], 4),
      IC_sup = round(ci[, 2], 4),
      t_valor = round(coef_summary[, "t value"], 3),
      p_valor = ifelse(coef_summary[, "Pr(>|t|)"] < 0.001, 
                       "< 0.001", 
                       sprintf("%.3f", coef_summary[, "Pr(>|t|)"])),
      stringsAsFactors = FALSE
    )
  } else {
    results_df <- data.frame(
      Variavel = rownames(coef_summary),
      Coeficiente = round(coef_summary[, "Estimate"], 4),
      Erro_Padrao = round(coef_summary[, "Std. Error"], 4),
      t_valor = round(coef_summary[, "t value"], 3),
      p_valor = ifelse(coef_summary[, "Pr(>|t|)"] < 0.001, 
                       "< 0.001", 
                       sprintf("%.3f", coef_summary[, "Pr(>|t|)"])),
      stringsAsFactors = FALSE
    )
  }
  
  # Criar tabela gt
  col_labels <- list(
    Variavel = "Variável",
    Coeficiente = "Coeficiente",
    Erro_Padrao = "Erro Padrão",
    t_valor = "t",
    p_valor = "p-valor"
  )
  
  if (show_ci) {
    col_labels$IC_inf <- "IC 95% (inf)"
    col_labels$IC_sup <- "IC 95% (sup)"
  }
  
  gt_table <- results_df %>%
    gt() %>%
    tab_header(
      title = title,
      subtitle = sprintf("R² = %.3f, R² ajustado = %.3f", 
                        summary(model)$r.squared, 
                        summary(model)$adj.r.squared)
    ) %>%
    cols_label(.list = col_labels) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(columns = p_valor, rows = p_valor == "< 0.001")
    )
  
  return(gt_table)
}

# Exemplo de uso das funções
# Descomente para testar:
# data(mtcars)
# create_descriptive_table(mtcars, c("mpg", "hp", "wt"))
# create_correlation_table(mtcars[, c("mpg", "hp", "wt", "qsec")])
# model <- lm(mpg ~ hp + wt, data = mtcars)
# create_regression_table(model)
