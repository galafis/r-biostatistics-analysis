# ============================================================================
# PLOTTING FUNCTIONS
# Funções auxiliares para visualização de dados e análise estatística
# Projeto: R Biostatistics Analysis
# Autor: [Nome do Autor]
# Data: 2024
# ============================================================================

# Dependências necessárias
require(ggplot2)
require(dplyr)
require(scales)
require(gridExtra)
require(RColorBrewer)

#' Histograma com densidade e estatísticas descritivas
#' 
#' @param data Data frame
#' @param variable Nome da variável (string)
#' @param bins Número de bins (default: 30)
#' @param title Título do gráfico
#' @param show_stats Mostrar estatísticas no gráfico (default: TRUE)
#' @return ggplot object
plot_histogram <- function(data, variable, bins = 30, title = NULL, show_stats = TRUE) {
  if (!variable %in% names(data)) {
    stop(paste("Variável", variable, "não encontrada no dataset"))
  }
  
  # Remover valores missing
  clean_data <- data[!is.na(data[[variable]]), ]
  
  # Calcular estatísticas
  stats <- list(
    mean = mean(clean_data[[variable]]),
    median = median(clean_data[[variable]]),
    sd = sd(clean_data[[variable]])
  )
  
  p <- ggplot(clean_data, aes_string(x = variable)) +
    geom_histogram(aes(y = ..density..), bins = bins, fill = "lightblue", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    geom_vline(xintercept = stats$mean, color = "blue", linetype = "dashed", size = 1) +
    geom_vline(xintercept = stats$median, color = "green", linetype = "dashed", size = 1) +
    labs(title = title %||% paste("Distribuição de", variable),
         x = variable,
         y = "Densidade") +
    theme_minimal()
  
  if (show_stats) {
    stats_text <- sprintf("Média: %.2f\nMediana: %.2f\nDP: %.2f", 
                         stats$mean, stats$median, stats$sd)
    p <- p + annotate("text", x = Inf, y = Inf, label = stats_text, 
                     hjust = 1.1, vjust = 1.1, size = 3)
  }
  
  return(p)
}

#' Boxplot comparativo entre grupos
#' 
#' @param data Data frame
#' @param x_var Variável categórica (string)
#' @param y_var Variável numérica (string)
#' @param title Título do gráfico
#' @param show_points Mostrar pontos individuais (default: TRUE)
#' @return ggplot object
plot_boxplot_comparison <- function(data, x_var, y_var, title = NULL, show_points = TRUE) {
  if (!all(c(x_var, y_var) %in% names(data))) {
    stop("Uma ou mais variáveis não encontradas no dataset")
  }
  
  p <- ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    labs(title = title %||% paste(y_var, "por", x_var),
         x = x_var,
         y = y_var) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (show_points) {
    p <- p + geom_jitter(width = 0.2, alpha = 0.5, color = "red")
  }
  
  return(p)
}

#' Scatter plot com linha de regressão
#' 
#' @param data Data frame
#' @param x_var Variável x (string)
#' @param y_var Variável y (string)
#' @param color_var Variável para colorir pontos (opcional)
#' @param title Título do gráfico
#' @param add_regression Adicionar linha de regressão (default: TRUE)
#' @return ggplot object
plot_scatter <- function(data, x_var, y_var, color_var = NULL, title = NULL, add_regression = TRUE) {
  if (!all(c(x_var, y_var) %in% names(data))) {
    stop("Uma ou mais variáveis não encontradas no dataset")
  }
  
  aes_mapping <- aes_string(x = x_var, y = y_var)
  if (!is.null(color_var)) {
    aes_mapping$colour <- as.name(color_var)
  }
  
  p <- ggplot(data, aes_mapping) +
    geom_point(alpha = 0.6) +
    labs(title = title %||% paste(y_var, "vs", x_var),
         x = x_var,
         y = y_var) +
    theme_minimal()
  
  if (add_regression) {
    p <- p + geom_smooth(method = "lm", se = TRUE, color = "red")
  }
  
  return(p)
}

#' Gráfico de barras para frequências
#' 
#' @param data Data frame
#' @param variable Variável categórica (string)
#' @param title Título do gráfico
#' @param show_percentages Mostrar percentuais (default: TRUE)
#' @return ggplot object
plot_frequency_bar <- function(data, variable, title = NULL, show_percentages = TRUE) {
  if (!variable %in% names(data)) {
    stop(paste("Variável", variable, "não encontrada no dataset"))
  }
  
  # Calcular frequências
  freq_data <- data %>%
    count(!!sym(variable)) %>%
    mutate(percentage = n / sum(n) * 100)
  
  p <- ggplot(freq_data, aes_string(x = variable, y = "n")) +
    geom_bar(stat = "identity", fill = "lightblue", alpha = 0.7) +
    labs(title = title %||% paste("Frequências de", variable),
         x = variable,
         y = "Frequência") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (show_percentages) {
    p <- p + geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")),
                      vjust = -0.5, size = 3)
  }
  
  return(p)
}

#' QQ-plot para verificação de normalidade
#' 
#' @param data Data frame
#' @param variable Variável numérica (string)
#' @param title Título do gráfico
#' @return ggplot object
plot_qq <- function(data, variable, title = NULL) {
  if (!variable %in% names(data)) {
    stop(paste("Variável", variable, "não encontrada no dataset"))
  }
  
  clean_data <- data[!is.na(data[[variable]]), ]
  
  p <- ggplot(clean_data, aes_string(sample = variable)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    labs(title = title %||% paste("QQ-Plot:", variable),
         x = "Quantis Teóricos",
         y = "Quantis Observados") +
    theme_minimal()
  
  return(p)
}

#' Heatmap de correlação
#' 
#' @param data Data frame com variáveis numéricas
#' @param title Título do gráfico
#' @param method Método de correlação ("pearson", "spearman", "kendall")
#' @return ggplot object
plot_correlation_heatmap <- function(data, title = NULL, method = "pearson") {
  # Selecionar apenas variáveis numéricas
  numeric_data <- data %>% select_if(is.numeric)
  
  if (ncol(numeric_data) < 2) {
    stop("Pelo menos duas variáveis numéricas são necessárias")
  }
  
  # Calcular matriz de correlação
  cor_matrix <- cor(numeric_data, use = "complete.obs", method = method)
  
  # Transformar em formato longo
  cor_data <- as.data.frame(as.table(cor_matrix))
  names(cor_data) <- c("Var1", "Var2", "Correlation")
  
  p <- ggplot(cor_data, aes(Var1, Var2, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = 0, limit = c(-1, 1)) +
    geom_text(aes(label = round(Correlation, 2)), size = 3) +
    labs(title = title %||% "Matriz de Correlação",
         x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Gráfico de sobrevivência (Kaplan-Meier)
#' 
#' @param data Data frame
#' @param time_var Variável de tempo (string)
#' @param event_var Variável de evento (string)
#' @param group_var Variável de agrupamento (opcional)
#' @param title Título do gráfico
#' @return ggplot object
plot_survival <- function(data, time_var, event_var, group_var = NULL, title = NULL) {
  require(survival)
  require(survminer)
  
  if (!all(c(time_var, event_var) %in% names(data))) {
    stop("Uma ou mais variáveis não encontradas no dataset")
  }
  
  if (is.null(group_var)) {
    surv_formula <- as.formula(paste("Surv(", time_var, ",", event_var, ") ~ 1"))
  } else {
    surv_formula <- as.formula(paste("Surv(", time_var, ",", event_var, ") ~", group_var))
  }
  
  surv_fit <- survfit(surv_formula, data = data)
  
  p <- ggsurvplot(surv_fit, data = data,
                  title = title %||% "Curva de Sobrevivência",
                  xlab = "Tempo",
                  ylab = "Probabilidade de Sobrevivência",
                  risk.table = TRUE,
                  conf.int = TRUE)
  
  return(p)
}

#' Salvar múltiplos gráficos em arquivo PDF
#' 
#' @param plots Lista de objetos ggplot
#' @param filename Nome do arquivo (sem extensão)
#' @param width Largura em polegadas (default: 8)
#' @param height Altura em polegadas (default: 6)
save_multiple_plots <- function(plots, filename, width = 8, height = 6) {
  pdf_name <- paste0(filename, ".pdf")
  pdf(pdf_name, width = width, height = height)
  
  for (i in seq_along(plots)) {
    print(plots[[i]])
  }
  
  dev.off()
  message(paste("Gráficos salvos em:", pdf_name))
}

# Operador auxiliar para valores padrão
`%||%` <- function(x, y) if (is.null(x)) y else x

cat("Funções de plotting carregadas com sucesso!\n")
cat("Funções disponíveis:\n")
cat("- plot_histogram()\n")
cat("- plot_boxplot_comparison()\n")
cat("- plot_scatter()\n")
cat("- plot_frequency_bar()\n")
cat("- plot_qq()\n")
cat("- plot_correlation_heatmap()\n")
cat("- plot_survival()\n")
cat("- save_multiple_plots()\n")
