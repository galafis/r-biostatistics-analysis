# mixed_models.R
# Modelos lineares mistos para dados longitudinais em pesquisa clínica/biomédica
# Pacotes: lme4, lmerTest, nlme, ggplot2
# Exemplo reprodutível ao final

#' Ajusta um modelo linear misto (efeitos aleatórios por sujeito)
#' @param data data.frame com dados longitudinais em formato longo
#' @param outcome string com nome do desfecho contínuo (ex.: biomarcador)
#' @param time string com nome da variável de tempo (numérica ou fator ordenado)
#' @param subject string com ID do participante/paciente
#' @param fixed_terms vetor de termos fixos adicionais (ex.: c("treatment","age","sex"))
#' @param random_slope lógico; se TRUE inclui inclinação aleatória de tempo
#' @return lista com: modelo (lmer), ANOVA de efeitos fixos, ICs, e função de predição
#' @examples
#' # Simulação de exemplo clínico
#' set.seed(1)
#' n <- 120; m <- 4
#' id <- rep(1:n, each=m)
#' time <- rep(0:3, times=n)
#' trt <- rep(sample(c("Placebo","Ativo"), n, TRUE), each=m)
#' u0 <- rnorm(n, 0, 1); u1 <- rnorm(n, 0, 0.4)
#' y <- 10 + 0.5*time + 0.8*(trt=="Ativo") + 0.6*time*(trt=="Ativo") +
#'      rep(u0, each=m) + rep(u1, each=m)*time + rnorm(n*m, 0, 1)
#' df <- data.frame(id, time, trt, y)
#' fit <- mixed_lmm(df, outcome="y", time="time", subject="id", fixed_terms=c("trt"), random_slope=TRUE)
#' summary(fit$model)
#' fit$anova
#' predict_df <- data.frame(id=1, time=0:3, trt="Ativo")
#' fit$predict_fn(predict_df)
mixed_lmm <- function(data, outcome, time, subject, fixed_terms=NULL, random_slope=TRUE){
  stopifnot(is.data.frame(data), all(c(outcome,time,subject) %in% names(data)))
  if (!is.null(fixed_terms)) stopifnot(all(fixed_terms %in% names(data)))
  if (!requireNamespace("lme4", quietly=TRUE)) stop("Package 'lme4' is required")
  if (!requireNamespace("lmerTest", quietly=TRUE)) warning("Package 'lmerTest' not found; p-values may be unavailable")
  if (!requireNamespace("ggplot2", quietly=TRUE)) warning("Package 'ggplot2' suggested for plotting")

  # construir fórmula de efeitos fixos
  fixed <- c(time, fixed_terms)
  fixed <- fixed[!is.na(fixed) & nzchar(fixed)]
  rhs <- if (length(fixed)) paste(fixed, collapse=" + ") else "1"

  # estrutura aleatória
  rand <- if (isTRUE(random_slope)) sprintf("(1 + %s | %s)", time, subject) else sprintf("(1 | %s)", subject)
  fml <- stats::as.formula(sprintf("%s ~ %s + %s", outcome, rhs, rand))

  # ajustar
  model <- lme4::lmer(fml, data=data, REML=TRUE)
  an <- tryCatch({
    if (requireNamespace("lmerTest", quietly=TRUE)) lmerTest::anova(model) else stats::anova(model)
  }, error=function(e) stats::anova(model))
  ci <- suppressMessages(stats::confint(model, method="Wald"))

  predict_fn <- function(newdata, re.form=NA){ # por padrão, efeitos populacionais
    stats::predict(model, newdata=newdata, re.form=re.form, allow.new.levels=TRUE)
  }

  list(model=model, anova=an, confint=ci, formula=fml, predict_fn=predict_fn)
}

#' Ajusta um modelo misto com correlação/variância no tempo pelo nlme
#' Útil para estruturas como AR(1), Compound Symmetry, etc.
#' @param corr estrutura de correlação: "ar1", "cs" ou NULL
#' @param var_func heterocedasticidade via varPower/varIdent (string opcional)
#' @examples
#' # nlme com AR(1)
#' # fit2 <- mixed_nlme(df, outcome="y", time="time", subject="id", fixed_terms=c("trt"), corr="ar1")
#' mixed_nlme <- function(data, outcome, time, subject, fixed_terms=NULL, corr=c("ar1","cs", "none"), var_func=NULL){
  corr <- match.arg(corr)
  stopifnot(is.data.frame(data), all(c(outcome,time,subject) %in% names(data)))
  if (!is.null(fixed_terms)) stopifnot(all(fixed_terms %in% names(data)))
  if (!requireNamespace("nlme", quietly=TRUE)) stop("Package 'nlme' is required")

  fixed <- c(time, fixed_terms)
  fixed <- fixed[!is.na(fixed) & nzchar(fixed)]
  rhs <- if (length(fixed)) paste(fixed, collapse=" + ") else "1"
  fml_fixed <- stats::as.formula(sprintf("%s ~ %s", outcome, rhs))
  rand <- stats::as.formula(sprintf("~ 1 + %s | %s", time, subject))

  corStruct <- switch(corr,
    ar1 = nlme::corAR1(form=stats::as.formula(sprintf("~ %s | %s", time, subject))),
    cs  = nlme::corCompSymm(form=stats::as.formula(sprintf("~ %s | %s", time, subject))),
    none = NULL
  )
  varStruct <- NULL
  if (!is.null(var_func)){
    if (grepl("power", var_func, ignore.case=TRUE)) varStruct <- nlme::varPower()
    if (grepl("ident", var_func, ignore.case=TRUE)) varStruct <- nlme::varIdent(form=stats::as.formula(sprintf("~1|%s", time)))
  }

  model <- nlme::lme(fixed=fml_fixed, random=rand, data=data, method="REML",
                      correlation=corStruct, weights=varStruct, na.action=na.omit)
  ci <- tryCatch(nlme::intervals(model), error=function(e) NULL)
  list(model=model, intervals=ci)
}

#' Gráfico de perfis e médias ajustadas por tratamento ao longo do tempo
#' @examples
#' # plot_profiles(df, outcome="y", time="time", group="trt")
plot_profiles <- function(data, outcome, time, group=NULL){
  if (!requireNamespace("ggplot2", quietly=TRUE)) stop("Package 'ggplot2' is required for plotting")
  library(ggplot2)
  p <- ggplot(data, aes_string(x=time, y=outcome, color=group, group=paste0("interaction(", group, ", ", "id)") )) +
    geom_line(alpha=0.3) +
    stat_summary(fun=mean, geom="line", size=1.2, aes(group=group)) +
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.1) +
    theme_minimal() + labs(x="Tempo", y=outcome, color=group)
  p
}
