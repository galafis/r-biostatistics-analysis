# 🇧🇷 Análises Bioestatísticas Avançadas | 🇺🇸 Advanced Biostatistical Analysis

<div align="center">

![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Biostatistics](https://img.shields.io/badge/Biostatistics-Medical%20Research-red?style=for-the-badge)
![Clinical](https://img.shields.io/badge/Clinical-Trials-green?style=for-the-badge)
![License](https://img.shields.io/badge/License-MIT-green.svg?style=for-the-badge)

**Plataforma completa para análises bioestatísticas e pesquisa médica**

[🧬 Análises](#-análises-disponíveis) • [📊 Métodos](#-métodos-estatísticos) • [⚡ Instalação](#-instalação) • [🏥 Aplicações](#-aplicações-clínicas)

</div>

---

## 🇧🇷 Português

### 🧬 Visão Geral

Plataforma abrangente de **análises bioestatísticas** desenvolvida em R, especializada em:

- 🏥 **Ensaios Clínicos**: Delineamento, análise e interpretação
- 📊 **Análise de Sobrevivência**: Kaplan-Meier, Cox, modelos paramétricos
- 🧬 **Epidemiologia**: Estudos observacionais, meta-análises
- 📈 **Biometria**: Análise de dados biomédicos complexos
- 📋 **Relatórios Regulatórios**: Conformidade com ICH, FDA, EMA

### 🎯 Objetivos da Plataforma

- **Facilitar análises** estatísticas em pesquisa médica
- **Padronizar métodos** bioestatísticos
- **Automatizar relatórios** para agências regulatórias
- **Validar resultados** com métodos robustos
- **Educar pesquisadores** em boas práticas

### 🛠️ Stack Tecnológico

#### Análise de Sobrevivência
- **survival**: Análise de sobrevivência clássica
- **survminer**: Visualização de curvas de sobrevivência
- **flexsurv**: Modelos paramétricos flexíveis
- **survivalROC**: Curvas ROC dependentes do tempo

#### Meta-análise
- **meta**: Meta-análises padrão
- **metafor**: Meta-análises avançadas
- **netmeta**: Meta-análises em rede
- **forestplot**: Gráficos forest profissionais

#### Ensaios Clínicos
- **gsDesign**: Delineamento de estudos sequenciais
- **rpact**: Análises adaptativas
- **PowerTOST**: Cálculos de poder para bioequivalência
- **Hmisc**: Utilitários para análise clínica

#### Epidemiologia
- **epiR**: Análises epidemiológicas
- **epitools**: Ferramentas epidemiológicas
- **tableone**: Tabelas descritivas
- **MatchIt**: Pareamento de casos

#### Visualização Médica
- **ggplot2**: Gráficos estatísticos
- **ggpubr**: Gráficos para publicação
- **survminer**: Visualização de sobrevivência
- **forestplot**: Gráficos forest

### 📋 Estrutura da Plataforma

```
r-biostatistics-analysis/
├── 📁 clinical_trials/            # Análises de ensaios clínicos
│   ├── 📄 rct_analysis.R         # Ensaios randomizados
│   ├── 📄 adaptive_designs.R     # Delineamentos adaptativos
│   ├── 📄 bioequivalence.R       # Estudos de bioequivalência
│   └── 📄 dose_finding.R         # Estudos dose-resposta
├── 📁 survival_analysis/          # Análise de sobrevivência
│   ├── 📄 kaplan_meier.R         # Estimador Kaplan-Meier
│   ├── 📄 cox_regression.R       # Regressão de Cox
│   ├── 📄 parametric_models.R    # Modelos paramétricos
│   ├── 📄 competing_risks.R      # Riscos competitivos
│   └── 📄 time_varying_effects.R # Efeitos tempo-dependentes
├── 📁 meta_analysis/              # Meta-análises
│   ├── 📄 fixed_effects.R        # Efeitos fixos
│   ├── 📄 random_effects.R       # Efeitos aleatórios
│   ├── 📄 network_meta.R         # Meta-análise em rede
│   ├── 📄 publication_bias.R     # Viés de publicação
│   └── 📄 sensitivity_analysis.R # Análise de sensibilidade
├── 📁 epidemiology/               # Estudos epidemiológicos
│   ├── 📄 cohort_studies.R       # Estudos de coorte
│   ├── 📄 case_control.R         # Caso-controle
│   ├── 📄 cross_sectional.R      # Estudos transversais
│   ├── 📄 propensity_score.R     # Escore de propensão
│   └── 📄 causal_inference.R     # Inferência causal
├── 📁 diagnostic_tests/           # Testes diagnósticos
│   ├── 📄 roc_analysis.R         # Análise ROC
│   ├── 📄 diagnostic_accuracy.R  # Acurácia diagnóstica
│   ├── 📄 agreement_studies.R    # Estudos de concordância
│   └── 📄 screening_tests.R      # Testes de rastreamento
├── 📁 longitudinal/               # Dados longitudinais
│   ├── 📄 mixed_models.R         # Modelos mistos
│   ├── 📄 gee_analysis.R         # Equações de estimação generalizadas
│   ├── 📄 growth_curves.R        # Curvas de crescimento
│   └── 📄 repeated_measures.R    # Medidas repetidas
├── 📁 genomics/                   # Análises genômicas
│   ├── 📄 gwas_analysis.R        # Estudos de associação genômica
│   ├── 📄 linkage_analysis.R     # Análise de ligação
│   ├── 📄 population_genetics.R  # Genética populacional
│   └── 📄 pharmacogenomics.R     # Farmacogenômica
├── 📁 reports/                    # Templates de relatórios
│   ├── 📄 clinical_study_report.Rmd # Relatório de estudo clínico
│   ├── 📄 statistical_analysis_plan.Rmd # Plano de análise estatística
│   ├── 📄 interim_analysis.Rmd   # Análise interina
│   └── 📄 final_report.Rmd       # Relatório final
├── 📁 data/                       # Datasets médicos
│   ├── 📁 clinical_trials/       # Dados de ensaios
│   ├── 📁 survival/              # Dados de sobrevivência
│   ├── 📁 epidemiological/       # Dados epidemiológicos
│   └── 📄 data_dictionary.md     # Dicionário de dados
├── 📁 functions/                  # Funções customizadas
│   ├── 📄 biostat_utils.R        # Utilitários bioestatísticos
│   ├── 📄 plotting_functions.R   # Funções de visualização
│   ├── 📄 table_functions.R      # Funções para tabelas
│   └── 📄 validation_functions.R # Funções de validação
├── 📁 validation/                 # Validação de métodos
│   ├── 📄 method_validation.R    # Validação de métodos
│   ├── 📄 simulation_studies.R   # Estudos de simulação
│   └── 📄 benchmark_tests.R      # Testes de benchmark
├── 📄 README.md                  # Este arquivo
├── 📄 LICENSE                    # Licença MIT
├── 📄 .gitignore                # Arquivos ignorados
└── 📄 renv.lock                 # Controle de dependências
```

### 🧬 Análises Disponíveis

#### 1. 🏥 Ensaios Clínicos Randomizados

**Análise de Eficácia**
```r
# Análise de eficácia primária
efficacy_analysis <- function(data, endpoint, treatment, covariates = NULL) {
  if (is.null(covariates)) {
    # Análise não ajustada
    model <- lm(get(endpoint) ~ get(treatment), data = data)
  } else {
    # Análise ajustada
    formula_str <- paste(endpoint, "~", treatment, "+", paste(covariates, collapse = " + "))
    model <- lm(as.formula(formula_str), data = data)
  }
  
  # Resultados
  list(
    model = model,
    treatment_effect = summary(model)$coefficients[2, 1],
    ci_lower = confint(model)[2, 1],
    ci_upper = confint(model)[2, 2],
    p_value = summary(model)$coefficients[2, 4]
  )
}
```

**Análise de Não-Inferioridade**
```r
# Teste de não-inferioridade
non_inferiority_test <- function(treatment_diff, margin, alpha = 0.025) {
  # H0: treatment_diff <= -margin (inferior)
  # H1: treatment_diff > -margin (não-inferior)
  
  t_stat <- (treatment_diff + margin) / se_diff
  p_value <- 1 - pt(t_stat, df)
  
  list(
    conclusion = ifelse(p_value < alpha, "Não-inferior", "Inconclusivo"),
    p_value = p_value,
    margin = margin
  )
}
```

#### 2. 📊 Análise de Sobrevivência

**Estimador Kaplan-Meier**
```r
# Análise de sobrevivência Kaplan-Meier
km_analysis <- function(time, event, group = NULL) {
  library(survival)
  library(survminer)
  
  if (is.null(group)) {
    # Sobrevivência geral
    fit <- survfit(Surv(time, event) ~ 1)
  } else {
    # Sobrevivência por grupo
    fit <- survfit(Surv(time, event) ~ group)
  }
  
  # Gráfico de sobrevivência
  plot <- ggsurvplot(
    fit,
    data = data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    ncensor.plot = TRUE
  )
  
  list(fit = fit, plot = plot)
}
```

**Regressão de Cox**
```r
# Modelo de Cox multivariado
cox_analysis <- function(data, time_var, event_var, covariates) {
  library(survival)
  
  # Fórmula do modelo
  formula_str <- paste("Surv(", time_var, ",", event_var, ") ~", 
                      paste(covariates, collapse = " + "))
  
  # Ajustar modelo
  cox_model <- coxph(as.formula(formula_str), data = data)
  
  # Teste de proporcionalidade
  ph_test <- cox.zph(cox_model)
  
  # Hazard ratios
  hr_table <- exp(cbind(
    HR = coef(cox_model),
    confint(cox_model)
  ))
  
  list(
    model = cox_model,
    hazard_ratios = hr_table,
    proportionality_test = ph_test
  )
}
```

#### 3. 🔬 Meta-análise

**Meta-análise de Efeitos Aleatórios**
```r
# Meta-análise com efeitos aleatórios
random_effects_meta <- function(effect_sizes, variances, study_names) {
  library(metafor)
  
  # Modelo de efeitos aleatórios
  meta_model <- rma(yi = effect_sizes, vi = variances, 
                   slab = study_names, method = "REML")
  
  # Forest plot
  forest_plot <- forest(meta_model,
                       showweights = TRUE,
                       header = TRUE)
  
  # Teste de heterogeneidade
  heterogeneity <- list(
    Q = meta_model$QE,
    p_value = meta_model$QEp,
    I2 = meta_model$I2,
    tau2 = meta_model$tau2
  )
  
  list(
    model = meta_model,
    heterogeneity = heterogeneity,
    plot = forest_plot
  )
}
```

**Análise de Viés de Publicação**
```r
# Teste de viés de publicação
publication_bias_test <- function(meta_model) {
  # Teste de Egger
  egger_test <- regtest(meta_model)
  
  # Funnel plot
  funnel_plot <- funnel(meta_model, main = "Funnel Plot")
  
  # Trim and fill
  trimfill <- trimfill(meta_model)
  
  list(
    egger_test = egger_test,
    funnel_plot = funnel_plot,
    trimfill = trimfill
  )
}
```

#### 4. 🧪 Testes Diagnósticos

**Análise ROC**
```r
# Análise de curva ROC
roc_analysis <- function(predictor, outcome) {
  library(pROC)
  
  # Curva ROC
  roc_curve <- roc(outcome, predictor)
  
  # Área sob a curva
  auc_value <- auc(roc_curve)
  auc_ci <- ci.auc(roc_curve)
  
  # Ponto de corte ótimo (Youden)
  optimal_cutoff <- coords(roc_curve, "best", ret = "threshold")
  
  # Sensibilidade e especificidade no ponto ótimo
  sens_spec <- coords(roc_curve, optimal_cutoff, ret = c("sensitivity", "specificity"))
  
  list(
    roc_curve = roc_curve,
    auc = auc_value,
    auc_ci = auc_ci,
    optimal_cutoff = optimal_cutoff,
    sensitivity = sens_spec$sensitivity,
    specificity = sens_spec$specificity
  )
}
```

#### 5. 📈 Modelos Longitudinais

**Modelos Lineares Mistos**
```r
# Modelo linear misto para dados longitudinais
mixed_model_analysis <- function(data, outcome, time, subject, treatment, covariates = NULL) {
  library(lme4)
  library(lmerTest)
  
  # Construir fórmula
  fixed_effects <- paste(c(time, treatment, paste0(time, "*", treatment), covariates), 
                        collapse = " + ")
  formula_str <- paste(outcome, "~", fixed_effects, "+ (", time, "|", subject, ")")
  
  # Ajustar modelo
  mixed_model <- lmer(as.formula(formula_str), data = data)
  
  # Testes de efeitos fixos
  fixed_effects_test <- anova(mixed_model)
  
  # Gráfico de perfis
  profile_plot <- ggplot(data, aes_string(x = time, y = outcome, color = treatment)) +
    geom_smooth(method = "loess", se = TRUE) +
    geom_point(alpha = 0.3) +
    theme_minimal()
  
  list(
    model = mixed_model,
    fixed_effects = fixed_effects_test,
    plot = profile_plot
  )
}
```

### 🏥 Aplicações Clínicas

#### 1. Oncologia
- **Análise de Sobrevivência**: Sobrevida global, livre de progressão
- **Ensaios Fase I**: Estudos de dose-escalação
- **Biomarcadores**: Análise de marcadores prognósticos
- **Qualidade de Vida**: Análise de desfechos reportados pelo paciente

#### 2. Cardiologia
- **Estudos de Desfechos Cardiovasculares**: MACE, mortalidade
- **Análise de Fatores de Risco**: Modelos preditivos
- **Ensaios de Prevenção**: Análise de eventos raros
- **Estudos de Imagem**: Análise de medidas repetidas

#### 3. Neurologia
- **Estudos de Demência**: Análise de declínio cognitivo
- **Ensaios em Esclerose Múltipla**: Análise de recidivas
- **Estudos de AVC**: Análise de desfechos funcionais
- **Biomarcadores Neurológicos**: Análise de biomarcadores

#### 4. Infectologia
- **Ensaios de Vacinas**: Análise de imunogenicidade
- **Estudos de Resistência**: Análise de mutações
- **Farmacocinética**: Modelos PK/PD
- **Estudos Epidemiológicos**: Análise de surtos

### 🎯 Competências Demonstradas

#### Bioestatística Clínica
- ✅ **Delineamento de Estudos**: RCT, observacionais, adaptativos
- ✅ **Análise de Sobrevivência**: Kaplan-Meier, Cox, paramétricos
- ✅ **Meta-análise**: Efeitos fixos/aleatórios, rede, viés
- ✅ **Dados Longitudinais**: Modelos mistos, GEE, curvas de crescimento

#### Métodos Estatísticos
- ✅ **Inferência Causal**: Escore de propensão, variáveis instrumentais
- ✅ **Testes Diagnósticos**: ROC, acurácia, concordância
- ✅ **Análise Multivariada**: Regressão múltipla, logística, Poisson
- ✅ **Métodos Bayesianos**: Análise bayesiana, MCMC

#### Regulatório e Compliance
- ✅ **ICH Guidelines**: E9, E10, E6 (GCP)
- ✅ **FDA Guidance**: Adaptive designs, missing data
- ✅ **EMA Guidelines**: Biostatistical methodology
- ✅ **Validação**: Validação de software estatístico

### 📊 Exemplos de Relatórios

#### Relatório de Estudo Clínico (CSR)
```r
# Template para CSR
csr_template <- function(study_data, primary_endpoint, treatment_var) {
  # Análise demográfica
  demographics <- tableone::CreateTableOne(
    vars = demographic_vars,
    strata = treatment_var,
    data = study_data
  )
  
  # Análise de eficácia primária
  primary_analysis <- efficacy_analysis(
    data = study_data,
    endpoint = primary_endpoint,
    treatment = treatment_var
  )
  
  # Análise de segurança
  safety_analysis <- safety_summary(study_data)
  
  # Gerar relatório
  rmarkdown::render("templates/csr_template.Rmd",
                   params = list(
                     demographics = demographics,
                     efficacy = primary_analysis,
                     safety = safety_analysis
                   ))
}
```

### 🔧 Configuração e Validação

#### Ambiente Validado
```r
# Configuração de ambiente validado
setup_validated_environment <- function() {
  # Versões específicas de pacotes
  required_packages <- c(
    "survival@3.2-13",
    "meta@5.2-0",
    "metafor@3.4-0",
    "lme4@1.1-29"
  )
  
  # Instalar versões específicas
  devtools::install_version(required_packages)
  
  # Documentar ambiente
  sessionInfo()
}
```

#### Validação de Resultados
```r
# Validação cruzada de resultados
validate_results <- function(analysis_function, test_data, reference_software) {
  # Executar análise em R
  r_results <- analysis_function(test_data)
  
  # Comparar com software de referência
  comparison <- compare_results(r_results, reference_software)
  
  # Documentar diferenças
  validation_report <- generate_validation_report(comparison)
  
  return(validation_report)
}
```

---

## 🇺🇸 English

### 🧬 Overview

Comprehensive **biostatistical analysis** platform developed in R, specialized in:

- 🏥 **Clinical Trials**: Design, analysis, and interpretation
- 📊 **Survival Analysis**: Kaplan-Meier, Cox, parametric models
- 🧬 **Epidemiology**: Observational studies, meta-analyses
- 📈 **Biometrics**: Complex biomedical data analysis
- 📋 **Regulatory Reports**: ICH, FDA, EMA compliance

### 🎯 Platform Objectives

- **Facilitate statistical analyses** in medical research
- **Standardize biostatistical methods**
- **Automate reports** for regulatory agencies
- **Validate results** with robust methods
- **Educate researchers** in best practices

### 🧬 Available Analyses

#### 1. 🏥 Randomized Clinical Trials
- Efficacy analysis (primary and secondary endpoints)
- Safety analysis and adverse events
- Non-inferiority and superiority testing
- Interim analysis and adaptive designs

#### 2. 📊 Survival Analysis
- Kaplan-Meier estimation
- Cox proportional hazards regression
- Parametric survival models
- Competing risks analysis

#### 3. 🔬 Meta-analysis
- Fixed and random effects models
- Network meta-analysis
- Publication bias assessment
- Sensitivity analysis

#### 4. 🧪 Diagnostic Tests
- ROC curve analysis
- Diagnostic accuracy measures
- Agreement studies
- Screening test evaluation

#### 5. 📈 Longitudinal Models
- Linear mixed models
- Generalized estimating equations
- Growth curve analysis
- Repeated measures ANOVA

### 🎯 Skills Demonstrated

#### Clinical Biostatistics
- ✅ **Study Design**: RCT, observational, adaptive
- ✅ **Survival Analysis**: Kaplan-Meier, Cox, parametric
- ✅ **Meta-analysis**: Fixed/random effects, network, bias
- ✅ **Longitudinal Data**: Mixed models, GEE, growth curves

#### Statistical Methods
- ✅ **Causal Inference**: Propensity score, instrumental variables
- ✅ **Diagnostic Tests**: ROC, accuracy, agreement
- ✅ **Multivariate Analysis**: Multiple, logistic, Poisson regression
- ✅ **Bayesian Methods**: Bayesian analysis, MCMC

#### Regulatory and Compliance
- ✅ **ICH Guidelines**: E9, E10, E6 (GCP)
- ✅ **FDA Guidance**: Adaptive designs, missing data
- ✅ **EMA Guidelines**: Biostatistical methodology
- ✅ **Validation**: Statistical software validation

---

## 📄 Licença | License

MIT License - veja o arquivo [LICENSE](LICENSE) para detalhes | see [LICENSE](LICENSE) file for details

## 📞 Contato | Contact

**GitHub**: [@galafis](https://github.com/galafis)  
**LinkedIn**: [Gabriel Demetrios Lafis](https://linkedin.com/in/galafis)  
**Email**: gabriel.lafis@example.com

---

<div align="center">

**Desenvolvido com ❤️ para Pesquisa Médica | Developed with ❤️ for Medical Research**

[![GitHub](https://img.shields.io/badge/GitHub-galafis-blue?style=flat-square&logo=github)](https://github.com/galafis)
[![R](https://img.shields.io/badge/R-276DC3?style=flat-square&logo=r&logoColor=white)](https://www.r-project.org/)

</div>

