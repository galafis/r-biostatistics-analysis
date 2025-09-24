# Clinical Trials Data | Dados de Ensaios Clínicos

## 🇧🇷 Português | 🇺🇸 English

Esta pasta contém datasets exemplo para análises de ensaios clínicos. | This folder contains example datasets for clinical trial analyses.

## 📊 Datasets Disponíveis | Available Datasets

### 🎆 Ensaios Clínicos Randomizados | Randomized Clinical Trials

#### `rct_data.csv`
Dataset exemplo para ensaios clínicos randomizados controlados com placebo.

**Características | Characteristics:**
- Número de participantes | Number of participants: 200
- Grupos de tratamento | Treatment groups: Placebo vs. Ativo
- Desfecho primário | Primary endpoint: Escore de eficácia (0-100)
- Período de seguimento | Follow-up period: 12 semanas

#### `rct_safety_data.csv`
Dados de segurança e eventos adversos.

**Conteúdo | Content:**
- Eventos adversos por severidade | Adverse events by severity
- Relacionamento com o fármaco | Relationship to study drug
- Ações tomadas | Actions taken
- Desfecho do evento | Event outcome

### ⚙️ Estudos de Bioequivalência | Bioequivalence Studies

#### `bioequivalence_data.csv`
Dataset para estudos de bioequivalência com delineamento cruzado.

**Características | Characteristics:**
- Delineamento | Design: Cruzado 2x2 (2x2 crossover)
- Participantes | Subjects: 24 voluntários sadios
- Produtos | Products: Teste vs. Referência
- Parâmetros PK | PK parameters: Cmax, AUC, Tmax

#### `pk_concentration_data.csv`
Concentrações plasmáticas ao longo do tempo.

**Detalhes | Details:**
- Tempos de coleta | Collection times: 0-48h
- Método analítico | Analytical method: LC-MS/MS
- Limite inferior de quantificação | LLOQ: 0.1 ng/mL

### 🔍 Estudos de Dose-Resposta | Dose-Response Studies

#### `dose_response_data.csv`
Dados para modelagem dose-resposta.

**Características | Characteristics:**
- Doses testadas | Doses tested: 5, 10, 20, 40 mg
- Modelo | Model: Sigmóide Emax
- Biomarcador | Biomarker: Inibição enzimática (%)
- Tempo de avaliação | Assessment time: 2h pós-dose

## 🔌 Estrutura dos Dados | Data Structure

### Variáveis Comuns | Common Variables
- `subject_id`: Identificador único do participante | Unique subject identifier
- `treatment`: Grupo de tratamento | Treatment group  
- `period`: Período do estudo | Study period
- `sequence`: Sequência de tratamento | Treatment sequence
- `baseline_*`: Variáveis basais | Baseline variables
- `endpoint_*`: Variáveis de desfecho | Endpoint variables

### Codificação | Coding
- Missing values: `NA`
- Treatment codes: `"Placebo"`, `"Active"`, `"Test"`, `"Reference"`
- Dates: ISO format `YYYY-MM-DD`
- Time points: Hours from first dose

## 🛠️ Ferramentas de Análise | Analysis Tools

### Pacotes R Recomendados | Recommended R Packages
- **PowerTOST**: Cálculos de poder para bioequivalência
- **nlme/lme4**: Modelos de efeitos mistos
- **survival**: Análise de tempo até o evento
- **gsDesign**: Delineamento de estudos adaptativos
- **rpact**: Análises interinas

### Scripts de Análise | Analysis Scripts
Ver pasta principal `/clinical_trials/` para scripts R completos.

## 📊 Exemplos de Análise | Analysis Examples

### Análise de Eficácia Primária
```r
# Carregar dados
data <- read.csv("rct_data.csv")

# Análise ANCOVA
model <- lm(endpoint_score ~ treatment + baseline_score, data = data)
summary(model)
```

### Análise de Bioequivalência
```r
# Carregar PowerTOST
library(PowerTOST)

# Análise 2x2 crossover
be_result <- be.ABE.CI(data = be_data, 
                       log.scale = TRUE,
                       alpha = 0.05)
```

## 📄 Documentação | Documentation

### Padrões Regulatórios | Regulatory Standards
- ICH E6 (GCP): Boas Práticas Clínicas
- ICH E9: Princípios Estatísticos para Ensaios Clínicos  
- FDA 21 CFR 320: Estudos de Bioequivalência
- EMA CHMP: Guidelines de Bioequivalência

### Validação de Dados | Data Validation
- Consistência entre datasets
- Verificação de outliers
- Completude de dados críticos
- Aderência aos critérios de inclusão/exclusão

---

**Nota | Note:** Todos os dados são sintéticos e criados para fins educacionais. Não devem ser utilizados para inferências clínicas reais. | All data is synthetic and created for educational purposes. Should not be used for real clinical inferences.

**Atualizado | Updated:** Setembro 2025  
**Versão | Version:** 1.0
