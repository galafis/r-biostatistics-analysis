# Data Dictionary | Dicionário de Dados

## 🇧🇷 Português | 🇺🇸 English

Este dicionário de dados descreve todos os datasets e variáveis utilizados na plataforma de análises bioestatísticas. | This data dictionary describes all datasets and variables used in the biostatistical analysis platform.

## 📁 Estrutura dos Dados | Data Structure

### 📊 Datasets de Ensaios Clínicos | Clinical Trials Datasets

#### `clinical_trials/rct_data.csv`
Dataset exemplo para ensaios clínicos randomizados.

| Variável | Tipo | Descrição | English Description |
|----------|------|-----------|--------------------|
| `subject_id` | integer | Identificador único do participante | Unique participant identifier |
| `treatment` | factor | Grupo de tratamento (Placebo/Ativo) | Treatment group (Placebo/Active) |
| `age` | numeric | Idade em anos | Age in years |
| `sex` | factor | Sexo (M/F) | Sex (M/F) |
| `baseline_score` | numeric | Escore basal (0-100) | Baseline score (0-100) |
| `endpoint_score` | numeric | Escore final (0-100) | Final endpoint score (0-100) |
| `adverse_events` | integer | Número de eventos adversos | Number of adverse events |
| `dropout` | logical | Descontinuação do estudo | Study dropout |
| `site` | factor | Centro de pesquisa | Research site |

#### `clinical_trials/bioequivalence_data.csv`
Dataset para estudos de bioequivalência.

| Variável | Tipo | Descrição | English Description |
|----------|------|-----------|--------------------|
| `subject` | integer | ID do participante | Participant ID |
| `sequence` | factor | Sequência de tratamento (TR/RT) | Treatment sequence (TR/RT) |
| `period` | integer | Período do estudo (1/2) | Study period (1/2) |
| `treatment` | factor | Tratamento (Test/Reference) | Treatment (Test/Reference) |
| `cmax` | numeric | Concentração máxima (ng/mL) | Maximum concentration (ng/mL) |
| `auc_last` | numeric | AUC até última medição | AUC to last measurement |
| `auc_inf` | numeric | AUC até infinito | AUC to infinity |
| `tmax` | numeric | Tempo para Cmax (h) | Time to Cmax (h) |

### ⏳ Datasets de Análise de Sobrevivência | Survival Analysis Datasets

#### `survival/oncology_survival.csv`
Dataset exemplo para análise de sobrevivência em oncologia.

| Variável | Tipo | Descrição | English Description |
|----------|------|-----------|--------------------|
| `patient_id` | integer | ID do paciente | Patient ID |
| `age` | numeric | Idade ao diagnóstico | Age at diagnosis |
| `sex` | factor | Sexo (M/F) | Sex (M/F) |
| `stage` | factor | Estágio tumoral (I-IV) | Tumor stage (I-IV) |
| `histology` | factor | Tipo histológico | Histological type |
| `treatment` | factor | Tratamento recebido | Treatment received |
| `time_months` | numeric | Tempo de seguimento (meses) | Follow-up time (months) |
| `event` | integer | Evento observado (0=censura, 1=morte) | Observed event (0=censored, 1=death) |
| `progression_time` | numeric | Tempo até progressão | Time to progression |
| `progression` | integer | Progressão (0/1) | Progression (0/1) |

#### `survival/cardiovascular_data.csv`
Dataset para estudos cardiovasculares.

| Variável | Tipo | Descrição | English Description |
|----------|------|-----------|--------------------|
| `subject_id` | integer | ID do participante | Subject ID |
| `age` | numeric | Idade (anos) | Age (years) |
| `sex` | factor | Sexo (M/F) | Sex (M/F) |
| `diabetes` | logical | Presença de diabetes | Diabetes presence |
| `hypertension` | logical | Presença de hipertensão | Hypertension presence |
| `smoking` | factor | Status tabágico | Smoking status |
| `cholesterol` | numeric | Colesterol total (mg/dL) | Total cholesterol (mg/dL) |
| `treatment` | factor | Grupo de tratamento | Treatment group |
| `mace_time` | numeric | Tempo até MACE (meses) | Time to MACE (months) |
| `mace_event` | integer | Evento MACE (0/1) | MACE event (0/1) |

### 📈 Datasets Epidemiológicos | Epidemiological Datasets

#### `epidemiological/cohort_study.csv`
Dataset exemplo para estudo de coorte.

| Variável | Tipo | Descrição | English Description |
|----------|------|-----------|--------------------|
| `participant_id` | integer | ID do participante | Participant ID |
| `baseline_age` | numeric | Idade na entrada | Age at entry |
| `sex` | factor | Sexo (M/F) | Sex (M/F) |
| `exposure` | factor | Status de exposição | Exposure status |
| `bmi` | numeric | Índice de massa corporal | Body mass index |
| `education` | factor | Nível educacional | Education level |
| `income` | factor | Faixa de renda | Income bracket |
| `follow_up_years` | numeric | Anos de seguimento | Years of follow-up |
| `outcome` | integer | Desfecho de interesse (0/1) | Outcome of interest (0/1) |
| `confounders` | list | Variáveis confundidoras | Confounding variables |

#### `epidemiological/case_control.csv`
Dataset para estudo caso-controle.

| Variável | Tipo | Descrição | English Description |
|----------|------|-----------|--------------------|
| `id` | integer | ID do participante | Participant ID |
| `case_control` | factor | Status caso/controle | Case/control status |
| `age` | numeric | Idade (anos) | Age (years) |
| `sex` | factor | Sexo (M/F) | Sex (M/F) |
| `exposure_main` | factor | Exposição principal | Main exposure |
| `exposure_duration` | numeric | Duração da exposição | Exposure duration |
| `matching_var1` | factor | Variável de pareamento 1 | Matching variable 1 |
| `matching_var2` | factor | Variável de pareamento 2 | Matching variable 2 |
| `covariate1` | numeric | Covariável 1 | Covariate 1 |
| `covariate2` | factor | Covariável 2 | Covariate 2 |

## 🔬 Metadados dos Datasets | Dataset Metadata

### Informações Gerais | General Information

- **Formato**: CSV (UTF-8)
- **Separador**: Vírgula (,)
- **Codificação de valores ausentes**: NA
- **Controle de qualidade**: Validação automática implementada

### Convenções de Nomenclatura | Naming Conventions

- **IDs**: Sempre em formato numérico sequencial
- **Fatores**: Níveis em ordem alfabética ou clinicamente relevante
- **Tempos**: Sempre em unidades consistentes (meses/anos)
- **Eventos**: Codificação binária (0/1) onde 1 = evento ocorrido

### Padrões de Codificação | Coding Standards

#### Variáveis Categóricas | Categorical Variables
- `sex`: "M" = Masculino/Male, "F" = Feminino/Female
- `treatment`: "Placebo", "Active", "Control", "Intervention"
- `stage`: "I", "II", "III", "IV" (quando aplicável)

#### Variáveis Temporais | Temporal Variables
- Tempos de seguimento em meses (precisão: 0.1)
- Datas em formato ISO (YYYY-MM-DD)
- Eventos definidos consistentemente através dos estudos

## 📋 Validação e Controle de Qualidade | Validation & Quality Control

### Verificações Implementadas | Implemented Checks
1. **Valores ausentes**: Identificação e tratamento
2. **Outliers**: Detecção estatística automática
3. **Consistência**: Validação cruzada entre variáveis
4. **Completude**: Relatório de missingness por variável

### Limitações dos Dados | Data Limitations
- Datasets são exemplos sintéticos para fins didáticos
- Baseados em distribuições realistas da literatura
- Não devem ser usados para inferências clínicas reais

## 📚 Referências | References

1. ICH E9 Statistical Principles for Clinical Trials
2. FDA Guidance for Industry: Good Clinical Practice
3. EMA Guideline on Clinical Trials in Small Populations
4. STROBE Statement for Observational Studies

---

**Última atualização**: Setembro 2025  
**Versão**: 1.0  
**Responsável**: Gabriel Demetrios Lafis  
**Contato**: gabriel.lafis@example.com
