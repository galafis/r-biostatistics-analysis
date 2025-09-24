# gwas_analysis.R
# Estudos de Associação Genômica Ampla (GWAS)
# Pacotes: GenABEL/GENESIS, snpStats, data.table, qqman
# Obs.: GenABEL está arquivado; alternativamente usar SNPRelate/GENESIS quando disponível

#' Executa pipeline básico de GWAS (qualidade, codificação, associação)
#' @param geno matriz/objeto genotípico (snpStats::SnpMatrix) ou path para PLINK bed/bim/fam
#' @param pheno data.frame com id, fenótipo e covariáveis
#' @param id_col nome da coluna de ID para ligar fenótipo e genótipo
#' @param trait nome do fenótipo (quantitativo ou binário 0/1)
#' @param covars vetor de covariáveis (ex.: c("age","sex","PC1","PC2"))
#' @param binary lógico, TRUE para traço binário (logístico)
#' @param maf_min frequência alélica mínima (default 0.01)
#' @param call_rate_min proporção mínima de chamadas por SNP (default 0.95)
#' @param hwe_pmin p-valor mínimo para HWE em controles ou global (default 1e-6)
#' @return lista com resultados: tabela de associação, Manhattan/QQ e estatísticas
#' @examples
#' # library(snpStats)
#' # g <- read.plink("data/genotypes.bed")
#' # res <- gwas_run(g$genotypes, pheno_df, id_col="IID", trait="BMI",
#' #                 covars=c("age","sex","PC1","PC2"), binary=FALSE)
#' # head(res$assoc)
#' gwas_run <- function(geno, pheno, id_col, trait, covars=NULL, binary=FALSE,
#'                      maf_min=0.01, call_rate_min=0.95, hwe_pmin=1e-6){
#'   if (!requireNamespace("snpStats", quietly=TRUE)) stop("Package 'snpStats' is required")
#'   if (!requireNamespace("data.table", quietly=TRUE)) stop("Package 'data.table' is required")
#'   # Importar genótipos se caminho PLINK foi fornecido
#'   if (inherits(geno, "character")){
#'     pl <- snpStats::read.plink(bed= paste0(geno, ".bed"),
#'                                bim= paste0(geno, ".bim"),
#'                                fam= paste0(geno, ".fam"))
#'     G <- pl$genotypes
#'     map <- data.frame(chr=pl$map$chromosome, snp=pl$map$snp.name,
#'                       pos=pl$map$position, stringsAsFactors=FALSE)
#'     fam <- pl$fam
#'     sample_ids <- rownames(G)
#'   } else {
#'     G <- geno
#'     if (is.null(rownames(G))) stop("Genotype matrix must have rownames as IDs")
#'     sample_ids <- rownames(G)
#'     map <- NULL
#'   }
#'   # Armonizar IDs
#'   stopifnot(id_col %in% names(pheno), trait %in% names(pheno))
#'   ph <- data.table::as.data.table(pheno)
#'   ph <- ph[get(id_col) %in% sample_ids]
#'   # Subconjunto de amostras
#'   keep <- match(ph[[id_col]], sample_ids)
#'   G <- G[keep, , drop=FALSE]
#'   rownames(G) <- ph[[id_col]]
#'   # Filtros de qualidade
#'   cr <- 1 - snpStats::col.summary(G)$Call.rate # errado, ajustar abaixo
#'   cs <- snpStats::col.summary(G)
#'   maf <- pmin(cs$MAF, 1-cs$MAF)
#'   callrate <- cs$Call.rate
#'   hwe <- cs$z.HWE # aproximação; usar pval se disponível
#'   pass <- which(maf >= maf_min & callrate >= call_rate_min)
#'   if (length(pass)==0) stop("No SNPs pass QC. Relax thresholds.")
#'   G <- G[, pass]
#'   if (!is.null(map)) map <- map[pass, , drop=FALSE]
#'   # Modelo
#'   df <- data.frame(y = ph[[trait]])
#'   if (!is.null(covars)){
#'     stopifnot(all(covars %in% names(ph)))
#'     df <- cbind(df, ph[, covars, with=FALSE])
#'   }
#'   # Função de teste por-SNP
#'   test_snp <- function(g){
#'     # g é vetor genotípico (0/1/2 ou NA)
#'     if (binary){
#'       fit <- try(stats::glm(df$y ~ g + ., data=df, family=stats::binomial()), silent=TRUE)
#'     } else {
#'       fit <- try(stats::lm(df$y ~ g + ., data=df), silent=TRUE)
#'     }
#'     if (inherits(fit, "try-error")) return(c(NA, NA, NA))
#'     co <- summary(fit)$coefficients
#'     if (!("g" %in% rownames(co))) return(c(NA, NA, NA))
#'     beta <- co["g", 1]
#'     se <- co["g", 2]
#'     p <- co["g", 4]
#'     c(beta, se, p)
#'   }
#'   # Aplicar por coluna (SNP)
#'   mat <- as.matrix(G)
#'   res <- apply(mat, 2, test_snp)
#'   res <- t(res)
#'   colnames(res) <- c("beta","se","p")
#'   res <- as.data.frame(res)
#'   res$snp <- colnames(G)
#'   if (!is.null(map)){
#'     res$chr <- map$chr; res$pos <- map$pos
#'   }
#'   # Ajuste múltiplo e lambda GC
#'   res$padj_bh <- stats::p.adjust(res$p, method="BH")
#'   chi2 <- stats::qchisq(1-res$p, df=1)
#'   lambda_gc <- stats::median(chi2, na.rm=TRUE)/0.456
#'   # Gráficos opcionais
#'   manh <- NULL; qq <- NULL
#'   if (requireNamespace("qqman", quietly=TRUE) && !is.null(map)){
#'     manh <- try(qqman::manhattan(res, chr="chr", bp="pos", snp="snp", p="p",
#'                                  main="GWAS Manhattan", genomewideline=-log10(5e-8), suggestiveline=TRUE), silent=TRUE)
#'     qq <- try(qqman::qq(res$p), silent=TRUE)
#'   }
#'   list(assoc=res, lambda_gc=lambda_gc, manhattan_plot=manh, qq_plot=qq)
#' }
