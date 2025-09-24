# population_genetics.R
# Funções de Genética de Populações: Fst, HWE, LD, estrutura populacional
# Pacotes: genetics, hierfstat, adegenet, SNPRelate (opcional)

#' Teste de Equilíbrio de Hardy-Weinberg por SNP
#' @param geno vetor genotípico (AA/AB/BB ou 0/1/2) para um SNP
#' @return lista com p-valor do teste e frequências
hwe_test <- function(geno){
  if (!requireNamespace("HardyWeinberg", quietly=TRUE)) stop("Package 'HardyWeinberg' is required")
  g <- geno
  if (is.numeric(g)){
    # transformar 0/1/2 em contagens
    AA <- sum(g==0, na.rm=TRUE); AB <- sum(g==1, na.rm=TRUE); BB <- sum(g==2, na.rm=TRUE)
  } else {
    g <- as.character(g)
    AA <- sum(g %in% c("AA","A/A","A A"), na.rm=TRUE)
    AB <- sum(g %in% c("AB","A/B","A B","BA","B/A","B A"), na.rm=TRUE)
    BB <- sum(g %in% c("BB","B/B","B B"), na.rm=TRUE)
  }
  res <- HardyWeinberg::HWChisq(c(AA, AB, BB))
  p <- res$pval
  n <- AA+AB+BB
  pA <- (2*AA + AB)/(2*n)
  list(p.value=p, freq_A=pA, counts=c(AA=AA,AB=AB,BB=BB))
}

#' Medida de diferenciação Fst (Weir & Cockerham) entre populações
#' @param genos lista de matrizes por população (indiv x SNP, 0/1/2)
#' @return Fst por SNP e global
fst_weir_cockerham <- function(genos){
  if (!requireNamespace("hierfstat", quietly=TRUE)) stop("Package 'hierfstat' is required")
  # Empilhar dados com labels populacionais
  mats <- lapply(seq_along(genos), function(i){
    m <- genos[[i]]; if (!is.matrix(m)) m <- as.matrix(m)
    data.frame(pop=i, m, check.names=FALSE)
  })
  df <- do.call(rbind, mats)
  # hierfstat espera alelos 0/1; transformar genótipo 0/1/2 em pares de alelos
  to_alleles <- function(v){
    # retorna duas colunas de alelos por SNP
    a1 <- ifelse(v==0, 0, ifelse(v==1, 0, 1))
    a2 <- ifelse(v==0, 0, ifelse(v==1, 1, 1))
    cbind(a1, a2)
  }
  pops <- df$pop; df$pop <- NULL
  alle <- do.call(cbind, lapply(df, to_alleles))
  dat <- data.frame(pop=pops, alle)
  fst <- hierfstat::wc(dat)
  list(Fst_global=fst$FST, Theta=fst$Theta)
}

#' Desequilíbrio de ligação (LD) r^2 entre dois SNPs
#' @param g1 vetor 0/1/2 SNP1, @param g2 vetor 0/1/2 SNP2
ld_r2 <- function(g1, g2){
  if (!requireNamespace("genetics", quietly=TRUE)) stop("Package 'genetics' is required")
  tab <- table(g1, g2)
  res <- genetics::LD(tab)
  as.numeric(res$`R^2`)
}
