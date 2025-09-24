# linkage_analysis.R
# Análise de Ligação (familial linkage analysis) para traços mendelianos
# Pacotes: genetics, paramlink, kinship2

#' Calcula LOD score para marcador(s) em pedigree(s)
#' @param ped objeto pedigree (kinship2::pedigree) ou data.frame com id, father, mother, sex, affect
#' @param markers data.frame/matriz de genótipos por indivíduo (marcadores codificados)
#' @param theta recombinação presumida (0-0.5)
#' @param model modelo paramétrico (penetrâncias) opcional
#' @return lista com LOD estimado por marcador
#' @examples
#' # ped <- kinship2::pedigree(id=1:4, dadid=c(0,0,1,1), momid=c(0,0,2,2), sex=c(1,2,1,2), affected=c(0,0,1,0))
#' # lod <- linkage_lod(ped, markers_df, theta=0.01)
linkage_lod <- function(ped, markers, theta=0.01, model=NULL){
  if (!requireNamespace("paramlink", quietly=TRUE)) stop("Package 'paramlink' is required")
  if (!requireNamespace("kinship2", quietly=TRUE)) stop("Package 'kinship2' is required")
  # Converter pedigree se necessário
  if (!inherits(ped, "pedigree")){
    stop("Provide kinship2::pedigree object for robust analysis")
  }
  # Construir objeto linkdat
  ld <- paramlink::linkdat(ped)
  # Adicionar marcadores
  if (is.data.frame(markers)) markers <- as.matrix(markers)
  if (is.null(colnames(markers))) colnames(markers) <- paste0("M", seq_len(ncol(markers)))
  for (j in seq_len(ncol(markers))){
    geno <- markers[, j]
    mk <- paramlink::marker(ld, alleles=sort(unique(na.omit(geno))), genotypes=geno)
    ld <- paramlink::setMarkers(ld, mk)
  }
  # Modelo paramétrico se fornecido
  if (!is.null(model)) ld <- paramlink::setModel(ld, model=model)
  # Calcular LOD por marcador
  lods <- sapply(seq_len(paramlink::nMarkers(ld)), function(i){
    paramlink::lod(ld, theta=theta, markers=i)
  })
  data.frame(marker=colnames(markers), LOD=as.numeric(lods))
}
