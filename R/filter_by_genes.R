#' Filter the multi-kTSP family by genes in the expression dataset
#'
#' @description
#' This function returns a multi-kTSP model containing only those kTSP
#' models whose genes are all found in the expression data. BCCS prediction
#' proceeds using only this subfamily.
#'
#' @param MP A multi-kTSP object, that is one of \code{MPred_erposneg},
#'   \code{MPred_erpos}, \code{MPred_erneg}
#'
#' @param dat An object of class matrix with samples as column names and gene symbols
#'   as row names, and entries the corresponding expression levels. The included tibble
#'   \code{bccs_classifier_genes} contains the gene symbols and ENTREZIDs used in the
#'   BCCSclassifier models, along with Affymetrix (hgu133a) and Illumina
#'   (illuminaHuman v3) probes that can be used to represent the associated genes.
#'   Expression values must \bold{not} be scaled or median-centered.
#'
#' @export
filter_mpred_by_genes <- function(MP, dat) {
  MGenes <- purrr::map(
    MP, ~ kpred_genes(.)
  )
  present_ind <- purrr::map_int(
    MGenes,
    ~ as.integer(length(setdiff(., rownames(dat))) == 0)
  )
  MP2 <- MP[present_ind == 1]
  MP2
}

# ----- Helper function ---------


kpred_genes <- function(KP) {
  gene_list <- purrr::map(KP, ~ rownames(.$ktspdat))
  genes <- unique(unlist(gene_list, use.names = F))
  genes
}
