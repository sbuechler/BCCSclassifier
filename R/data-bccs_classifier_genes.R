#' BCCS Classifier Genes Dataset
#'
#' The dataset \code{bccs_classifier_genes} is a tibble (data.frame) with one record
#' for each of the 430 genes used in the Multi-kTSP models presented here.
#' This data.frame can be used to create the expression matrix needed for
#' BCCS classification from whole-transcriptome datasets created with various
#' technologies.
#'
#' @format A tibble with 430 records and the following columns.
#' \describe{
#' \item{SYMBOL}{<chr> gene symbol}
#' \item{ENTREZID}{<chr> ENTREZID of the gene}
#'  \item{affy_feature}{<chr> preferred Affymetrix hgu133a probe to
#'  represent the gene}
#'  \item{illumina_feature}{<chr> preferred IlluminaHuman v3 probe to
#'  represent the gene}
#'  \item{ensembl_gene}{<chr> Ensembl format gene identifier}
#' }
#'
#' To create the needed expression matrix from, e.g., an Affymetrix
#' expression dataset, the features should be restricted to the probes in the
#' \code{affy_feature} column, row names set to the corresponding \code{SYMBOL}
#' and column names set to the sample identifiers.
#'
#' @rdname bccs_classifier_genes
"bccs_classifier_genes"

