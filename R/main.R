#' BCCS predictor
#'
#' A function to predict the subtypes of breast cancers with respect to a
#' BCCSclassifier model
#'
#' Learn more in the vignette by executing `vignette("BCCSclassifier")`
#'
#' @param dat An object of class matrix with samples as column names, gene symbols
#'   as row names, and entries the corresponding expression levels. The included tibble
#'   \code{bccs_classifier_genes} contains the gene symbols and ENTREZIDs used in the
#'   BCCSclassifier models, along with Affymetrix (hgu133a) and Illumina
#'   (illuminaHuman v3) probes that can be used to represent the associated genes.
#'   Expression values must \bold{not} be scaled or median-centered.
#'
#' @param model A character string that is one of "erposneg", "erpos", "erneg",
#'   to determine whether BCCSclassifier(ER+/-), BCCSclassifier(ER+), or
#'   BCCSclassifier(ER-) should be used for subtyping.
#'
#' @return The output is a tibble (data.frame) with columns sample, subtype_prediction,
#'   frequency. There is one record for each sample in the column names of \code{dat}.
#'   \code{frequency} is the fraction of kTSP models in the family that predicted the
#'   sample was in the finally chosen subtype.
#'
#' @export
predict_bccs <- function(dat, model) {
  PPI <- make_pairwise_predictions(dat, model)
  SPI <- make_subtype_predictions(PPI)
  SPF <- fequency_of_subtype(SPI)
  final_bccs_pred <- select_final_prediction(SPF)
  final_bccs_pred
}

#' Frequency of BCCS subtype prediction
#'
#' @description
#' A function to predict for each sample and subtype, the fraction of
#' kTSP models in the family that predict the sample is in the subtype.
#'
#' Learn more in `vignette("BCCSclassifier")`
#'
#' @param dat An object of class matrix with samples as column names and gene symbols
#'   as row names, and entries the corresponding expression levels. The included tibble
#'   \code{bccs_classifier_genes} contains the gene symbols and ENTREZIDs used in the
#'   BCCSclassifier models, along with Affymetrix (hgu133a) and Illumina
#'   (illuminaHuman v3) probes that can be used to represent the associated genes.
#'   Expression values must \bold{not} be scaled or median-centered.
#'
#' @param model A character string that is one of "erposneg", "erpos", "erneg",
#'   to determine whether BCCSclassifier(ER+/-), BCCSclassifier(ER+), or
#'   BCCSclassifier(ER-) should be used for subtyping.
#'
#' @return The output is a tibble (data.frame) with columns sample, subtype_prediction,
#'   frequency. There is one record for each sample in the column names of \code{dat} and
#'   each possible subtype. \code{frequency} is the fraction of kTSP models in
#'   the family that predicted the sample was in the corresponding subtype.
#'
#' @export
predict_bccs_frequency <- function(dat, model) {
  PPI <- make_pairwise_predictions(dat, model)
  SPI <- make_subtype_predictions(PPI)
  SPF <- fequency_of_subtype(SPI)
  SPF
}
