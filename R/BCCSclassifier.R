#' Breast Cancer Consensus Subtypes Classification Package
#'
#' @description
#' A package for computing the Breast Cancer Consensus Subtypes from
#' whole-transcriptome data using the multi-kTSP method. For details, please see the vignette:
#'  `browseVignettes(package = "BCCSclassifier")`
#'
#' @section Top-level BCCSclassifier functions:
#'
#' \code{predict_bccs} predicts the subtype of each sample with respect to one
#' of three models: all breast cancer samples (erposneg),
#' estrogen-receptor positive samples (erpos), or estrogen-receptor negative
#' samples (erneg).
#'
#' \code{predict_bccs_frequency} reports, for each sample and subtype, the fraction of
#' kTSP models in the family that predict the sample is in that subtype.
#'
#'
#'
#' @importFrom tibble tibble
#'
#'
#' @docType package
#'
"_PACKAGE"
#> NULL
