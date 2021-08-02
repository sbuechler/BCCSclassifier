#' Sample Expression Dataset
#'
#' The dataset \code{exp_data_sample} is a matrix of expression data that can be used
#' in test runs of the package functions.
#'
#' @format An object of class "matrix" with 430 rows and 200 columns.
#' \describe{
#' \item{colnames}{ are codes for samples}
#' \item{rownames}{ are gene symbols}
#'  \item{entries}{ expression values}
#' }
#'
#' From a microarray dataset we selected a representative probe of each gene
#' using \code{bccs_classifier_genes}, and randomly sampled 200 samples.
#'
#' @rdname exp_data_sample
"exp_data_sample"

