#' Multi-kTSP model objects
#'
#' Classification of breast cancer into BCCS subtypes uses
#'  Multi-kTSP models, introduced in *citation*.
#'
#' @format A Multi-kTSP model is coded in R as what we call an \emph{MPred} object.
#'  An MPred object is a large list, each component of which is a kTSP model.
#'  A kTSP model is itself a list, each component of which is an object of class `ktspair`.
#'  The `ktspair` object tests, for a particular pair of possible subtypes,
#'  whether a sample is more likely to be in one subtype or the other.
#'  The kTSP model includes a component for each possible pair of subtypes.
#'
#' The BCCSclassifier package includes three Multi-kTSP models:
#' \describe{
#'   \item{MPred_erposneg}{For classifying breast cancer samples into
#'   one of BCS1-5}
#'   \item{MPred_erpos}{For classifying ER+ breast cancer samples into
#'   one of PCS1-4}
#'   \item{MPred_erneg}{For classifying ER- breast cancer samples into
#'   one of NCS1-3}
#' }
#'
#' @rdname MPred_erposneg
"MPred_erposneg"

