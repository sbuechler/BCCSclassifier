% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BCCSclassifier.R
\docType{package}
\name{BCCSclassifier-package}
\alias{BCCSclassifier}
\alias{BCCSclassifier-package}
\title{Breast Cancer Consensus Subtypes Classification Package}
\description{
A package for computing the Breast Cancer Consensus Subtypes from
whole-transcriptome data using the multi-kTSP method. For details, please see the vignette:
\code{browseVignettes(package = "BCCSclassifier")}
}
\section{Top-level BCCSclassifier functions}{


\code{predict_bccs} predicts the subtype of each sample with respect to one
of three models: all breast cancer samples (erposneg),
estrogen-receptor positive samples (erpos), or estrogen-receptor negative
samples (erneg).

\code{predict_bccs_frequency} reports, for each sample and subtype, the fraction of
kTSP models in the family that predict the sample is in that subtype.
}

\seealso{
Useful links:
\itemize{
  \item Report bugs at \url{https://github.com/sbuechler/BCCSclassifier/issues}
}

}
\author{
\strong{Maintainer}: Steven Buechler \email{steve@nd.edu} (\href{https://orcid.org/0000-0003-3018-9206}{ORCID})

Authors:
\itemize{
  \item Christina Horr \email{chorr@alumni.nd.edu} (\href{https://orcid.org/0000-0001-9684-2840}{ORCID})
}

}
