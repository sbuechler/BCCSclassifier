% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_by_genes.R
\name{filter_mpred_by_genes}
\alias{filter_mpred_by_genes}
\title{Filter the multi-kTSP family by genes in the expression dataset}
\usage{
filter_mpred_by_genes(MP, dat)
}
\arguments{
\item{MP}{A multi-kTSP object, that is one of \code{MPred_erposneg},
\code{MPred_erpos}, \code{MPred_erneg}}

\item{dat}{An object of class matrix with samples as column names and gene symbols
as row names, and entries the corresponding expression levels. The included tibble
\code{bccs_classifier_genes} contains the gene symbols and ENTREZIDs used in the
BCCSclassifier models, along with Affymetrix (hgu133a) and Illumina
(illuminaHuman v3) probes that can be used to represent the associated genes.
Expression values must \bold{not} be scaled or median-centered.}
}
\description{
This function returns a multi-kTSP model containing only those kTSP
models whose genes are all found in the expression data. BCCS prediction
proceeds using only this subfamily.
}
