
# BCCSclassifier

## Overview
BCCSclassifier identifies the Breast Cancer Consensus Subtypes (BCCS) of
breast cancer samples based on the expression values of genes as
published in Horr, C., Buechler, S. Breast Cancer Consensus Subtypes: A system for subtyping breast cancer tumors based on gene expression. *npj Breast Cancer*. (2021). The package is for research use only, and
expects that expression values were obtained from fresh-frozen primary
tumor samples.

## Installation
After installing the devtools package, BCCSclassifier can be installed
by

```r 
devtools::install_github("sbuechler/BCCSclassifier", build_vignettes = TRUE)
```

## Basic usage
The BCCS classifications of a dataset of breast cancer samples can be
computed with the `predict_bccs` function.  `predict_bccs` takes as
arguments:

**dat** - a matrix of gene expression values for a dataset of breast
cancer samples, with sample identifiers as column names and gene symbols
as row names, and

**model** - one of the strings "erpos", "erneg", "erposneg",  depending
on the subtyping model with which the samples should be analyzed.

`predict_bccs` can be used for a large dataset of samples or a single
sample. Other requirements on the expression matrix are described in
`vignette("BCCSclassifier")`.

`predict_bccs` returns a data.frame (more precisely, a tibble) with
columns sample, subtype_prediction (the predicted classification of the
sample), and frequency (the fraction of kTSP models in the family that
classified the sample into this subtype).

Much more on the functions underlying `predict_bccs` can be found in
`vignette("BCCSclassifier")` and documentation for these component
functions.
