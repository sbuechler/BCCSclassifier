#' Compute the frequency of subtype prediction
#'
#' @description
#' This function computes, for each sample and each possible subtype, the fraction of
#' kTSP models in the family that predict this subtype assignment
#'
#' @param SPI A tibble with columns sample, subtype_prediction, index, as produced
#'   by \code{make_subtype_predictions}
#'
#' @return A tibble with columns sample, subtype_prediction, frequency. We call
#'   this object SPF, for \emph{subtype prediction frequency}
#'
#' @export
fequency_of_subtype <- function(SPI) {
  SPI <- SPI %>% dplyr::ungroup()
  numb_family <- length(unique(SPI$index)) # size of Multi-kTSP family
  # Compute fraction of sample-subtype pairs
  bccs_subtype_frac1 <- SPI %>% dplyr::group_by(sample, subtype_prediction) %>%
    dplyr::summarise(frequency = dplyr::n() / numb_family, .groups = "drop")
  # pad with rows with 0 frequency for completeness
  subtypes <- bccs_subtype_frac1  %>% dplyr::select(subtype_prediction) %>%
    dplyr::distinct()
  samps <- bccs_subtype_frac1 %>% dplyr::select(sample) %>% dplyr::distinct()
  subs_samps <- dplyr::full_join(samps, subtypes, by = character()) # all possible combinations
  bccs_subtype_frac2 <- bccs_subtype_frac1 %>%
    dplyr::right_join(subs_samps, by = c("sample", "subtype_prediction"))
  # 0 frequency sample-subtype pairs got NA for frequency by right join
  bccs_subtype_frac2$frequency[is.na(bccs_subtype_frac2$frequency)] <- 0 # convert
  bccs_subtype_frac2
}
