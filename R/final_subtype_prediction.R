#' BCCS prediction from frequency
#'
#' @description
#' This function identifies, for each sample, the subtype that is predicted
#' most frequently across the multi-kTSP family, or reports a "tie" if the
#' maximum isn't unique. The result returns the maximum frequency as a column
#' in the output.
#'
#' @param SPF A tibble with columns sample, subtype_prediction, frequency, as
#'   produced by \code{fequency_of_subtype}.
#'
#' @return A tibble with columns sample, subtype_prediction, frequency, with one
#'   record per sample
#'
#' @export
select_final_prediction <- function(SPF) {
  # Filter to subtype(s) with maximal frequency for each sample
  final_prediction <- SPF %>%
    dplyr::group_by(sample) %>% dplyr::mutate(max_frac = max(frequency)) %>%
    dplyr::filter(frequency >= max_frac) %>% dplyr::ungroup()
  final_prediction <- final_prediction %>%
    dplyr::select(sample, subtype_prediction, frequency)
  # Look for ties in the final prediction
  final_prediction <- final_prediction %>% dplyr::group_by(sample) %>%
    dplyr::mutate(tie_cnt = dplyr::n()) %>% dplyr::ungroup()
  final_pred1 <- final_prediction %>% dplyr::filter(tie_cnt == 1)
  final_pred2 <- final_prediction %>%
    dplyr::filter(tie_cnt > 1) # Multiple subtypes with max frequency
  final_pred2 <- final_pred2 %>% dplyr::mutate(subtype_prediction = "tie")
  final_pred2 <- final_pred2 %>% dplyr::distinct()
  final_prediction <- dplyr::bind_rows(final_pred1, final_pred2) %>%
    dplyr::select(-tie_cnt)
  final_prediction
}

