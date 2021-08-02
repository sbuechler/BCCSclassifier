#' Predict subtype assignments, by index
#'
#' Predicts, for each sample and kTSP model, the best subtype using the
#' pairwise test results in the given kTSP model
#'
#'
#' @description
#' This function selects the subtype that wins the most pairwise competitions, for each
#'   kTSP model in the family
#'
#' @param PPI A tibble with columns sample, pair, pairwise_prediction, index,
#'   as produced by the \code{make_pairwise_predictions} function.
#'
#' @return A tibble with columns sample, subtype_prediction, index, that we
#'   consider an SPI object, short for \emph{subtype prediction - indexed}
#'
#' @export
make_subtype_predictions <- function(PPI) {
  predict_subtp <- PPI %>% dplyr::group_by(sample, index) %>% tidyr::nest() %>%
    dplyr::mutate(
      subtype_prediction = purrr::map_chr(
        data,
        ~ get_max(.)
      )
    ) %>% dplyr::select(-data)
  predict_subtp <- predict_subtp %>% dplyr::ungroup() %>%
    dplyr::select(sample, subtype_prediction, index)
  predict_subtp
}

# -----  inner function  ------


get_max <- function(df) {
  tt <- table(df$pairwise_prediction)
  tt1 <- tt[tt == max(tt)]
  if(length(tt1) > 1) {out <- "tie"}
  else {out <- names(tt1)}
  out
}

