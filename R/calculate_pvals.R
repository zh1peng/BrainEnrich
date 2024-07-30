#' Calculate P-Values
#'
#' This function calculates p-values based on the provided true statistics and null statistics lists.
#' It supports two methods for p-value calculation: 'standard' and 'split_pos_neg'.
#'
#' @param statList.true A named list of true statistics.
#' @param statList.null A named list of null statistics corresponding to the true statistics.
#' @param method The method to be used for p-value calculation. Either 'standard' or 'split_pos_neg'. Default is 'standard'.
#' @return A list of calculated p-values.
#' @importFrom purrr map2
#' @export
calculate_pvals <- function(statList.true, statList.null, method = c("standard", "split_pos_neg")) {
  method <- match.arg(method)

  if (!identical(names(statList.true), names(statList.null))) {
    stop("statList.true and statList.null are not matched")
  }

  pvals <- map2(statList.true, statList.null, function(true_stat, null_stat) {
    if (method == "standard") {
      (sum(abs(null_stat) >= abs(true_stat)) + 1) / (length(null_stat) + 1)
    } else if (method == "split_pos_neg") {
      if (true_stat >= 0) {
        (sum(null_stat >= true_stat) + 1) / (sum(null_stat >= 0) + 1)
      } else {
        (sum(null_stat <= true_stat) + 1) / (sum(null_stat < 0) + 1)
      }
    }
  })

  return(unlist(pvals))
}
