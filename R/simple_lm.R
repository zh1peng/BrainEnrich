#' Perform Linear Regression with Multiple Predictors and Covariates
#'
#' This function fits linear models for specified dependent variables using given predictors and covariates.
#' It returns a data frame containing model summaries.
#'
#' @param dependent_df A data frame containing the dependent variables.
#' @param pred_df A data frame containing the predictor variables.
#' @param cov_df A data frame containing the covariate variables.
#' @param stat2return A character string specifying which statistic to return ("statistic", "p.value", or "full"). Default is "full". "statistic" returns only the t-value for permutation purposes, "p.value" returns only the p-value for simulation analysis, and "full" returns all information for the parametric test.
#' @return A data frame containing model summaries. Depending on `stat2return`, the output can include different statistics:
#' \itemize{
#'   \item If `stat2return` is "all", the output includes unstandardized and standardized coefficients, standard errors, t-values, confidence intervals, p-values, adjusted p-values, and significance markers.
#'   \item If `stat2return` is "tval", the output includes only the t-values.
#'   \item If `stat2return` is "tval", the output includes only the t-values as a list.
#'   \item If `stat2return` is "pval", the output includes only the p-values.
#' }
#' @importFrom dplyr filter mutate select ungroup %>% all_of case_when
#' @importFrom tidyr pivot_longer unnest
#' @importFrom purrr map
#' @importFrom broom tidy
#' @importFrom parameters standardize_parameters
#' @importFrom stats p.adjust relevel
#' @importFrom tibble deframe
#' @importFrom rlang .data 
#' @export
simple_lm <- function(dependent_df,
                      pred_df,
                      cov_df,
                      stat2return = c("all", "tval", "pval", "tval_list")) {
  stat2return <- match.arg(stat2return)
  # Check if the specified variables exist in the data frame
  df <- cbind(dependent_df, pred_df, cov_df)
  dependent_vars <- colnames(dependent_df)
  pred_var <- colnames(pred_df)
  cov_vars <- colnames(cov_df)

  if (is.numeric(pred_df[[pred_var]])) {
    # If the column is a continuous variable
    var2extract <- pred_var
  } else if (is.factor(pred_df[[pred_var]])) {
    levels_var <- levels(pred_df[[pred_var]])
    if (length(levels_var) == 2) {
      # Ensure the first level is the reference
      pred_df[[pred_var]] <- relevel(pred_df[[pred_var]], ref = levels_var[1])
      # Create var2extract with the second level
      var2extract <- paste0(pred_var, levels_var[2])
    } else {
      # Stop if the factor has more than two levels
      stop("If you have more than two levels in pred_df, consider analysis directly with brain score data.")
    }
  }

  # Pivot longer, group, nest, and fit models

  if (stat2return == "all") {
    res <- df %>%
      tidyr::pivot_longer(cols = all_of(dependent_vars), names_to = "Dependent_vars", values_to = "Dependent_value") %>%
      dplyr::group_by(.data$Dependent_vars) %>%
      tidyr::nest() %>%
      dplyr::mutate(
        lm_model = purrr::map(.data$data, ~ eval(bquote(lm(.(as.formula(paste("Dependent_value ~", paste(c(pred_var, cov_vars), collapse = "+")))), data = .x)))),
        tidy_model = purrr::map(.data$lm_model, broom::tidy),
        std_coefs = purrr::map(.data$lm_model, ~ parameters::standardize_parameters(.x, method = "refit"))
      ) %>%
      tidyr::unnest(.data$tidy_model) %>%
      dplyr::filter(.data$term == var2extract) %>%
      tidyr::unnest(.data$std_coefs) %>%
      dplyr::filter(.data$Parameter == var2extract) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        p.adj = p.adjust(.data$p.value, method = "fdr")
        # ifsig = dplyr::case_when(
        #   .data$p.adj < 0.001 ~ "***",
        #   .data$p.adj < 0.01 ~ "**",
        #   .data$p.adj < 0.05 ~ "*",
        #   TRUE ~ "n.s."
        # )
      ) %>%
      dplyr::select(.data$Dependent_vars, .data$term, .data$estimate, .data$std.error, .data$statistic, .data$Std_Coefficient, .data$CI_low, .data$CI_high, .data$p.value, .data$p.adj) %>%
      dplyr::rename(
        Predictor = .data$term,
        Unstandardized_Coefficient = .data$estimate,
        Standard_Error = .data$std.error,
        t_Value = .data$statistic,
        Standardized_Coefficient = .data$Std_Coefficient,
        CI_95_Lower = .data$CI_low,
        CI_95_Upper = .data$CI_high,
        p.val = .data$p.value
      )
  } else if (stat2return == "tval") {
    res <- df %>%
      tidyr::pivot_longer(cols = all_of(dependent_vars), names_to = "Dependent_vars", values_to = "Dependent_value") %>%
      dplyr::group_by(.data$Dependent_vars) %>%
      tidyr::nest() %>%
      dplyr::mutate(
        lm_model = purrr::map(.data$data, ~ eval(bquote(lm(.(as.formula(paste("Dependent_value ~", paste(c(pred_var, cov_vars), collapse = "+")))), data = .x)))),
        tidy_model = purrr::map(.data$lm_model, broom::tidy)
      ) %>%
      tidyr::unnest(.data$tidy_model) %>%
      dplyr::filter(.data$term == var2extract) %>%
      dplyr::select(.data$Dependent_vars, .data$term, .data$statistic) %>%
      dplyr::rename(Predictor = .data$term, tval = .data$statistic) %>%
      dplyr::ungroup()
  } else if (stat2return == "tval_list") {
    res <- df %>%
      tidyr::pivot_longer(cols = all_of(dependent_vars), names_to = "Dependent_vars", values_to = "Dependent_value") %>%
      dplyr::group_by(.data$Dependent_vars) %>%
      tidyr::nest() %>%
      dplyr::mutate(
        lm_model = purrr::map(.data$data, ~ eval(bquote(lm(.(as.formula(paste("Dependent_value ~", paste(c(pred_var, cov_vars), collapse = "+")))), data = .x)))),
        tidy_model = purrr::map(.data$lm_model, broom::tidy)
      ) %>%
      tidyr::unnest(.data$tidy_model) %>%
      dplyr::filter(.data$term == var2extract) %>%
      dplyr::select(.data$Dependent_vars, .data$statistic) %>%
      tibble::deframe() %>%
      as.list()
  } else if (stat2return == "pval") {
    res <- df %>%
      tidyr::pivot_longer(cols = all_of(dependent_vars), names_to = "Dependent_vars", values_to = "Dependent_value") %>%
      dplyr::group_by(.data$Dependent_vars) %>%
      tidyr::nest() %>%
      dplyr::mutate(
        lm_model = purrr::map(.data$data, ~ eval(bquote(lm(.(as.formula(paste("Dependent_value ~", paste(c(pred_var, cov_vars), collapse = "+")))), data = .x)))),
        tidy_model = purrr::map(.data$lm_model, broom::tidy)
      ) %>%
      tidyr::unnest(.data$tidy_model) %>%
      dplyr::filter(.data$term == var2extract) %>%
      dplyr::select(.data$Dependent_vars, .data$term, .data$p.value) %>%
      dplyr::rename(Predictor = .data$term, pval = .data$p.value) %>%
      dplyr::ungroup()
  }
  return(res)
}
