#' Calculate Brain Scores for Gene Sets
#'
#' This function calculates scores for gene sets based on brain data.
#'
#' @param brain_data A data frame of brain data. Region by 1 column.
#' @param gene_data A data frame of gene expression data.
#' @param annoData An environment containing annotation data.
#' @param cor_method A character string specifying the correlation method.
#'                   Default is 'pearson'. Other options include 'spearman', 'pls1c',
#'                   'pls1w', 'custom'.
#' @param aggre_method A character string specifying the aggregation method.
#'                     Default is 'mean'. Other options include 'median', 'meanabs',
#'                     'meansqr', 'maxmean', 'ks_orig', 'ks_weighted', 'ks_pos_neg_sum', 'sign_test', 'rank_sum', 'custom'.
#' @param prefix A character string to be prefixed to the column names of the result. Default is NULL.
#' @param minGSSize An integer specifying the minimum gene set size. Default is 10.
#' @param maxGSSize An integer specifying the maximum gene set size. Default is 200.
#' @return A data frame containing the gene set scores with regions as rows and gene sets as columns.
#' @export
brainscore <- function(brain_data,
                       gene_data,
                       annoData,
                       cor_method = c("pearson", "spearman", "pls1c", "pls1w", "custom"),
                       aggre_method = c(
                         "mean", "median", "meanabs", "meansqr", "maxmean",
                         "ks_orig", "ks_weighted", "ks_pos_neg_sum", "sign_test", "rank_sum", "custom"
                       ),
                       prefix = NULL,
                       minGSSize = 10,
                       maxGSSize = 200,
                       n_cores = 0) {
  # Match arguments to ensure valid inputs
  cor_method <- match.arg(cor_method)
  aggre_method <- match.arg(aggre_method)

  # Input checks
  stopifnot(is.environment(annoData))

  if (!identical(rownames(gene_data), rownames(brain_data))) {
    stop("Rownames of 'gene_data' and 'brain_data' must be identical.")
  }

  # Calculate gene-brain correlations
  message("Calculating gene-brain correlations...")
  geneList <- corr_brain_gene(gene_data = gene_data, brain_data = brain_data, method = cor_method)

  # Generate gene set list from annotation data
  message("Generating gene set list from annotation data...")
  geneSetList <- get_geneSetList(annoData)

  # Filter gene set list
  message("Filtering gene set list...")
  selected.gs <- filter_geneSetList(rownames(geneList), geneSetList, minGSSize = minGSSize, maxGSSize = maxGSSize)

  # Aggregate gene set scores
  message("Aggregating gene set scores...")
  gs.score <- aggregate_geneSetList(geneList, selected.gs, method = aggre_method, n_cores = n_cores)

  # Prepare result data frame
  res <- data.frame(gs.score)
  rownames(res) <- colnames(brain_data)

  # Add prefix to column names if specified
  if (!is.null(prefix) && length(prefix) > 0 && is.character(prefix)) {
    colnames(res) <- paste0(prefix, names(selected.gs))
  } else {
    colnames(res) <- names(selected.gs)
  }

  return(res)
}


#' Perform Linear Regression with Multiple Predictors and Covariates
#'
#' This function fits linear models for specified dependent variables using given predictors and covariates.
#' It returns a data frame containing model summaries.
#'
#' @param df A data frame containing the data.
#' @param dependent_vars A character vector of dependent variable names.
#' @param pred_var A character vector of predictor variable names.
#' @param cov_vars A character vector of covariate variable names.
#' @param var2extract A character string or vector specifying which predictor's coefficient to extract. If NULL, defaults to pred_var. For a factor variable, specify the level to extract.
#' @param stat2return A character string specifying which statistic to return ("statistic", "p.value", or "full"). Default is "full". "statistic" returns only the t-value for permutation purpose, "p.value" returns only the p-value for simulation analysis, and "full" returns all information for the parametric test.
#' @return A data frame containing model summaries.
#' @import dplyr tidyr purrr broom parameters
#' @export
simple_lm <- function(df,
                      dependent_vars,
                      pred_var,
                      cov_vars,
                      var2extract = NULL,
                      stat2return = c("full", "statistic", "p.value")) {
  if (is.null(var2extract)) {
    var2extract <- pred_var
  }

  # Check if the specified variables exist in the data frame
  all_vars <- c(dependent_vars, pred_var, cov_vars)
  missing_vars <- setdiff(all_vars, colnames(df))
  if (length(missing_vars) > 0) {
    stop("The following variables are not found in the data frame: ", paste(missing_vars, collapse = ", "))
  }

  df_model <- df %>%
    dplyr::select(all_of(pred_var), all_of(cov_vars), all_of(dependent_vars)) %>%
    pivot_longer(cols = all_of(dependent_vars), names_to = "Dependent_vars", values_to = "Dependent_value") %>%
    group_by(Dependent_vars) %>%
    nest() %>%
    mutate(
      lm_model = map(data, ~ lm(paste("Dependent_value", "~", paste(c(pred_var, cov_vars), collapse = "+")), data = .x)),
      tidy_model = map(lm_model, tidy)
    )

  # Conditionally calculate std_coefs only if stat2return is "full"
  if (stat2return == "full") {
    df_model <- df_model %>%
      mutate(std_coefs = map(lm_model, ~ standardize_parameters(.x, method = "refit") %>%
        as_tibble() %>%
        filter(Parameter %in% var2extract)))
  }

  df_model <- df_model %>%
    unnest(tidy_model) %>%
    filter(term == var2extract)

  if (stat2return == "full") {
    df_model <- df_model %>%
      unnest(std_coefs) %>%
      filter(Parameter == var2extract) %>%
      mutate(
        p.adj = p.adjust(p.value, method = "fdr"),
        ifsig = case_when(
          p.adj < 0.001 ~ "***",
          p.adj < 0.01 ~ "**",
          p.adj < 0.05 ~ "*",
          TRUE ~ "n.s."
        )
      ) %>%
      dplyr::select(Dependent_vars, term, estimate, std.error, statistic, Std_Coefficient, CI_low, CI_high, p.value, p.adj, ifsig) %>%
      rename(
        Predictor = term,
        Unstandardized_Coefficient = estimate,
        Standard_Error = std.error,
        t_Value = statistic,
        Standardized_Coefficient = Std_Coefficient,
        CI_95_Lower = CI_low,
        CI_95_Upper = CI_high,
        p.val = p.value
      )
  } else if (stat2return == "statistic") {
    df_model <- df_model %>%
      dplyr::select(Dependent_vars, term, statistic) %>%
      rename(Predictor = term, t_Value = statistic)
  } else if (stat2return == "p.value") {
    df_model <- df_model %>%
      dplyr::select(Dependent_vars, term, p.value) %>%
      rename(Predictor = term, p.val = p.value)
  } else {
    stop("Invalid value for stat2return. Choose from 'statistic', 'p.value', or 'full'.")
  }

  return(df_model)
}


brainscore.comp
