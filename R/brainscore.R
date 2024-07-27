#' Calculate Brain Scores for Gene Sets
#'
#' This function calculates scores for gene sets based on brain data. It includes options for different null models.
#'
#' @param brain_data A data frame of brain data. Region by 1 column.
#' @param gene_data A data frame of gene expression data.
#' @param annoData An environment containing annotation data.
#' @param cor_method A character string specifying the correlation method. Default is 'pearson'. Other options include 'spearman', 'pls1c', 'pls1w', 'custom'.
#' @param aggre_method A character string specifying the aggregation method. Default is 'mean'. Other options include 'median', 'meanabs', 'meansqr', 'maxmean', 'ks_orig', 'ks_weighted', 'ks_pos_neg_sum', 'sign_test', 'rank_sum', 'custom'.
#' @param null_model A character string specifying the null model. Default is 'none'. Other options include 'spin_brain', 'resample_gene', 'coexp_matched'.
#' @param minGSSize An integer specifying the minimum gene set size. Default is 10.
#' @param maxGSSize An integer specifying the maximum gene set size. Default is 200.
#' @param n_cores An integer specifying the number of cores to use for parallel processing. Default is 0 (no parallel processing).
#' @param n_perm An integer specifying the number of permutations for null models. Default is 5000.
#' @param perm_id A matrix of permutation indices for 'spin_brain' null model. Default is NULL. Either perm_id or any of coord.l or coord.r must be provided if choosing spin_brain mode.
#' @param coord.l A matrix of coordinates for the left hemisphere for 'spin_brain' null model. Default is NULL.
#' @param coord.r A matrix of coordinates for the right hemisphere for 'spin_brain' null model. Default is NULL.
#' @param seed An integer specifying the seed for reproducibility of spinning brain. Default is NULL.
#' @param matchcoexp_tol A numeric value specifying the tolerance for matching co-expression in 'coexp_matched' null model. Default is 0.05.
#' @param matchcoexp_max_iter An integer specifying the maximum iterations for matching co-expression in 'coexp_matched' null model. Default is 1000000.
#' @return A data frame containing the gene set scores with regions as rows and gene sets as columns.
#' @export
brainscore <- function(brain_data,
                       gene_data,
                       annoData,
                       cor_method = c("pearson", "spearman", "pls1c", "pls1w", "custom"),
                       aggre_method = c("mean", "median", "meanabs", "meansqr", "maxmean", "ks_orig", "ks_weighted", "ks_pos_neg_sum", "sign_test", "rank_sum", "custom"),
                       null_model = c("none", "spin_brain", "resample_gene", "coexp_matched"),
                       minGSSize = 10,
                       maxGSSize = 200,
                       n_cores = 0,
                       n_perm = 5000,
                       perm_id = NULL,
                       coord.l = NULL,
                       coord.r = NULL,
                       seed = NULL,
                       matchcoexp_tol = 0.05,
                       matchcoexp_max_iter = 1000000) {
  
  # Check inputs
  stopifnot(is.environment(annoData))
  stopifnot(identical(rownames(gene_data), rownames(brain_data)))
  
  cor_method <- match.arg(cor_method)
  aggre_method <- match.arg(aggre_method)
  null_model <- match.arg(null_model)
  
  if (null_model == "spin_brain" && is.null(perm_id) && is.null(coord.l) && is.null(coord.r)) {
    stop("For null_model 'spin_brain', 'perm_id' or at least one of 'coord.l' or 'coord.r' must be provided.")
  }
  
  if (null_model == "spin_brain" && !is.null(perm_id) && n_perm > ncol(perm_id)) {
    stop("The number of permutations must be less than or equal to the number of columns in 'perm_id'.")
  } else if (null_model == "spin_brain" && !is.null(perm_id) && n_perm < ncol(perm_id)) {
    perm_id <- perm_id[, 1:n_perm]
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
  
  if (null_model == 'none') {
    message("Aggregating gene set scores...")
    gs.score <- aggregate_geneSetList(geneList, selected.gs, method = aggre_method, n_cores = n_cores)
    
  } else if (null_model == "spin_brain") {
    message("Generating null brain data with spin_brain model...")
    if (is.null(perm_id)) {
      perm_id <- rotate_parcellation(coord.l = coord.l, coord.r = coord.r, nrot = n_perm, seed=seed)
    }
    message("Aggregating gene set scores in spin_brain mode...")
    progress_interval <- max(1, round(n_perm / 10))
    gs.score <- lapply(1:n_perm, function(idx) {
      if (idx %% progress_interval == 0) {
        message(paste("Processing permutation", idx, "of", n_perm, "..."))
      }
      null_brain_data <- brain_data[perm_id[, idx], , drop = FALSE]
      rownames(null_brain_data) <- rownames(brain_data)
      geneList.null <- corr_brain_gene(gene_data = gene_data, brain_data = null_brain_data, method = cor_method)
      gs_score.null <- aggregate_geneSetList(geneList.null, selected.gs, method = aggre_method, n_cores = n_cores)
      return(gs_score.null)
    })
    
  } else if (null_model == "resample_gene") {
    message("Aggregating gene set scores in resample_gene mode...")
    progress_interval <- max(1, round(n_perm / 10))
    gs.score <- lapply(1:n_perm, function(idx) {
      if (idx %% progress_interval == 0) {
        message(paste("Processing permutation", idx, "of", n_perm, "..."))
      }
      geneList.null <- geneList[sample(1:nrow(geneList), size = nrow(geneList), replace = FALSE), ]
      rownames(geneList.null) <- rownames(geneList)
      gs_score.null <- aggregate_geneSetList(geneList.null, selected.gs, method = aggre_method, n_cores = n_cores)
      return(gs_score.null)
    })
    
  } else if (null_model == "coexp_matched") {
    message("Generating null gene list with coexp_matched model...")
    sampled_gs <- resample_geneSetList_matching_coexp(gene_data, selected.gs, tol = matchcoexp_tol, max_iter = matchcoexp_max_iter, n_perm = n_perm, n_cores = n_cores)
    
    if (is.null(sampled_gs)) {
      stop("NULL sampled_gs returned.")
    }
    
    if (!identical(names(geneSetList), names(sampled_gs))) {
      stop("geneSetList and sampled_geneSetList are not matched.")
    }
    message("Aggregating gene set scores in coexp_matched mode...")
    progress_interval <- max(1, round(n_perm / 10))
    gs.score <- lapply(1:n_perm, function(idx) {
      if (idx %% progress_interval == 0) {
        message(paste("Processing permutation", idx, "of", n_perm, "..."))
      }
      sampled_gs_iter <- lapply(sampled_gs, function(x) x[[idx]])
      gs_score.null <- aggregate_geneSetList(geneList, sampled_gs_iter, method = aggre_method, n_cores = n_cores)
      return(gs_score.null)
    })
  }
  
  # Add name to the list
  if (null_model != 'none') {
    names(gs.score) <- paste0('null_', 1:n_perm)
  }
  
  return(gs.score)
}


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
#' @import dplyr tidyr purrr broom parameters tibble
#' @export
simple_lm <- function(dependent_df,
                      pred_df,
                      cov_df,
                      stat2return = c("all", "tval", "pval",'tval_list')) {
  
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
  res <- df %>%
    tidyr::pivot_longer(cols = all_of(dependent_vars), names_to = "Dependent_vars", values_to = "Dependent_value") %>%
    dplyr::group_by(Dependent_vars) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      lm_model = map(data, ~ eval(bquote(lm(.(as.formula(paste("Dependent_value ~", paste(c(pred_var, cov_vars), collapse = "+")))), data = .x)))),
      tidy_model = map(lm_model, tidy)
    ) %>%
    {if (stat2return == "all") {
         dplyr::mutate(., std_coefs = map(lm_model, ~ standardize_parameters(.x, method = "refit")))
      } else {
        .
      }
    } %>%
    tidyr::unnest(tidy_model) %>%
    dplyr::filter(term == var2extract) %>%
    {if (stat2return == "all") {
          unnest(., std_coefs) %>%
          dplyr::filter(Parameter == var2extract) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(
            p.adj = p.adjust(p.value, method = "fdr"),
            ifsig = case_when(
              p.adj < 0.001 ~ "***",
              p.adj < 0.01 ~ "**",
              p.adj < 0.05 ~ "*",
              TRUE ~ "n.s."
            )
          ) %>%
          dplyr::select(Dependent_vars, term, estimate, std.error, statistic, Std_Coefficient, CI_low, CI_high, p.value, p.adj, ifsig) %>%
          dplyr::rename(
            Predictor = term,
            Unstandardized_Coefficient = estimate,
            Standard_Error = std.error,
            t_Value = statistic,
            Standardized_Coefficient = Std_Coefficient,
            CI_95_Lower = CI_low,
            CI_95_Upper = CI_high,
            p.val = p.value
          )
      } else if (stat2return == "tval") {
          dplyr::select(., Dependent_vars, term, statistic) %>%
          dplyr::rename(Predictor = term, tval = statistic) %>%
          dplyr::ungroup()
      } else if (stat2return == "tval_list") {
          dplyr::select(., Dependent_vars, statistic) %>%
          tibble::deframe() %>%
          as.list()
      }
       else if (stat2return == "pval") {
          dplyr::select(., Dependent_vars, term, p.value) %>%
          dplyr::rename(Predictor = term, pval = p.value) %>%
          dplyr::ungroup()
      } 
    }
  return(res)
}