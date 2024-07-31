#' Find Core Genes Influencing Aggregated Score or LM Coefficients between molecular profile and behavioral data
#'
#' This function performs a Leave-One-Out (LOO) analysis on gene sets to determine core genes
#' that influence the aggregated score. It can utilize parallel processing
#' to enhance computation efficiency and supports two types of analysis: one that considers
#' only gene sets and another that includes predictor and covariate data frames.
#'
#' @param geneList A matrix of genes by subs, each column representing a subject / a group-level result.
#' @param geneSetList A list of gene sets, each containing names of genes.
#' @param pred_df Optional data frame of a predictor. If NULL, it is perfomred for group-level enrichment.
#' @param cov_df Optional data frame of covariates. If NULL, it is perfomred for group-level enrichment.
#' @param aggre_method The aggregation method used to compute the scores.
#' @param n_cores The number of cores to use for parallel processing; defaults to 1.
#'                 Uses all available cores minus one if set to 0.
#' @param threshold_type The method to determine significance ('sd' for standard deviation, 'percentile' for percentile threshold).
#' @param threshold_value Numeric value specifying the threshold level; meaning depends on `threshold_type`.
#' @import pbapply parallel
#' @return A list of core genes for each gene set.
#' @export
find_core_genes <- function(geneList, geneSetList, pred_df = NULL, cov_df = NULL, aggre_method, n_cores = 1, threshold_type = c("sd", "percentile"), threshold_value = 1) {
  if (is.null(pred_df) & is.null(cov_df)) {
    type1_analysis <- TRUE # find core genes contribute to group enrichment results
  } else {
    type1_analysis <- FALSE
  }

  # Validate the threshold_type argument
  threshold_type <- match.arg(threshold_type)

  if (threshold_type == "sd") {
    if (threshold_value < 0 || threshold_value > 3) {
      stop("For 'sd', threshold_value should be between 0 and 3.")
    }
  } else if (threshold_type == "percentile") {
    if (threshold_value < 1 || threshold_value > 99 || threshold_value %% 10 != 0) {
      stop("For 'percentile', threshold_value should be a multiple of 10 and between 1 and 99.")
    }
  }

  # Determine the number of cores to use
  if (n_cores == 0 | n_cores > detectCores() - 1) {
    n_cores <- detectCores() - 1
  }

  # Initialize a cluster of workers
  cl <- makeCluster(n_cores)

  if (type1_analysis) {
    # Export necessary variables to the cluster
    clusterExport(cl, varlist = c("geneList", "aggregate_geneSet", "aggre_method"), envir = environment())
    # Parallelize the processing using pblapply for progress bar
    loo_changes <- pblapply(seq_along(geneSetList), function(i) {
      gs <- geneSetList[[i]]
      full_score <- aggregate_geneSet(geneList, gs, method = aggre_method)
      # Perform Leave-One-Out Analysis for each gene in the gene set
      loo_results <- sapply(seq_along(gs), function(gs_i) {
        modified_geneSet <- setdiff(gs, gs[gs_i])
        loo_score <- aggregate_geneSet(geneList, modified_geneSet, method = aggre_method)
        adiff <- abs(full_score - loo_score)
        return(adiff)
      })
      names(loo_results) <- gs
      loo_results
    }, cl = cl)
  } else if (!type1_analysis & is.data.frame(pred_df) & is.data.frame(cov_df)) {
    # Export necessary variables to the cluster
    clusterExport(cl, varlist = c("geneList", "aggregate_geneSet", "aggre_method", "pred_df", "cov_df", "simple_lm"), envir = environment())

    loo_changes <- pblapply(seq_along(geneSetList), function(i) {
      gs <- geneSetList[[i]]
      gs_name <- names(geneSetList)[i]
      full_score <- aggregate_geneSet(geneList, gs, method = aggre_method)
      dependent_df.full <- data.frame(full_score)
      colnames(dependent_df.full) <- gs_name
      lm_res.full <- simple_lm(dependent_df = dependent_df.full, pred_df = pred_df, cov_df = cov_df, stat2return = "all")

      # Perform Leave-One-Out Analysis for each gene in the gene set
      loo_results <- sapply(seq_along(gs), function(gs_i) {
        modified_geneSet <- setdiff(gs, gs[gs_i])
        loo_score <- aggregate_geneSet(geneList, modified_geneSet, method = aggre_method)
        dependent_df.loo <- data.frame(loo_score)
        colnames(dependent_df.loo) <- gs_name
        lm_res.loo <- simple_lm(dependent_df = dependent_df.loo, pred_df = pred_df, cov_df = cov_df, stat2return = "all")
        adiff <- abs(lm_res.full$Standardized_Coefficient - lm_res.loo$Standardized_Coefficient)
        return(adiff)
      })
      names(loo_results) <- gs
      loo_results
    }, cl = cl)
  }

  # Stop the cluster after processing
  stopCluster(cl)
  core_genes <- lapply(loo_changes, identify_core_genes, threshold_type = threshold_type, threshold_value = threshold_value)
  names(core_genes) <- names(geneSetList)
  return(core_genes)
}


#' Identify core genes based on a specified threshold method
#'
#' @param changes Named vector of changes from LOO analysis.
#' @param threshold_type Character string indicating the method to determine the threshold ("sd" or "percentile").
#' @param threshold_value Numeric value indicating the percentile (if method is "percentile") or the number of standard deviations (if method is "sd").
#' @importFrom stats quantile sd
#' @return Vector of core genes or NA if no core genes are identified.
identify_core_genes <- function(changes, threshold_type = c("sd", "percentile"), threshold_value = 1) {
  threshold_type <- match.arg(threshold_type)
  # Validate threshold_value based on method
  if (threshold_type == "percentile") {
    if (threshold_value < 1 || threshold_value > 99 || threshold_value %% 10 != 0) {
      stop("For 'percentile', threshold_value should be a multiple of 10 and between 1 and 99.")
    }
    threshold <- quantile(changes, probs = threshold_value / 100)
  } else if (threshold_type == "sd") {
    if (threshold_value < 0 || threshold_value > 4) {
      stop("For 'sd', threshold_value should be between 0 and 4.")
    }
    mean_change <- mean(changes)
    sd_change <- sd(changes)
    threshold <- mean_change + threshold_value * sd_change
  }
  core_genes <- names(changes[changes > threshold])
  # Remove core genes with 0 changes
  core_genes <- core_genes[changes[core_genes] != 0]
  # If no core genes left, return NA
  if (length(core_genes) == 0) {
    core_genes <- "NOT_FOUND"
  }
  return(core_genes)
}
