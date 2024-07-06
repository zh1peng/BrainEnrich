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
#'                     'meansqr', 'maxmean', 'ks_orig', 'ks_weighted', 'ks_pos_neg_sum',
#'                     'local_fdr', 'sign_test', 'rank_sum', 'custom'.
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
                         "ks_orig", "ks_weighted", "ks_pos_neg_sum",
                         "local_fdr", "sign_test", "rank_sum", "custom"
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
