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
#' @param verbose A logical specifying whether to print messages. Default is TRUE.
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
                       matchcoexp_max_iter = 1000000,
                       verbose = TRUE) {
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
  if (verbose) {
    message("Calculating gene-brain correlations...")
  }
  geneList <- corr_brain_gene(gene_data = gene_data, brain_data = brain_data, method = cor_method)

  # Generate gene set list from annotation data
  if (verbose) {
    message("Generating gene set list from annotation data...")
  }
  geneSetList <- get_geneSetList(annoData)

  # Filter gene set list
  if (verbose) {
    message("Filtering gene set list...")
  }
  selected.gs <- filter_geneSetList(rownames(geneList), geneSetList, minGSSize = minGSSize, maxGSSize = maxGSSize)
  if (verbose) {
    message("Number of gene sets left: ", length(selected.gs))
  }
  if (null_model == "none") {
    if (verbose) {
      message("Aggregating gene set scores...")
    }
    gs.score <- aggregate_geneSetList(geneList, selected.gs, method = aggre_method, n_cores = n_cores)
  } else if (null_model == "spin_brain") {
    if (verbose) {
      message("Generating null brain data with spin_brain model...")
    }
    if (is.null(perm_id)) {
      perm_id <- rotate_parcellation(coord.l = coord.l, coord.r = coord.r, nrot = n_perm, seed = seed)
    }
    if (verbose) {
      message("Aggregating gene set scores in spin_brain mode...")
    }
    progress_interval <- max(1, round(n_perm / 10))
    gs.score <- lapply(1:n_perm, function(idx) {
      if (idx %% progress_interval == 0) {
        if (verbose) {
          message(paste("Processing permutation", idx, "of", n_perm, "..."))
        }
      }
      null_brain_data <- brain_data[perm_id[, idx], , drop = FALSE]
      rownames(null_brain_data) <- rownames(brain_data)
      geneList.null <- corr_brain_gene(gene_data = gene_data, brain_data = null_brain_data, method = cor_method)
      gs_score.null <- aggregate_geneSetList(geneList.null, selected.gs, method = aggre_method, n_cores = n_cores)
      return(gs_score.null)
    })
  } else if (null_model == "resample_gene") {
    if (verbose) {
      message("Aggregating gene set scores in resample_gene mode...")
    }
    progress_interval <- max(1, round(n_perm / 10))
    gs.score <- lapply(1:n_perm, function(idx) {
      if (idx %% progress_interval == 0) {
        if (verbose) {
          message(paste("Processing permutation", idx, "of", n_perm, "..."))
        }
      }
      geneList.null <- geneList[sample(1:nrow(geneList), size = nrow(geneList), replace = FALSE), ]
      rownames(geneList.null) <- rownames(geneList)
      gs_score.null <- aggregate_geneSetList(geneList.null, selected.gs, method = aggre_method, n_cores = n_cores)
      return(gs_score.null)
    })
  } else if (null_model == "coexp_matched") {
    if (verbose) {
      message("Generating null gene list with coexp_matched model...")
    }
    sampled_gs <- resample_geneSetList_matching_coexp(gene_data, selected.gs, tol = matchcoexp_tol, max_iter = matchcoexp_max_iter, n_perm = n_perm, n_cores = n_cores)

    if (is.null(sampled_gs)) {
      stop("NULL sampled_gs returned.")
    }

    if (!identical(names(geneSetList), names(sampled_gs))) {
      stop("geneSetList and sampled_geneSetList are not matched.")
    }
    if (verbose) {
      message("Aggregating gene set scores in coexp_matched mode...")
    }
    progress_interval <- max(1, round(n_perm / 10))
    gs.score <- lapply(1:n_perm, function(idx) {
      if (idx %% progress_interval == 0) {
        if (verbose) {
          message(paste("Processing permutation", idx, "of", n_perm, "..."))
        }
      }
      sampled_gs_iter <- lapply(sampled_gs, function(x) x[[idx]])
      gs_score.null <- aggregate_geneSetList(geneList, sampled_gs_iter, method = aggre_method, n_cores = n_cores)
      return(gs_score.null)
    })
  }

  # Add name to the list
  if (null_model != "none") {
    names(gs.score) <- paste0("null_", 1:n_perm)
  }

  # Add attributes to geneList
  attr(gs.score, "type") <- null_model
  attr(gs.score, "cor_method") <- cor_method
  attr(gs.score, "aggre_method") <- aggre_method
  return(gs.score)
}
