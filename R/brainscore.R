#' Calculate Brain Scores for Gene Sets
#'
#' This function calculates scores for gene sets based on brain data. It supports different null models.
#' If a null model is specified, the function calculates null scores based on the chosen model.
#' If the null model is 'none', the function calculates the raw (empirical) scores.
#'
#' @param brain_data A data frame of brain data with regions as rows and subjects as columns. The row names (i.e., region names) must match those in gene_data.
#' @param gene_data A data frame of gene expression data with regions as rows and genes as columns. The row names (i.e., region names) must match those in brain_data.
#' @param annoData An environment containing annotation data. See `get_annoData` for more details.
#' @param cor_method A character string specifying the correlation method.
#'                   Default is 'pearson'. Other options include 'spearman', 'pls1c', and 'pls1w'. If a custom function that takes (gene_data, brain_data) as input is provided, the function uses the custom correlation method and sets `cor_method` to 'custom'.
#' @param aggre_method A character string specifying the aggregation method.
#'                     Default is 'mean'. Other options include 'median', 'meanabs',
#'                     'meansqr', 'maxmean', 'ks_orig', 'ks_weighted', 'ks_pos_neg_sum',
#'                     'sign_test', and 'rank_sum'. If a custom function that takes (geneList, geneSet) as input is provided, the function uses the custom aggregation method and sets `aggre_method` to 'custom'.
#' @param null_model A character string specifying the null model. Default is 'none', which calculates raw (empirical) scores. Other options include 'spin_brain', 'resample_gene', and 'coexp_matched'.
#' @param minGSSize An integer specifying the minimum gene set size after intersecting with the genes in gene_data. Default is 10.
#' @param maxGSSize An integer specifying the maximum gene set size after intersecting with the genes in gene_data. Default is 200.
#' @param n_cores An integer specifying the number of cores to use for parallel processing during permutation. Default is 0 (uses all available cores minus one).
#' @param n_perm An integer specifying the number of permutations for null models. Default is NULL (used when `null_model` is 'none').
#' @param perm_id A matrix of permutation indices for the 'spin_brain' null model. Default is NULL. Either perm_id or `coord.l`/`coord.r` must be provided when using the 'spin_brain' model.
#' @param coord.l A matrix of coordinates for the left hemisphere, used in the 'spin_brain' null model. Default is NULL. It can be NULL if `coord.r` or `perm_id` is provided.
#' @param coord.r A matrix of coordinates for the right hemisphere, used in the 'spin_brain' null model. Default is NULL. It can be NULL if `coord.l` or `perm_id` is provided.
#' @param seed An integer specifying the seed for reproducibility when using the 'spin_brain' model. Default is NULL.
#' @param matchcoexp_tol A numeric value specifying the tolerance for matched co-expression. Lower values result in better matching but require more iterations. Default is 0.05. See `resample_geneSetList_matching_coexp` for more details.
#' @param matchcoexp_max_iter An integer specifying the maximum number of iterations for matched co-expression. Default is 1,000,000. See `resample_geneSetList_matching_coexp` for more details.
#' @param verbose A logical indicating whether to print messages during processing. Default is TRUE.

#' @return A data frame containing the gene set scores with regions as rows and gene sets as columns.
#' @import pbapply
#' @import parallel
#' @export
brainscore <- function(brain_data,
                       gene_data,
                       annoData,
                       cor_method = c("pearson", "spearman", "pls1c", "pls1w"),
                       aggre_method = c("mean", "median", "meanabs", "meansqr", "maxmean", "ks_orig", "ks_weighted", "ks_pos_neg_sum", "sign_test", "rank_sum"),
                       null_model = c("none", "spin_brain", "resample_gene", "coexp_matched"),
                       minGSSize = 10,
                       maxGSSize = 200,
                       n_cores = 0,
                       n_perm = NULL,
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

  if (null_model == "none") {
    n_perm <- NULL
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

  if (null_model != "none") {
    if (n_cores == 0) {
      n_cores <- max(detectCores() - 1, 1) # Use all cores minus one, but ensure at least 1 core is used
    } else {
      n_cores <- min(n_cores, detectCores()) # Ensure n_cores does not exceed the number of available cores
    }
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

    # Initialize a cluster of workers
    cl <- if (n_cores > 1) makeCluster(n_cores) else NULL
    if (!is.null(cl)) {
      clusterExport(cl, c("brain_data", "gene_data", "n_perm", "perm_id", "aggregate_geneSet", "corr_brain_gene", "aggre_method", "selected.gs", "cor_method"), envir = environment())
    }

    gs.score <- pblapply(1:n_perm, function(idx) {
      null_brain_data <- brain_data[perm_id[, idx], , drop = FALSE]
      rownames(null_brain_data) <- rownames(brain_data)
      geneList.null <- corr_brain_gene(gene_data = gene_data, brain_data = null_brain_data, method = cor_method)
      gs_score.null <- aggregate_geneSetList(geneList.null, selected.gs, method = aggre_method, n_cores = 1)
      return(gs_score.null)
    }, cl = cl)

    if (!is.null(cl)) stopCluster(cl)
  } else if (null_model == "resample_gene") {
    if (verbose) {
      message("Aggregating gene set scores in resample_gene mode...")
    }

    # Initialize a cluster of workers
    cl <- if (n_cores > 1) makeCluster(n_cores) else NULL
    if (!is.null(cl)) {
      clusterExport(cl, c("n_perm", "geneList", "aggregate_geneSet", "aggre_method", "selected.gs"), envir = environment())
    }

    gs.score <- pblapply(1:n_perm, function(idx) {
      geneList.null <- geneList[sample(1:nrow(geneList), size = nrow(geneList), replace = FALSE), ]
      rownames(geneList.null) <- rownames(geneList)
      gs_score.null <- aggregate_geneSetList(geneList.null, selected.gs, method = aggre_method, n_cores = 1)
      return(gs_score.null)
    }, cl = cl)

    if (!is.null(cl)) stopCluster(cl)
  } else if (null_model == "coexp_matched") {
    if (verbose) {
      message("Generating null gene list with coexp_matched model...")
    }
    sampled_gs <- resample_geneSetList_matching_coexp(gene_data, selected.gs, tol = matchcoexp_tol, max_iter = matchcoexp_max_iter, n_perm = n_perm)

    if (is.null(sampled_gs)) stop("NULL sampled_gs returned.")
    if (!identical(names(selected.gs), names(sampled_gs))) stop("geneSetList and sampled_geneSetList are not matched.")

    if (verbose) {
      message("Aggregating gene set scores in coexp_matched mode...")
    }

    # Initialize a cluster of workers
    cl <- if (n_cores > 1) makeCluster(n_cores) else NULL
    if (!is.null(cl)) {
      clusterExport(cl, c("n_perm", "sampled_gs", "aggregate_geneSet", "aggre_method"), envir = environment())
    }

    gs.score <- pblapply(1:n_perm, function(idx) {
      sampled_gs_iter <- lapply(sampled_gs, function(x) x[[idx]])
      gs_score.null <- aggregate_geneSetList(geneList, sampled_gs_iter, method = aggre_method, n_cores = 1)
      return(gs_score.null)
    }, cl = cl)

    if (!is.null(cl)) stopCluster(cl)
  }

  # Add name to the list
  if (null_model != "none") {
    names(gs.score) <- paste0("null_", 1:n_perm)
  }

  # Add attributes to geneList
  attr(gs.score, "null_model") <- null_model
  attr(gs.score, "cor_method") <- cor_method
  attr(gs.score, "aggre_method") <- aggre_method
  attr(gs.score, "minGSSize") <- minGSSize
  attr(gs.score, "maxGSSize") <- maxGSSize
  attr(gs.score, "n_perm") <- n_perm
  return(gs.score)
}
