#' Resample Gene List
#'
#' Generates a set of permuted gene lists from the original gene list,
#' ensuring uniqueness in the permuted sets.
#'
#' @param geneList.true A matrix of gene expression data.
#' @param n_perm Number of permutations to generate.
#' @return A matrix of permuted gene lists.
#' @export
resample_gene <- function(geneList.true, n_perm = 5000) {
  # Check if the number of permutations requested exceeds 10,000. If it does, stop the function and display an error message.

  if (!is.matrix(geneList.true) || ncol(geneList.true) != 1) {
    stop("geneList.true should be a m x 1 matrix. Please include drop=FALSE when subsetting.")
  }

  # Ensure the gene names are present
  gnames <- rownames(geneList.true)
  if (is.null(gnames)) {
    stop("geneList.true must have row names.")
  }
  # Create a matrix 'geneList.null' by replicating gene lists. This is done 'n_perm + 100' times. Each replication involves shuffling the original gene list without replacement.
  # 'simplify = TRUE' ensures the result is a matrix.
  geneList.null <- replicate(n = n_perm + 100, expr = {
    sample(geneList.true, size = nrow(geneList.true), replace = FALSE)
  }, simplify = TRUE)

  # Remove duplicate columns from the transposed 'geneList.null' matrix. This is to ensure that all permutations are unique.
  geneList.null <- geneList.null[, !duplicated(t(geneList.null))]

  # Trim the matrix to keep only the first 'n_perm' columns, in case more than 'n_perm' unique permutations were generated initially.
  geneList.null <- geneList.null[, 1:n_perm]

  # Set the row names of 'geneList.null' to match those of 'geneList.true' to maintain consistency in identifiers.
  rownames(geneList.null) <- rownames(geneList.true)

  # Create column names for 'geneList.null' that indicate these are null permutations, e.g., 'null_1', 'null_2', ..., 'null_n_perm'.
  colnames(geneList.null) <- paste0("null_", seq_len(n_perm))

  # Copy attributes from 'geneList.true' to 'geneList.null' to retain metadata, such as dimension names, class, or other user-defined attributes.
  attr(geneList.null, "is_fisherz") <- attr(geneList.true, "is_fisherz")
  attr(geneList.null, "n.region") <- attr(geneList.true, "n.region")

  # Return the matrix of null permutations.
  return(geneList.null)
}


#' Resample Gene Sets with Specified Constraints
#'
#' This function resamples gene sets based on specific constraints like matching
#' co-expression patterns. The methodology implemented is informed by Wei et al. (2022) on statistical testing in transcriptomic-neuroimaging studies.
#' It is important to note that restricting null models to a subset of genes can be problematic. The empirical statistics sampled from
#' the full gene pool differ from those derived from a restricted pool. Therefore, usage of this approach should be
#' with caution.
#'
#' @param gene_data A matrix or data frame representing gene expression data.
#' @param geneSetList A list of gene sets to be resampled.
#' @param tol A numeric value indicating the tolerance for matching co-expression patterns (default = 0.01).
#' @param max_iter An integer indicating the maximum number of iterations for the sampling process (default = 1000000).
#' @param n_perm An integer indicating the number of permutations to generate (default = 5000).
#' @param n_cores An integer indicating the number of cores to use for parallel processing (default = 1).
#' @return A list of resampled gene sets based on the specified constraints.
#' @import pbapply
#' @import parallel
#' @importFrom stats cor
#' @export
#' @references
#' Wei, Y., de Lange, S. C., Pijnenburg, R., Scholtens, L. H., Ardesch, D. J., Watanabe, K., Posthuma, D., & van den Heuvel, M. P. (2022).
#' Statistical testing in transcriptomic-neuroimaging studies: A how-to and evaluation of methods assessing spatial and gene specificity.
#' Human Brain Mapping, 43(3), 885–901.
resample_geneSetList_matching_coexp <- function(gene_data, geneSetList, tol = 0.01, max_iter = 1000000, n_perm = 5000, n_cores = 1) {
  if (!ask_user_continue("Resampling gene sets with coexpression matched will take a long time.")) {
    cat("Operation aborted by the user.\n")
    return(NULL)
  }
  # Calculate the co-expression matrix
  coexp_matrix <- cor(gene_data)
  total_gs <- length(geneSetList)

  # Determine the number of cores to use
  if (n_cores == 0) {
    n_cores <- max(detectCores() - 1, 1) # Use all cores minus one, but ensure at least 1 core is used
  } else {
    n_cores <- min(n_cores, detectCores()) # Ensure n_cores does not exceed the number of available cores
  }

  # Initialize a cluster of workers
  cl <- makeCluster(n_cores)

  # Export necessary variables and functions to the cluster
  clusterExport(cl, c(
    "geneSetList", "coexp_matrix", "tol", "max_iter",
    "n_perm", "sample_gs_matching_coexp", "total_gs"
  ),
  envir = environment()
  )

  # Parallelize the processing using pblapply for progress bar
  sampled_geneSetList <- pblapply(seq_along(geneSetList), function(i) {
    gs <- geneSetList[[i]]
    cat(sprintf("Sampling gene set %d/%d: %s (gs size: %d) \n", i, total_gs, names(geneSetList)[i], length(gs)))
    sample_gs_matching_coexp(gs = gs, coexp_matrix = coexp_matrix, tol = tol, max_iter = max_iter, n_target = n_perm)
  }, cl = cl)

  # Stop the cluster after processing
  stopCluster(cl)

  names(sampled_geneSetList) <- names(geneSetList)

  return(sampled_geneSetList)
}




#' Sample Gene Sets Matching Co-Expression
#'
#' This function samples gene sets that closely match the co-expression profile of a target gene set.
#' The methodology implemented is informed by Wei et al. (2022) on statistical testing in transcriptomic-neuroimaging studies.
#'
#' @param gs The target gene set for which similar co-expression profiles are sought.
#' @param coexp_matrix A co-expression matrix, typically calculated as the correlation matrix of gene expression data.
#' @param tol Tolerance for the difference between the co-expression of the target and the sampled gene sets.
#' @param max_iter Maximum number of iterations to attempt finding matches.
#' @param n_target Number of gene sets to sample.
#' @return A list of gene sets that closely match the target gene set's co-expression profile.
#' @references
#' Wei, Y., de Lange, S. C., Pijnenburg, R., Scholtens, L. H., Ardesch, D. J., Watanabe, K., Posthuma, D., & van den Heuvel, M. P. (2022).
#' Statistical testing in transcriptomic-neuroimaging studies: A how-to and evaluation of methods assessing spatial and gene specificity.
#' Human Brain Mapping, 43(3), 885–901.


sample_gs_matching_coexp <- function(gs, coexp_matrix, tol = 0.01, max_iter = 1000000, n_target = 5000) {
  # Calculate the mean co-expression value for the lower triangle of the co-expression matrix of the target gene set



  gs_coexp_matrix <- coexp_matrix[gs, gs]
  gs_coexp_lower <- gs_coexp_matrix[lower.tri(gs_coexp_matrix)]
  target_coexp <- mean(gs_coexp_lower)

  # Initialize a list to hold sampled gene sets that meet the criteria
  sampled_gs <- list()
  cat(sprintf("Target coexp: %f \n", target_coexp))
  for (i in 1:max_iter) {
    if (i %% 10000 == 0) {
      cat(sprintf("Iteration %d \n", i))
    }
    sampled_genes <- sample(colnames(coexp_matrix), length(gs), replace = FALSE)
    sampled_coexp_matrix <- coexp_matrix[sampled_genes, sampled_genes]
    sampled_coexp_lower <- sampled_coexp_matrix[lower.tri(sampled_coexp_matrix)]
    sampled_coexp_avg <- mean(sampled_coexp_lower)

    if (abs(sampled_coexp_avg - target_coexp) < tol) {
      sampled_gs[[length(sampled_gs) + 1]] <- sampled_genes
    }

    # Remove duplicates and check if the target number of sets has been reached
    sampled_gs <- unique(sampled_gs)
    if (length(sampled_gs) == n_target) {
      break
    }
  }

  if (length(sampled_gs) < n_target) {
    stop("n_target not reached, consider increasing max_iter or decreasing tol")
  }

  return(sampled_gs)
}



#' Swap Gene List
#'
#' This function swaps the original values of a gene set (`orig_gs`) with the sampled gene sets (`sampled_gs`) within the given `geneList.true`.
#' The resulting `geneList.null` is in the same format as generated by other approaches.
#'
#' @param geneList.true A matrix representing the true gene list, with dimensions m x 1.
#' @param orig_gs A vector of original gene set identifiers.
#' @param sampled_gs A list of sampled gene sets, each being a vector of gene identifiers.
#' @return A matrix where each column represents a null gene list generated by swapping `orig_gs` with each set in `sampled_gs`.

swap_geneList <- function(geneList.true, orig_gs, sampled_gs) {
  # Check if geneList.true is a matrix with one column
  if (!is.matrix(geneList.true) || ncol(geneList.true) != 1) {
    stop("geneList.true should be a m x 1 matrix. Please include drop=FALSE when subsetting.")
  }

  # Ensure the gene names are present
  gnames <- rownames(geneList.true)
  if (is.null(gnames)) {
    stop("geneList.true must have row names.")
  }

  # Extract original values from geneList.true for the original gene set
  origVals <- as.numeric(geneList.true[gnames %in% orig_gs, ])

  # Create a list of null gene lists by swapping original values with sampled gene sets
  null_list <- lapply(sampled_gs, function(gs_i) {
    tmp_null <- geneList.true
    swapVals <- as.numeric(geneList.true[gnames %in% gs_i, ])
    tmp_null[gnames %in% gs_i, ] <- origVals
    tmp_null[gnames %in% orig_gs, ] <- swapVals
    return(tmp_null)
  })

  # Combine the list into a matrix
  null_geneList <- do.call(cbind, null_list)
  colnames(null_geneList) <- paste0("null_", 1:ncol(null_geneList))

  # Preserve attributes
  attr(null_geneList, "is_fisherz") <- attr(geneList.true, "is_fisherz")
  attr(null_geneList, "n.region") <- attr(geneList.true, "n.region")

  return(null_geneList)
}
