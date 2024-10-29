#' Perform Brain Score Linear Model Test
#'
#' This function performs a linear model test on brain score data with the option to use various null models for comparison.
#' It calculates gene set scores, performs linear modeling, calculates p-values, and identifies core genes.
#'
#' @param pred_df Data frame of predictor variables.
#' @param cov_df Data frame of covariate variables.
#' @param brain_data Data frame of brain imaging data.
#' @param gene_data Data frame of gene expression data.
#' @param annoData Environment containing annotation data.
#' @param gsScoreList.null Precomputed list of gene set scores for the null model by brainscore/brainscore.hpc function. Default is NULL.
#' @param cor_method Character string specifying the correlation method. Default is 'pearson'.
#'                   Other options include 'spearman', 'pls1c', 'pls1w', 'custom'.
#' @param aggre_method Character string specifying the aggregation method. Default is 'mean'.
#'                     Other options include 'median', 'meanabs', 'meansqr', 'maxmean', 'ks_orig', 'ks_weighted',
#'                     'ks_pos_neg_sum', 'sign_test', 'rank_sum', 'custom'.
#' @param null_model Character string specifying the null model method. Default is 'spin_brain'.
#'                   Other options include 'resample_gene', 'coexp_matched', 'none'.
#' @param minGSSize Integer specifying the minimum gene set size. Default is 10.
#' @param maxGSSize Integer specifying the maximum gene set size. Default is 200.
#' @param n_cores Integer specifying the number of cores to use for parallel processing. Default is 0.
#' @param n_perm Integer specifying the number of permutations. Default is 5000.
#' @param perm_id Optional permutation ID.
#' @param coord.l Optional left hemisphere coordinates.
#' @param coord.r Optional right hemisphere coordinates.
#' @param seed Optional random seed for generating perm_id.
#' @param threshold_type Character string specifying the threshold type for core genes. Default is 'sd'.
#'                       Other options include 'percentile'.
#' @param threshold_value Numeric value specifying the threshold level. Default is 1.
#' @param pvalueCutoff Numeric value specifying the p-value cutoff for significant results. Default is 0.05.
#' @param pAdjustMethod Character string specifying the method ("fdr","holm", "hochberg", "hommel", "bonferroni", "BH", "BY",  "none") for p-value adjustment. Default is 'fdr'. see p.adjust for more details.
#' @param padjCutoff Numeric value specifying the adjusted p-value cutoff for significant results. Default is NULL.
#' @param matchcoexp_tol Numeric value specifying the tolerance for matched coexpression. Default is 0.05.
#' @param matchcoexp_max_iter Integer specifying the maximum number of iterations for matched coexpression. Default is 1000000.
#' @param gsea_obj Logical specifying whether to return a GSEA object otherwise only a table will be returned. Default is TRUE.
#' @importFrom stats p.adjust
#' @importFrom utils getFromNamespace
#' @importFrom purrr list_transpose
#' @import parallel
#' @import pbapply
#' @importFrom dplyr select rename %>% everything
#' @importClassesFrom DOSE gseaResult
#' @return A data frame containing the results of the linear model test, including p-values, adjusted p-values,
#'         q-values, descriptions, and core genes.
#' @export
brainscore.lm_test <- function(pred_df,
                               cov_df,
                               brain_data,
                               gene_data,
                               annoData,
                               gsScoreList.null = NULL,
                               cor_method = c("pearson", "spearman", "pls1c", "pls1w", "custom"),
                               aggre_method = c("mean", "median", "meanabs", "meansqr", "maxmean", "ks_orig", "ks_weighted", "ks_pos_neg_sum", "sign_test", "rank_sum", "custom"),
                               null_model = c("spin_brain", "resample_gene", "coexp_matched", "none"),
                               minGSSize = 10,
                               maxGSSize = 200,
                               n_cores = 0,
                               n_perm = 5000,
                               perm_id = NULL,
                               coord.l = NULL,
                               coord.r = NULL,
                               seed = NULL,
                               threshold_type = c("sd", "percentile", "none"),
                               threshold_value = 1,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none"),
                               padjCutoff = NULL,
                               matchcoexp_tol = 0.05,
                               matchcoexp_max_iter = 1000000,
                               gsea_obj = TRUE) {
  # Validate arguments
  cor_method <- match.arg(cor_method)
  aggre_method <- match.arg(aggre_method)
  null_model <- match.arg(null_model)
  threshold_type <- match.arg(threshold_type)
  pAdjustMethod <- match.arg(pAdjustMethod)

  message("=========Empirical model======")
  # Generate true gene set scores
  gsScore.true <- brainscore(
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    cor_method = cor_method,
    aggre_method = aggre_method,
    null_model = "none",
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    n_cores = n_cores
  )
  dependent_df.true <- data.frame(gsScore.true, check.names = FALSE)
  message("Performing linear modeling using empirical gene set scores...")
  res <- simple_lm(dependent_df = dependent_df.true, pred_df = pred_df, cov_df = cov_df, stat2return = "all")
  if (null_model == "none") {
    message("Setting null model to none will return the empirical model results only. This is not recommended and only for quick testing purposes.")
    return(res)
  }
  stat.true <- simple_lm(dependent_df = dependent_df.true, pred_df = pred_df, cov_df = cov_df, stat2return = "tval_list")

  message("=========Null model======")
  # Generate null gene set scores
  if (is.null(gsScoreList.null)) {
    message("Computing null model...")
    gsScoreList.null <- brainscore(
      brain_data = brain_data,
      gene_data = gene_data,
      annoData = annoData,
      cor_method = cor_method,
      aggre_method = aggre_method,
      null_model = null_model,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      n_cores = n_cores,
      n_perm = n_perm,
      perm_id = perm_id,
      coord.l = coord.l,
      coord.r = coord.r,
      seed = seed,
      matchcoexp_tol = matchcoexp_tol,
      matchcoexp_max_iter = matchcoexp_max_iter
    )
  } else {
    # Check attributes of the precomputed null model
    null_model.precomp <- attr(gsScoreList.null, "null_model")
    cor_method.precomp <- attr(gsScoreList.null, "cor_method")
    aggre_method.precomp <- attr(gsScoreList.null, "aggre_method")
    minGSSize.precomp <- attr(gsScoreList.null, "minGSSize")
    maxGSSize.precomp <- attr(gsScoreList.null, "maxGSSize")
    n_perm.precomp <- attr(gsScoreList.null, "n_perm")

    # Check all attributes at once
    if (!((null_model.precomp == null_model) &&
      (cor_method.precomp == cor_method) &&
      (aggre_method.precomp == aggre_method) &&
      (minGSSize.precomp == minGSSize) &&
      (maxGSSize.precomp == maxGSSize) &&
      (n_perm.precomp == n_perm))) {
      message("Mismatches found between precomputed attributes and input variables.")
      message("Please check the following variables: null_model, cor_method, aggre_method, minGSSize, maxGSSize, n_perm.")
      stop("Please review the mismatches above.")
    } else {
      message("Using precomputed null model")
    }
  }

  message("Performing linear modeling using null gene set scores...")
  # Detect operating system and set number of cores appropriately

  if (n_cores == 0) {
    n_cores <- max(detectCores() - 1, 1)
  } else {
    n_cores <- min(n_cores, detectCores())
  }


  # Check for Windows and apply mclapply if not on Windows
  if (.Platform$OS.type != "windows") {
    message("Unix-like system detected. Using mclapply for parallel processing.")
    # Use mclapply for multi-core processing on Unix-like systems
    stat.tmp <- mclapply(seq_along(gsScoreList.null), function(i) {
      dependent_df.null <- data.frame(gsScoreList.null[[i]], check.names = FALSE)
      simple_lm(
        dependent_df = dependent_df.null,
        pred_df = pred_df,
        cov_df = cov_df,
        stat2return = "tval_list"
      )
    }, mc.cores = n_cores)
  } else {
    message("Windows detected. Using pblapply and makeCluster to do parallel processing. If large null GS score list is used, consider using a Unix-like system for more efficient processing.")

    cl <- if (n_cores > 1) makeCluster(n_cores) else NULL

    if (!is.null(cl)) {
      clusterExport(cl, list("gsScoreList.null", "pred_df", "cov_df", "simple_lm"),
        envir = environment()
      )
    }
    # Use parLapply for parallel processing
    stat.tmp <- pbapply(seq_along(gsScoreList.null), function(i) {
      dependent_df.null <- data.frame(gsScoreList.null[[i]], check.names = FALSE)
      simple_lm(
        dependent_df = dependent_df.null,
        pred_df = pred_df,
        cov_df = cov_df,
        stat2return = "tval_list"
      )
    }, cl = cl)
    # Stop the cluster after processing
    if (!is.null(cl)) stopCluster(cl)
  }

  # Combine results
  stat.null <- list_transpose(stat.tmp)

  # Calculate p-values
  message("Calculating p-values...")
  pvals <- calculate_pvals(stat.true, stat.null, method = "standard")
  calculate_qvalue <- getFromNamespace("calculate_qvalue", "DOSE")
  pvals.adj <- p.adjust(pvals, method = pAdjustMethod)
  qvals <- calculate_qvalue(pvals)

  # Prepare results
  message("Preparing results...")
  check_names <- all(
    names(stat.true) == names(stat.null),
    names(stat.null) == names(pvals),
    names(pvals) == names(pvals.adj),
    names(pvals.adj) == names(qvals),
    names(qvals) == res$Dependent_vars
  )

  if (!check_names) {
    stop("The names of the results are not consistent.")
  } else {
    res$Description <- get_termDescription(res$Dependent_vars, annoData)
    res$np.pval <- pvals
    res$np.padj <- pvals.adj
    res$np.qval <- qvals
    res$null_model <- null_model
  }

  # Filter significant results
  message("Filtering significant results...")
  res <- res[!is.na(res$np.pval), ]
  res <- res[res$np.pval <= pvalueCutoff, ]
  if (!is.null(padjCutoff)) {
    res <- res[res$np.padj <= padjCutoff, ]
  }
  res <- res[order(res$np.pval), ]

  if (nrow(res) == 0) {
    message("None of the gene sets are significant at the given p-value cutoff.")
  } else {
    geneList <- corr_brain_gene(brain_data = brain_data, gene_data = gene_data, method = cor_method)
    geneSetList <- get_geneSetList(annoData)
    selected.gs <- filter_geneSetList(rownames(geneList), geneSetList, minGSSize = minGSSize, maxGSSize = maxGSSize)
    survived.gs <- selected.gs[res$Dependent_vars]
    res$setSize <- sapply(survived.gs, length)
    if (threshold_type != "none" | gsea_obj) {
      message("Identifying core genes...")
      core_genes <- find_core_genes(geneList, survived.gs, pred_df = pred_df, cov_df = cov_df, aggre_method = aggre_method, n_cores = n_cores, threshold_type = threshold_type, threshold_value = threshold_value)
      res$core_genes <- sapply(core_genes, paste0, collapse = "/")
    }
  }

  if (gsea_obj) {
    params <- list(
      pvalueCutoff = pvalueCutoff,
      nPerm = n_perm,
      pAdjustMethod = pAdjustMethod,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      corMethod = cor_method,
      aggreMethod = aggre_method,
      nullType = null_model,
      thresType = threshold_type,
      thresVal = threshold_value
    )

    if (nrow(res) != 0) {
      res <- res %>%
        dplyr::rename(ID = .data$Dependent_vars) %>%
        dplyr::select(.data$ID, .data$Description, .data$setSize, everything()) %>%
        dplyr::select(-.data$p.val, -.data$p.adj) %>%
        dplyr::rename(
          pvalue = .data$np.pval,
          p.adjust = .data$np.padj,
          qvalue = .data$np.qval,
          core_enrichment = .data$core_genes
        )
    }
    res <- as.data.frame(res)
    rownames(res) <- res$ID
    message("Analysis complete.")
    return(new("gseaResult",
      result = res,
      geneSets = selected.gs,
      geneList = unlist(stat.true),
      permScores = as.matrix(do.call(rbind, stat.null)),
      params = params,
      keytype = "SYMBOL",
      readable = TRUE
    ))
  } else {
    return(res)
  }
}
