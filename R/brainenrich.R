#' Perform Brain Gene Set Analysis
#'
#' This function performs a gene set analysis using brain data.
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
#'                    'sign_test', 'rank_sum', 'custom'.
#' @param null_model A character string specifying the null model.
#'                   Default is 'spin_brain'. Other options include 'resample_gene',
#'                   'coexp_matched'.
#' @param n_perm An integer specifying the number of permutations. Default is 5000.
#' @param perm_id A matrix of permutation IDs. Default is NULL.
#' @param coord.l A matrix of left hemisphere coordinates. Default is NULL.
#' @param coord.r A matrix of right hemisphere coordinates. Default is NULL.
#' @param seed An integer specifying the seed for reproducibility of spinning brain. Default is NULL.
#' @param n_cores An integer specifying the number of cores to use. Default is 0 (use all cores - 1).
#' @param minGSSize An integer specifying the minimum gene set size. Default is 10.
#' @param maxGSSize An integer specifying the maximum gene set size. Default is 200.
#' @param threshold_type A character string specifying the threshold type for core genes.
#'                   Default is 'sd'. Other option is 'percentile'.
#' @param threshold_value A numeric value specifying the threshold value for core genes. Default is 1.
#' @param pvalueCutoff A numeric value specifying the p-value cutoff for output. Default is 0.05.
#' @param pAdjustMethod A character string specifying the method ("fdr","holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none") for p-value adjustment.
#'                      Default is 'fdr'. see p.adjust for more details.
#' @param matchcoexp_tol A numeric value specifying the tolerance for matched co-expression.
#'                       Lower value means better matching but will take much more iterations. Default is 0.05.
#' @param matchcoexp_max_iter An integer specifying the maximum number of iterations
#'                            for matched co-expression. Default is 1000000.
#' @importClassesFrom DOSE gseaResult
#' @import DOSE
#' @importFrom utils getFromNamespace
#' @return A gseaResult object containing the enrichment results.
#' @export
brainenrich <- function(brain_data,
                        gene_data,
                        annoData,
                        cor_method = c("pearson", "spearman", "pls1c", "pls1w", "custom"),
                        aggre_method = c("mean", "median", "meanabs", "meansqr", "maxmean", "ks_orig", "ks_weighted", "ks_pos_neg_sum", "sign_test", "rank_sum", "custom"),
                        null_model = c("spin_brain", "resample_gene", "coexp_matched"),
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
                        matchcoexp_tol = 0.05,
                        matchcoexp_max_iter = 1000000) {
  # Check inputs
  stopifnot(is.environment(annoData))

  if (null_model == "spin_brain" && is.null(perm_id) && is.null(coord.l) && is.null(coord.r)) {
    stop("For null_model 'spin_brain', 'perm_id' or at least one of 'coord.l' or 'coord.r' must be provided.")
  }

  if (null_model == "spin_brain" && !is.null(perm_id)) {
    if (n_perm > ncol(perm_id)) {
      stop("The number of permutations must be less than or equal to the number of columns in 'perm_id'.")
    } else if (n_perm < ncol(perm_id)) {{ perm_id <- perm_id[, 1:n_perm] }}
  }

  if (!identical(rownames(gene_data), rownames(brain_data))) {
    stop("Rownames of 'gene_data' and 'brain_data' must be identical.")
  }

  if (ncol(brain_data) != 1) {
    stop("This function is designed for group-level analysis and supports one column of brain data. Consider using s/lapply.")
  }

  cor_method <- match.arg(cor_method)
  aggre_method <- match.arg(aggre_method)
  null_model <- match.arg(null_model)
  threshold_type <- match.arg(threshold_type)
  pAdjustMethod <- match.arg(pAdjustMethod)

  # Perform analysis
  message("Calculating true gene-brain correlations...")
  geneList.true <- corr_brain_gene(gene_data = gene_data, brain_data = brain_data, method = cor_method)
  message("Generating gene set list from annotation data...")
  geneSetList <- get_geneSetList(annoData)
  message("Filtering gene set list...")
  selected.gs <- filter_geneSetList(rownames(geneList.true), geneSetList, minGSSize = minGSSize, maxGSSize = maxGSSize)
  message("Number of gene sets left: ", length(selected.gs))
  message("Aggregating true gene set scores...")
  gs_score.true <- aggregate_geneSetList(geneList.true, selected.gs, method = aggre_method, n_cores = n_cores)

  if (null_model == "spin_brain") {
    message("Generating null brain data with spin_brain model...")
    if (is.null(perm_id)) {
      perm_id <- rotate_parcellation(coord.l = coord.l, coord.r = coord.r, nrot = n_perm, seed = seed)
    }
    null_brain_data <- generate_null_brain_data(brain_data, perm_id)
    geneList.null <- corr_brain_gene(gene_data = gene_data, brain_data = null_brain_data, method = cor_method)
    message("Aggregating true gene set scores...")
    gs_score.null <- aggregate_geneSetList(geneList.null, selected.gs, method = aggre_method, n_cores = n_cores)
  } else if (null_model == "resample_gene") {
    message("Generating null gene list with resample_gene model...")
    geneList.null <- resample_gene(geneList.true, n_perm = n_perm)
    message("Aggregating true gene set scores...")
    gs_score.null <- aggregate_geneSetList(geneList.null, selected.gs, method = aggre_method, n_cores = n_cores)
  } else if (null_model == "coexp_matched") {
    message("Generating null gene list with coexp_matched model...")
    sampled_gs <- resample_geneSetList_matching_coexp(gene_data, selected.gs, tol = matchcoexp_tol, max_iter = matchcoexp_max_iter, n_perm = n_perm, n_cores = n_cores)
    message("Aggregating true gene set scores...")
    gs_score.null <- aggregate_geneSetList_matching_coexp(geneList.true, selected.gs, sampled_gs, method = aggre_method, n_cores = n_cores)
  }

  message("Calculating p-values...")
  if (aggre_method %in% c("ks_orig", "ks_weighted")) {
    pvals <- calculate_pvals(gs_score.true, gs_score.null, method = c("split_pos_neg"))
  } else {
    pvals <- calculate_pvals(gs_score.true, gs_score.null, method = c("standard"))
  }
  calculate_qvalue <- getFromNamespace("calculate_qvalue", "DOSE")
  pvals.adj <- p.adjust(pvals, method = pAdjustMethod)
  qvals <- calculate_qvalue(pvals)

  message("Preparing results...")
  gs.name <- names(selected.gs)
  TERM2NAME <- getFromNamespace("TERM2NAME", "DOSE")
  Description <- TERM2NAME(gs.name, annoData)

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

  res <- data.frame(
    ID = as.character(gs.name),
    Description = as.character(Description),
    setSize = sapply(selected.gs, length),
    gsScore = unlist(gs_score.true),
    pvalue = pvals,
    p.adjust = pvals.adj,
    qvalue = qvals,
    stringsAsFactors = FALSE
  )

  message("Filtering significant results...")
  res <- res[!is.na(res$pvalue), ]
  res <- res[res$pvalue <= pvalueCutoff, ]
  res <- res[res$p.adjust <= pvalueCutoff, ]
  res <- res[order(res$pvalue), ]

  if (nrow(res) == 0) {
    message("No term enriched under specific pvalueCutoff...")
    return(new("gseaResult",
      result = res,
      geneSets = selected.gs,
      geneList = geneList.true[, 1],
      permScores = as.matrix(do.call(rbind, gs_score.null)),
      params = params,
      keytype = "SYMBOL",
      readable = FALSE
    ))
  } else {
    if (threshold_type != "none") {
      message("Identifying core genes...")
      survived.gs <- selected.gs[res$ID]
      core_genes <- find_core_genes(geneList.true, survived.gs, aggre_method = aggre_method, n_cores = n_cores, threshold_type = threshold_type, threshold_value = threshold_value)
      res$core_enrichment <- sapply(core_genes, paste0, collapse = "/")
    }

    message("Analysis complete.")
    return(new("gseaResult",
      result = res,
      geneSets = selected.gs,
      geneList = geneList.true[, 1],
      permScores = as.matrix(do.call(rbind, gs_score.null)),
      params = params,
      keytype = "SYMBOL",
      readable = TRUE
    ))
  }
}
