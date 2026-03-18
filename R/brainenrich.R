#' Perform brain gene set analysis
#'
#' This function performs group-level gene set analysis by aggregating
#' gene-brain association values within predefined gene sets and comparing
#' empirical scores against null-model scores.
#'
#' @param brain_data A data frame/matrix of brain measurements with regions as rows.
#' @param gene_data A data frame/matrix of gene expression with regions as rows and genes as columns.
#' @param annoData An annotation object of class `EnrichAnno`.
#' @param cor_method Correlation/association method name or a custom function.
#' @param aggre_method Aggregation method name or a custom function.
#' @param null_model Null model type.
#' @param n_perm Number of permutations.
#' @param perm_id Optional permutation index matrix for `spin_brain`.
#' @param coord.l Optional left hemisphere coordinates.
#' @param coord.r Optional right hemisphere coordinates.
#' @param seed Optional seed.
#' @param n_cores Number of cores.
#' @param minGSSize Minimum gene set size.
#' @param maxGSSize Maximum gene set size.
#' @param threshold_type Core-gene threshold type.
#' @param threshold_value Core-gene threshold value.
#' @param pvalueCutoff P-value cutoff for reporting.
#' @param pAdjustMethod P-value adjust method.
#' @param padjCutoff Optional adjusted p-value cutoff.
#' @param matchcoexp_tol Co-expression matching tolerance.
#' @param matchcoexp_max_iter Co-expression matching max iterations.
#' @return An object of class `EnrichRes` with `analysis_type = "brainenrich"`.
#' @export
brainenrich <- function(brain_data,
                        gene_data,
                        annoData,
                        cor_method = c("pearson", "spearman", "pls1c", "pls1w"),
                        aggre_method = c("mean", "median", "meanabs", "meansqr", "maxmean", "ks_orig", "ks_weighted", "ks_pos_neg_sum", "sign_test", "rank_sum"),
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
                        padjCutoff = NULL,
                        matchcoexp_tol = 0.05,
                        matchcoexp_max_iter = 1000000) {
  if (is.environment(annoData)) {
    annoData <- anno_env_to_EnrichAnno(annoData)
  }
  if (!is_EnrichAnno(annoData)) {
    stop("annoData must be an EnrichAnno object.")
  }

  if (!identical(rownames(gene_data), rownames(brain_data))) {
    stop("Rownames of 'gene_data' and 'brain_data' must be identical.")
  }
  if (ncol(brain_data) != 1) {
    stop("This function is for group-level data and supports one brain data column.")
  }

  null_model <- match.arg(null_model)
  threshold_type <- match.arg(threshold_type)
  pAdjustMethod <- match.arg(pAdjustMethod)

  cor_method_label <- NULL
  aggre_method_label <- NULL
  if (is.function(cor_method)) {
    cor_method_label <- "custom"
  } else {
    cor_method <- match.arg(cor_method)
    cor_method_label <- cor_method
  }
  if (is.function(aggre_method)) {
    aggre_method_label <- "custom"
  } else {
    aggre_method <- match.arg(aggre_method)
    aggre_method_label <- aggre_method
  }

  if (null_model == "spin_brain" && is.null(perm_id) && is.null(coord.l) && is.null(coord.r)) {
    stop("For null_model 'spin_brain', provide perm_id or at least one of coord.l/coord.r.")
  }
  if (null_model == "spin_brain" && !is.null(perm_id)) {
    if (n_perm > ncol(perm_id)) {
      stop("n_perm must not exceed number of columns in perm_id.")
    } else if (n_perm < ncol(perm_id)) {
      perm_id <- perm_id[, 1:n_perm, drop = FALSE]
    }
  }

  message("Calculating true gene-brain correlations...")
  geneList.true <- corr_brain_gene(gene_data = gene_data, brain_data = brain_data, method = cor_method)
  geneSetList <- get_geneSetList(annoData)
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
    message("Aggregating null gene set scores...")
    gs_score.null <- aggregate_geneSetList(geneList.null, selected.gs, method = aggre_method, n_cores = n_cores)
  } else if (null_model == "resample_gene") {
    message("Generating null gene list with resample_gene model...")
    geneList.null <- resample_gene(geneList.true, n_perm = n_perm)
    message("Aggregating null gene set scores...")
    gs_score.null <- aggregate_geneSetList(geneList.null, selected.gs, method = aggre_method, n_cores = n_cores)
  } else {
    message("Generating null gene list with coexp_matched model...")
    sampled_gs <- resample_geneSetList_matching_coexp(
      gene_data, selected.gs,
      tol = matchcoexp_tol, max_iter = matchcoexp_max_iter,
      n_perm = n_perm, n_cores = n_cores
    )
    message("Aggregating null gene set scores...")
    gs_score.null <- aggregate_geneSetList_matching_coexp(
      geneList.true, selected.gs, sampled_gs,
      method = aggre_method, n_cores = n_cores
    )
  }

  message("Calculating p-values...")
  if (!is.function(aggre_method) && aggre_method %in% c("ks_orig", "ks_weighted")) {
    pvals <- calculate_pvals(gs_score.true, gs_score.null, method = "split_pos_neg")
  } else {
    pvals <- calculate_pvals(gs_score.true, gs_score.null, method = "standard")
  }
  pvals.adj <- stats::p.adjust(pvals, method = pAdjustMethod)
  qvals <- calculate_qvalue_local(pvals)

  gs.name <- names(selected.gs)
  Description <- get_termDescription(gs.name, annoData)
  params <- list(
    pvalueCutoff = pvalueCutoff,
    nPerm = n_perm,
    pAdjustMethod = pAdjustMethod,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    corMethod = cor_method_label,
    aggreMethod = aggre_method_label,
    nullType = null_model,
    thresType = threshold_type,
    thresVal = threshold_value
  )

  res <- data.frame(
    ID = as.character(gs.name),
    Description = as.character(Description),
    setSize = vapply(selected.gs, length, FUN.VALUE = integer(1)),
    gsScore = unlist(gs_score.true),
    pvalue = pvals,
    p.adjust = pvals.adj,
    qvalue = qvals,
    stringsAsFactors = FALSE
  )

  res <- res[!is.na(res$pvalue), , drop = FALSE]
  res <- res[res$pvalue <= pvalueCutoff, , drop = FALSE]
  if (!is.null(padjCutoff)) {
    res <- res[res$p.adjust <= padjCutoff, , drop = FALSE]
  }
  res <- res[order(res$pvalue), , drop = FALSE]

  core_genes <- NULL
  if (nrow(res) > 0 && threshold_type != "none") {
    message("Identifying core genes...")
    survived.gs <- selected.gs[res$ID]
    core_genes <- find_core_genes(
      geneList.true, survived.gs,
      aggre_method = aggre_method, n_cores = n_cores,
      threshold_type = threshold_type, threshold_value = threshold_value
    )
    res$core_enrichment <- vapply(core_genes, paste0, FUN.VALUE = character(1), collapse = "/")
  }

  message("Analysis complete.")
  new_EnrichRes(
    analysis_type = "brainenrich",
    table = res,
    gene_sets = selected.gs,
    perm_scores = as.matrix(do.call(rbind, gs_score.null)),
    params = params,
    diagnostics = list(
      n_gene_sets = length(selected.gs),
      n_results = nrow(res),
      null_model = null_model
    ),
    core_genes = core_genes,
    gene_scores = geneList.true
  )
}
