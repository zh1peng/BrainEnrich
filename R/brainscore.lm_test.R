#' Perform linear model testing on brain scores
#'
#' This function computes empirical and null-model brain gene-set scores, runs
#' linear models at the term level, and performs non-parametric significance
#' testing using permutation-derived null statistics.
#'
#' @param pred_df Data frame of predictor variables.
#' @param cov_df Data frame of covariate variables.
#' @param brain_data Data frame/matrix of brain imaging data.
#' @param gene_data Data frame/matrix of gene expression data.
#' @param annoData Annotation object of class `EnrichAnno`.
#' @param gsScoreList.null Optional precomputed null gene-set score list.
#' @param cor_method Correlation method name or custom function.
#' @param aggre_method Aggregation method name or custom function.
#' @param null_model Null model method.
#' @param minGSSize Minimum retained gene set size.
#' @param maxGSSize Maximum retained gene set size.
#' @param n_cores Number of cores.
#' @param n_perm Number of permutations.
#' @param perm_id Optional permutation ID matrix.
#' @param coord.l Optional left hemisphere coordinates.
#' @param coord.r Optional right hemisphere coordinates.
#' @param seed Optional seed.
#' @param threshold_type Core-gene threshold type.
#' @param threshold_value Core-gene threshold value.
#' @param pvalueCutoff Non-parametric p-value cutoff.
#' @param pAdjustMethod P-value adjustment method.
#' @param padjCutoff Optional adjusted p-value cutoff.
#' @param matchcoexp_tol Co-expression matching tolerance.
#' @param matchcoexp_max_iter Co-expression matching max iterations.
#' @param gsea_obj Logical; if TRUE return `EnrichRes`, else return filtered table.
#' @importFrom stats p.adjust
#' @importFrom purrr list_transpose
#' @import parallel
#' @import pbapply
#' @importFrom dplyr select rename %>% everything
#' @importFrom rlang .data
#' @return An `EnrichRes` object when `gsea_obj = TRUE`, otherwise a data frame.
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
  if (is.environment(annoData)) {
    annoData <- anno_env_to_EnrichAnno(annoData)
  }
  if (!is_EnrichAnno(annoData)) {
    stop("annoData must be an EnrichAnno object.")
  }

  if (!is.function(cor_method)) {
    cor_method <- match.arg(cor_method)
  }
  if (!is.function(aggre_method)) {
    aggre_method <- match.arg(aggre_method)
  }
  null_model <- match.arg(null_model)
  threshold_type <- match.arg(threshold_type)
  pAdjustMethod <- match.arg(pAdjustMethod)

  message("=========Empirical model======")
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
  model_table <- simple_lm(dependent_df = dependent_df.true, pred_df = pred_df, cov_df = cov_df, stat2return = "all")

  if (null_model == "none") {
    if (!gsea_obj) {
      return(model_table)
    }
    geneList <- corr_brain_gene(brain_data = brain_data, gene_data = gene_data, method = cor_method)
    geneSetList <- get_geneSetList(annoData)
    selected.gs <- filter_geneSetList(rownames(geneList), geneSetList, minGSSize = minGSSize, maxGSSize = maxGSSize)
    return(new_EnrichRes(
      analysis_type = "brainscore_lm_test",
      table = model_table,
      gene_sets = selected.gs,
      perm_scores = NULL,
      params = list(
        pvalueCutoff = pvalueCutoff,
        nPerm = n_perm,
        pAdjustMethod = pAdjustMethod,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        corMethod = if (is.function(cor_method)) "custom" else cor_method,
        aggreMethod = if (is.function(aggre_method)) "custom" else aggre_method,
        nullType = null_model
      ),
      diagnostics = list(n_results = nrow(model_table), null_model = null_model),
      model_table = model_table
    ))
  }

  stat.true <- simple_lm(dependent_df = dependent_df.true, pred_df = pred_df, cov_df = cov_df, stat2return = "tval_list")

  message("=========Null model======")
  if (is.null(gsScoreList.null)) {
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
    null_model.precomp <- attr(gsScoreList.null, "null_model")
    cor_method.precomp <- attr(gsScoreList.null, "cor_method")
    aggre_method.precomp <- attr(gsScoreList.null, "aggre_method")
    minGSSize.precomp <- attr(gsScoreList.null, "minGSSize")
    maxGSSize.precomp <- attr(gsScoreList.null, "maxGSSize")
    n_perm.precomp <- attr(gsScoreList.null, "n_perm")
    cor_label <- if (is.function(cor_method)) "custom" else cor_method
    aggre_label <- if (is.function(aggre_method)) "custom" else aggre_method

    if (!((null_model.precomp == null_model) &&
      (cor_method.precomp == cor_label) &&
      (aggre_method.precomp == aggre_label) &&
      (minGSSize.precomp == minGSSize) &&
      (maxGSSize.precomp == maxGSSize) &&
      (n_perm.precomp == n_perm))) {
      stop("Precomputed gsScoreList.null attributes do not match current inputs.")
    }
  }

  message("Performing linear modeling using null gene set scores...")
  n_cores <- be_resolve_n_cores(n_cores)
  stat.tmp <- be_parallel_lapply(seq_along(gsScoreList.null), function(i) {
    dependent_df.null <- data.frame(gsScoreList.null[[i]], check.names = FALSE)
    simple_lm(dependent_df = dependent_df.null, pred_df = pred_df, cov_df = cov_df, stat2return = "tval_list")
  },
  n_cores = n_cores,
  seed = seed,
  export = c("gsScoreList.null", "pred_df", "cov_df", "simple_lm"),
  envir = environment()
  )

  stat.null <- list_transpose(stat.tmp)
  pvals <- calculate_pvals(stat.true, stat.null, method = "standard")
  pvals.adj <- p.adjust(pvals, method = pAdjustMethod)
  qvals <- calculate_qvalue_local(pvals)

  check_names <- all(
    names(stat.true) == names(stat.null),
    names(stat.null) == names(pvals),
    names(pvals) == names(pvals.adj),
    names(pvals.adj) == names(qvals),
    names(qvals) == model_table$Dependent_vars
  )
  if (!check_names) {
    stop("The names of linear-model and null-model outputs are inconsistent.")
  }

  model_table$Description <- get_termDescription(model_table$Dependent_vars, annoData)
  model_table$np.pval <- pvals
  model_table$np.padj <- pvals.adj
  model_table$np.qval <- qvals
  model_table$null_model <- null_model

  res <- model_table[!is.na(model_table$np.pval), , drop = FALSE]
  res <- res[res$np.pval <= pvalueCutoff, , drop = FALSE]
  if (!is.null(padjCutoff)) {
    res <- res[res$np.padj <= padjCutoff, , drop = FALSE]
  }
  res <- res[order(res$np.pval), , drop = FALSE]

  geneList <- corr_brain_gene(brain_data = brain_data, gene_data = gene_data, method = cor_method)
  geneSetList <- get_geneSetList(annoData)
  selected.gs <- filter_geneSetList(rownames(geneList), geneSetList, minGSSize = minGSSize, maxGSSize = maxGSSize)

  core_genes_details <- NULL
  if (nrow(res) > 0) {
    survived.gs <- selected.gs[res$Dependent_vars]
    res$setSize <- vapply(survived.gs, length, FUN.VALUE = integer(1))
    if (threshold_type != "none") {
      core_genes_details <- find_core_genes(
        geneList, survived.gs,
        pred_df = pred_df, cov_df = cov_df,
        aggre_method = aggre_method, n_cores = n_cores,
        threshold_type = threshold_type, threshold_value = threshold_value
      )
      res$core_genes <- vapply(core_genes_details, paste0, FUN.VALUE = character(1), collapse = "/")
    }
  } else {
    message("None of the gene sets are significant at the given p-value cutoff.")
  }

  if (!gsea_obj) {
    return(res)
  }

  # Keep table naming aligned with enrichment output conventions.
  out_table <- res
  if (nrow(out_table) > 0) {
    out_table <- out_table %>%
      dplyr::rename(ID = .data$Dependent_vars) %>%
      dplyr::select(.data$ID, .data$Description, .data$setSize, dplyr::everything()) %>%
      dplyr::select(-.data$p.val, -.data$p.adj)
    if ("core_genes" %in% names(out_table)) {
      out_table <- out_table %>%
        dplyr::rename(
          pvalue = .data$np.pval,
          p.adjust = .data$np.padj,
          qvalue = .data$np.qval,
          core_enrichment = .data$core_genes
        )
    } else {
      out_table <- out_table %>%
        dplyr::rename(
          pvalue = .data$np.pval,
          p.adjust = .data$np.padj,
          qvalue = .data$np.qval
        )
    }
  } else {
    out_table <- data.frame()
  }

  params <- list(
    pvalueCutoff = pvalueCutoff,
    nPerm = n_perm,
    pAdjustMethod = pAdjustMethod,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    corMethod = if (is.function(cor_method)) "custom" else cor_method,
    aggreMethod = if (is.function(aggre_method)) "custom" else aggre_method,
    nullType = null_model,
    thresType = threshold_type,
    thresVal = threshold_value
  )

  new_EnrichRes(
    analysis_type = "brainscore_lm_test",
    table = out_table,
    gene_sets = selected.gs,
    perm_scores = as.matrix(do.call(rbind, stat.null)),
    params = params,
    diagnostics = list(n_results = nrow(out_table), null_model = null_model),
    model_table = model_table,
    np_table = data.frame(ID = names(pvals), np.pval = pvals, np.padj = pvals.adj, np.qval = qvals, stringsAsFactors = FALSE),
    core_genes = core_genes_details,
    gene_scores = unlist(stat.true)
  )
}
