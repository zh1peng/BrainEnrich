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
#' @import DOSE
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

  if (null_model == "none") {
    message("Aggregating gene set scores...")
    gs.score <- aggregate_geneSetList(geneList, selected.gs, method = aggre_method, n_cores = n_cores)
  } else if (null_model == "spin_brain") {
    message("Generating null brain data with spin_brain model...")
    if (is.null(perm_id)) {
      perm_id <- rotate_parcellation(coord.l = coord.l, coord.r = coord.r, nrot = n_perm, seed = seed)
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
  if (null_model != "none") {
    names(gs.score) <- paste0("null_", 1:n_perm)
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
                      stat2return = c("all", "tval", "pval", "tval_list")) {
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
    {
      if (stat2return == "all") {
        dplyr::mutate(., std_coefs = map(lm_model, ~ standardize_parameters(.x, method = "refit")))
      } else {
        .
      }
    } %>%
    tidyr::unnest(tidy_model) %>%
    dplyr::filter(term == var2extract) %>%
    {
      if (stat2return == "all") {
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
      } else if (stat2return == "pval") {
        dplyr::select(., Dependent_vars, term, p.value) %>%
          dplyr::rename(Predictor = term, pval = p.value) %>%
          dplyr::ungroup()
      }
    }
  return(res)
}





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
#' @param cor_method Character string specifying the correlation method. Default is 'pearson'.
#'                   Other options include 'spearman', 'pls1c', 'pls1w', 'custom'.
#' @param aggre_method Character string specifying the aggregation method. Default is 'mean'.
#'                     Other options include 'median', 'meanabs', 'meansqr', 'maxmean', 'ks_orig', 'ks_weighted',
#'                     'ks_pos_neg_sum', 'sign_test', 'rank_sum', 'custom'.
#' @param null_model Character string specifying the null model method. Default is 'spin_brain'.
#'                   Other options include 'resample_gene', 'coexp_matched'.
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
#' @param pAdjustMethod Character string specifying the method for p-value adjustment. Default is 'fdr'.
#' @param matchcoexp_tol Numeric value specifying the tolerance for matched coexpression. Default is 0.05.
#' @param matchcoexp_max_iter Integer specifying the maximum number of iterations for matched coexpression. Default is 1000000.
#' @return A data frame containing the results of the linear model test, including p-values, adjusted p-values,
#'         q-values, descriptions, and core genes.
#' @export
brainscore.lm_test <- function(pred_df,
                               cov_df,
                               brain_data,
                               gene_data,
                               annoData,
                               cor_method = c("pearson", "spearman", "pls1c", "pls1w", "custom"),
                               aggre_method = c(
                                 "mean", "median", "meanabs", "meansqr", "maxmean",
                                 "ks_orig", "ks_weighted", "ks_pos_neg_sum", "sign_test", "rank_sum", "custom"
                               ),
                               null_model = c("spin_brain", "resample_gene", "coexp_matched"), # score in null setting?
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
                               pAdjustMethod = "fdr",
                               matchcoexp_tol = 0.05,
                               matchcoexp_max_iter = 1000000) {
  # Validate arguments
  cor_method <- match.arg(cor_method)
  aggre_method <- match.arg(aggre_method)
  null_model <- match.arg(null_model)
  threshold_type <- match.arg(threshold_type)
  
  message("=========Emprical model======")
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
    n_cores = n_cores)
  dependent_df.true <- data.frame(gsScore.true, check.names = FALSE)
  res <- simple_lm(dependent_df = dependent_df.true, pred_df = pred_df, cov_df = cov_df, stat2return = "all")
  stat.true <- simple_lm(dependent_df = dependent_df.true, pred_df = pred_df, cov_df = cov_df, stat2return = "tval_list")

  message("=========Null model======")
  # Generate null gene set scores
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

  stat.tmp <- list()
  for (i in 1:length(gsScoreList.null)) {
    dependent_df.null <- data.frame(gsScoreList.null[[i]], check.names = FALSE)
    stat.tmp[[i]] <- simple_lm(
      dependent_df = dependent_df.null,
      pred_df = pred_df,
      cov_df = cov_df, stat2return = "tval_list"
    )
  }
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
  res <- res[res$np.padj <= pvalueCutoff, ]
  res <- res[order(res$np.pval), ]

  if (nrow(res) == 0) {
    message("None of the gene sets are significant at the given p-value cutoff.")
  } else {
    if (threshold_type != "none") {
      message("Identifying core genes...")
      geneList <- corr_brain_gene(brain_data = brain_data, gene_data = gene_data, method = cor_method)
      geneSetList <- get_geneSetList(annoData)
      selected.gs <- filter_geneSetList(rownames(geneList), geneSetList, minGSSize = minGSSize, maxGSSize = maxGSSize)

      survived.gs <- selected.gs[res$Dependent_vars]
      core_genes <- find_core_genes(geneList, survived.gs, pred_df = pred_df, cov_df = cov_df, aggre_method = aggre_method, n_cores = n_cores, threshold_type = threshold_type, threshold_value = threshold_value)
      res$core_genes <- sapply(core_genes, paste0, collapse = "/")
    }
    res <- res %>% rename(ID=Dependent_vars) %>%
    dplyr::select(ID, Description, everything())
    message("Analysis complete.")
  }
  return(res)
}

#' Perform Brain Score Simulation
#'
#' This function performs simulations on brain score data using different methods for comparison.
#' It calculates gene set scores, performs linear modeling, and returns the simulation results.
#'
#' @param pred_df Data frame of predictor variables.
#' @param cov_df Data frame of covariate variables.
#' @param brain_data Data frame of brain imaging data.
#' @param gene_data Data frame of gene expression data.
#' @param annoData Environment containing annotation data.
#' @param sim_n Integer specifying the number of simulations. Default is 1000.
#' @param subsample_size Integer or vector specifying the subsample sizes. Default is 100.
#' @param sim_type Character string specifying the simulation type. Default is 'randomize_pred'.
#'                 Other options include 'spin_brain', 'resample_gene', 'coexp_matched'.
#' @param cor_method Character string specifying the correlation method. Default is 'pearson'.
#'                   Other options include 'spearman', 'pls1c', 'pls1w', 'custom'.
#' @param aggre_method Character string specifying the aggregation method. Default is 'mean'.
#'                     Other options include 'median', 'meanabs', 'meansqr', 'maxmean',
#'                     'ks_orig', 'ks_weighted', 'ks_pos_neg_sum', 'sign_test', 'rank_sum', 'custom'.
#' @param null_model Character string specifying the null model method. Default is 'spin_brain'.
#'                   Other option is 'resample_gene'.
#' @param minGSSize Integer specifying the minimum gene set size. Default is 10.
#' @param maxGSSize Integer specifying the maximum gene set size. Default is 200.
#' @param n_cores Integer specifying the number of cores to use for parallel processing. Default is 0.
#' @param n_perm Integer specifying the number of permutations. Default is 5000.
#' @param perm_id Optional permutation ID.
#' @return A list of data frames containing the results of the simulations.
#' @export
brainscore.simulate <- function(pred_df,
                                cov_df,
                                brain_data,
                                gene_data,
                                annoData,
                                sim_n = 1000,
                                subsample_size = 100,
                                sim_type = c("randomize_pred", "spin_brain", "resample_gene"),
                                cor_method = c("pearson", "spearman", "pls1c", "pls1w", "custom"),
                                aggre_method = c(
                                  "mean", "median", "meanabs", "meansqr", "maxmean",
                                  "ks_orig", "ks_weighted", "ks_pos_neg_sum", "sign_test", "rank_sum", "custom"),
                                minGSSize = 10,
                                maxGSSize = 200,
                                n_cores = 0,
                                n_perm = 5000,
                                perm_id = NULL) {
  # Match arguments
  sim_type <- match.arg(sim_type)
  cor_method <- match.arg(cor_method)
  aggre_method <- match.arg(aggre_method)

  # Initialize results list
  results_list <- list()

  if (sim_type == "randomize_pred") {
    pred_df.sim <- pred_df
    gsScore <- brainscore(
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
    dependent_df <- data.frame(gsScore, check.names = FALSE)

    for (sim_i in 1:sim_n) {
      pred_df.sim[[1]] <- sample(pred_df[[1]])
      for (size_i in seq_along(subsample_size)) {
        size2use <- subsample_size[size_i]
        sampled_idx <- sample(1:nrow(dependent_df), size = size2use, replace = FALSE)
        sampled_res <- simple_lm(
          dependent_df = dependent_df[sampled_idx, , drop = FALSE],
          pred_df = pred_df.sim[sampled_idx, , drop = FALSE],
          cov_df = cov_df[sampled_idx, , drop = FALSE],
          stat2return = "all"
        )
        sampled_res <- sampled_res %>%
          mutate(
            nofdr_ifsig = case_when(
              p.val < 0.05 ~ 1,
              TRUE ~ 0),
            fdr_ifsig = case_when(
              p.adj < 0.05 ~ 1,
              TRUE ~ 0)) %>%
          dplyr::select(Dependent_vars, nofdr_ifsig, fdr_ifsig) %>%
          rename(
            !!paste0("nofdr_sim_", sim_i, "_subsample_", size2use) := nofdr_ifsig,
            !!paste0("fdr_sim_", sim_i, "_subsample_", size2use) := fdr_ifsig)
        results_list[[paste0("sim_", sim_i, "_subsample_", size2use)]] <- sampled_res
      }
    }
  } else if (sim_type == "spin_brain") {
    if (dim(perm_id)[2]< sim_n){
      stop("The number of simulation must be less than or equal to the number of columns in 'perm_id'.")
    }
    for (sim_i in 1:sim_n) {
      sim.brain_data <- brain_data[perm_id[, sim_i], , drop = FALSE]
      rownames(sim.brain_data) <- rownames(brain_data)
      gsScore.true <- brainscore(
        brain_data = sim.brain_data,
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

      gsScoreList.null <- brainscore(
        brain_data = sim.brain_data,
        gene_data = gene_data,
        annoData = annoData,
        cor_method = cor_method,
        aggre_method = aggre_method,
        null_model = "spin_brain",
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        n_cores = n_cores,
        n_perm = n_perm,
        perm_id = perm_id
      )

      for (size_i in seq_along(subsample_size)) {
        size2use <- subsample_size[size_i]
        sampled_idx <- sample(1:nrow(dependent_df.true), size = size2use, replace = FALSE)
        res <- simple_lm(
          dependent_df = dependent_df.true[sampled_idx, , drop = FALSE],
          pred_df = pred_df[sampled_idx, , drop = FALSE],
          cov_df = cov_df[sampled_idx, , drop = FALSE],
          stat2return = "pval"
        )

        stat.true <- simple_lm(
          dependent_df = dependent_df.true[sampled_idx, , drop = FALSE],
          pred_df = pred_df[sampled_idx, , drop = FALSE],
          cov_df = cov_df[sampled_idx, , drop = FALSE],
          stat2return = "tval_list"
        )
        stat.tmp <- list()
        for (i in 1:length(gsScoreList.null)) {
          dependent_df.null <- data.frame(gsScoreList.null[[i]], check.names = FALSE)
          stat.tmp[[i]] <- simple_lm(
            dependent_df = dependent_df.null[sampled_idx, , drop = FALSE],
            pred_df = pred_df[sampled_idx, , drop = FALSE],
            cov_df = cov_df[sampled_idx, , drop = FALSE],
            stat2return = "tval_list"
          )
        }
        stat.null <- list_transpose(stat.tmp)
        np_pval <- calculate_pvals(stat.true, stat.null, method = "standard")
        np_p.adj <- p.adjust(np_pval, method = 'fdr')

        check_names <- all(
          names(stat.true) == names(stat.null),
          names(stat.null) == names(np_pval),
          names(np_pval) == names(np_p.adj),
          names(np_p.adj) == res$Dependent_vars
        )
        if (!check_names) {
          stop("The names of the results are not consistent.")
        } else {
          res$p.adj <- p.adjust(res$pval, method = "fdr")
          res$np_pval <- np_pval
          res$np_p.adj <- np_p.adj
        }
        sampled_res <- res %>%
          mutate(
            nofdr_ifsig = case_when(
              pval < 0.05 ~ 1,
              TRUE ~ 0
            ),
            fdr_ifsig = case_when(
              p.adj < 0.05 ~ 1,
              TRUE ~ 0
            ),
            np_nofdr_ifsig = case_when(
              np_pval < 0.05 ~ 1,
              TRUE ~ 0
            ),
            np_fdr_ifsig = case_when(
              np_p.adj < 0.05 ~ 1,
              TRUE ~ 0
            )
          ) %>%
          dplyr::select(Dependent_vars, nofdr_ifsig, fdr_ifsig, np_nofdr_ifsig, np_fdr_ifsig) %>%
          rename(
            !!paste0("nofdr_sim_", sim_i, "_subsample_", size2use) := nofdr_ifsig,
            !!paste0("fdr_sim_", sim_i, "_subsample_", size2use) := fdr_ifsig,
            !!paste0("np_nofdr_sim_", sim_i, "_subsample_", size2use) := np_nofdr_ifsig,
            !!paste0("np_fdr_sim_", sim_i, "_subsample_", size2use) := np_fdr_ifsig
          )
        results_list[[paste0("sim_", sim_i, "_subsample_", size2use)]] <- sampled_res
      }
    }
  } else if (sim_type == "resample_gene") {
    geneList <- corr_brain_gene(gene_data = gene_data, brain_data = brain_data, method = cor_method)
    geneSetList <- get_geneSetList(annoData)
    selected.gs <- filter_geneSetList(rownames(geneList), geneSetList, minGSSize = minGSSize, maxGSSize = maxGSSize)

    for (sim_i in 1:sim_n) {
      sim.geneList <- geneList[sample(1:nrow(geneList), size = nrow(geneList), replace = FALSE), ]
      rownames(sim.geneList) <- rownames(geneList)
      sim.gsScore <- aggregate_geneSetList(sim.geneList, selected.gs, method = aggre_method, n_cores = n_cores)
      dependent_df.true <- data.frame(sim.gsScore, check.names = FALSE)
       message("Aggregating gene set scores in resample_gene mode...")
      progress_interval <- max(1, round(n_perm / 10))     
      gsScoreList.null <- lapply(1:n_perm, function(idx) {
         if (idx %% progress_interval == 0) {
        message(paste("Processing permutation", idx, "of", n_perm, "..."))
      }
        geneList.null <- geneList[sample(1:nrow(geneList), size = nrow(geneList), replace = FALSE), ]
        rownames(geneList.null) <- rownames(geneList)
        gs_score.null <- aggregate_geneSetList(geneList.null, selected.gs, method = aggre_method, n_cores = n_cores)
        return(gs_score.null)
      })

      for (size_i in seq_along(subsample_size)) {
        size2use <- subsample_size[size_i]
        sampled_idx <- sample(1:nrow(dependent_df.true), size = size2use, replace = FALSE)
        res <- simple_lm(
          dependent_df = dependent_df.true[sampled_idx, , drop = FALSE],
          pred_df = pred_df[sampled_idx, , drop = FALSE],
          cov_df = cov_df[sampled_idx, , drop = FALSE],
          stat2return = "pval"
        )
        stat.true <- simple_lm(
          dependent_df = dependent_df.true[sampled_idx, , drop = FALSE],
          pred_df = pred_df[sampled_idx, , drop = FALSE],
          cov_df = cov_df[sampled_idx, , drop = FALSE],
          stat2return = "tval_list"
        )
        stat.tmp <- list()
        for (i in seq_along(gsScoreList.null)) {
          dependent_df.null <- data.frame(gsScoreList.null[[i]], check.names = FALSE)
          stat.tmp[[i]] <- simple_lm(
            dependent_df = dependent_df.null[sampled_idx, , drop = FALSE],
            pred_df = pred_df[sampled_idx, , drop = FALSE],
            cov_df = cov_df[sampled_idx, , drop = FALSE],
            stat2return = "tval_list"
          )
        }
        stat.null <- list_transpose(stat.tmp)
        np_pval <- calculate_pvals(stat.true, stat.null, method = "standard")
        np_p.adj <- p.adjust(np_pval, method = 'fdr')

        check_names <- all(
          names(stat.true) == names(stat.null),
          names(stat.null) == names(np_pval),
          names(np_pval) == names(np_p.adj),
          names(np_p.adj) == res$Dependent_vars
        )
        if (!check_names) {
          stop("The names of the results are not consistent.")
        } else {
          res$p.adj <- p.adjust(res$pval, method = "fdr")
          res$np_pval <- np_pval
          res$np_p.adj <- np_p.adj
        }
        sampled_res <- res %>%
          mutate(nofdr_ifsig = case_when(
                        pval < 0.05 ~ 1,
                        TRUE ~ 0),
            fdr_ifsig = case_when(
                        p.adj < 0.05 ~ 1,
                        TRUE ~ 0),
            np_nofdr_ifsig = case_when(
                        np_pval < 0.05 ~ 1,
                        TRUE ~ 0),
            np_fdr_ifsig = case_when(
                      np_p.adj < 0.05 ~ 1,
                      TRUE ~ 0)) %>%
          dplyr::select(Dependent_vars, nofdr_ifsig, fdr_ifsig, np_nofdr_ifsig, np_fdr_ifsig) %>%
          rename(
            !!paste0("nofdr_sim_", sim_i, "_subsample_", size2use) := nofdr_ifsig,
            !!paste0("fdr_sim_", sim_i, "_subsample_", size2use) := fdr_ifsig,
            !!paste0("np_nofdr_sim_", sim_i, "_subsample_", size2use) := np_nofdr_ifsig,
            !!paste0("np_fdr_sim_", sim_i, "_subsample_", size2use) := np_fdr_ifsig
          )
        results_list[[paste0("sim_", sim_i, "_subsample_", size2use)]] <- sampled_res
      }
    }
  }

  return(results_list)
}
