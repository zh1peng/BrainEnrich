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
#' @param gsScoreList.null Precomputed list of gene set scores for the null model by brainscore/brainscore.hpc function. Default is NULL.
#' @param sim_n Integer specifying the number of simulations. Default is 1000.
#' @param subsample_size Integer or vector specifying the subsample sizes. Default is 100.
#' @param sim_setting Character string specifying the simulation setting. "type1": use shuffled data; "power": use original data. Default is 'type1'.
#' @param sim_type Character string specifying the simulation type. Default is 'randomize_pred'.
#'                 Other options include 'spin_brain', 'resample_gene', 'coexp_matched'.
#' @param cor_method Character string specifying the correlation method. Default is 'pearson'.
#'                   Other options include 'spearman', 'pls1c', 'pls1w', 'custom'.
#' @param aggre_method Character string specifying the aggregation method. Default is 'mean'.
#'                     Other options include 'median', 'meanabs', 'meansqr', 'maxmean',
#'                     'ks_orig', 'ks_weighted', 'ks_pos_neg_sum', 'sign_test', 'rank_sum', 'custom'.
#' @param minGSSize Integer specifying the minimum gene set size. Default is 10.
#' @param maxGSSize Integer specifying the maximum gene set size. Default is 200.
#' @param n_cores Integer specifying the number of cores to use for parallel processing. Default is 0 (use all cores -1).
#' @param n_perm Integer specifying the number of permutations. Default is 5000.
#' @param perm_id Optional permutation ID.
#' @importFrom dplyr case_when select rename %>%
#' @importFrom data.table :=
#' @return A list of data frames containing the results of the simulations.
#' @export
brainscore.simulate <- function(pred_df,
                                cov_df,
                                brain_data,
                                gene_data,
                                annoData,
                                gsScoreList.null = NULL,
                                sim_n = 1000,
                                subsample_size = 100,
                                sim_setting = c("type1", "power"),
                                sim_type = c("randomize_pred", "spin_brain", "resample_gene"),
                                cor_method = c("pearson", "spearman", "pls1c", "pls1w", "custom"),
                                aggre_method = c(
                                  "mean", "median", "meanabs", "meansqr", "maxmean",
                                  "ks_orig", "ks_weighted", "ks_pos_neg_sum", "sign_test", "rank_sum", "custom"
                                ),
                                minGSSize = 10,
                                maxGSSize = 200,
                                n_cores = 0,
                                perm_id = NULL,
                                n_perm = 1000) {
  # Match arguments
  sim_type <- match.arg(sim_type)
  cor_method <- match.arg(cor_method)
  aggre_method <- match.arg(aggre_method)
  sim_setting <- match.arg(sim_setting)

  if (sim_type == "randomize_pred") {
    message("Running randomize_pred simulation.")
    gsScore <- brainscore(
      brain_data = brain_data,
      gene_data = gene_data,
      annoData = annoData,
      cor_method = cor_method,
      aggre_method = aggre_method,
      null_model = "none",
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      n_cores = n_cores,
      verbose = FALSE
    )
    dependent_df <- data.frame(gsScore, check.names = FALSE)

    if (n_cores == 0) {
      n_cores <- max(detectCores() - 1, 1) # Use all cores minus one, but ensure at least 1 core is used
    } else {
      n_cores <- min(n_cores, detectCores()) # Ensure n_cores does not exceed the number of available cores
    }

    cl <- if (n_cores > 1) makeCluster(n_cores) else NULL
    if (!is.null(cl)) {
      clusterExport(cl, c("subsample_size", "pred_df", "dependent_df", "cov_df", "simple_lm", "sim_n", "sim_setting"), envir = environment())
    }
    message("Simulation with randomize_pred model...")
    results_list <- pblapply(1:sim_n, function(sim_i) {
      if (sim_setting == "power") {
        pred_df.sim <- pred_df
      } else if (sim_setting == "type1") {
        pred_df.sim[[1]] <- sample(pred_df[[1]])
      }
      sim_results <- list()
      for (size_i in seq_along(subsample_size)) {
        size2use <- subsample_size[size_i]
        sampled_idx <- sample(1:nrow(dependent_df), size = size2use, replace = FALSE)
        sampled_res <- simple_lm(
          dependent_df = dependent_df[sampled_idx, , drop = FALSE],
          pred_df = pred_df.sim[sampled_idx, , drop = FALSE],
          cov_df = cov_df[sampled_idx, , drop = FALSE],
          stat2return = "pval"
        )
        sampled_res <- sampled_res %>%
          dplyr::select(.data$Dependent_vars, .data$p.val) %>%
          dplyr::rename(
            !!paste0("pval_nofdr_sim_", sim_i, "_subsample_", size2use) := .data$p.val
          )
        sim_results[[paste0("sim_", sim_i, "_subsample_", size2use)]] <- sampled_res
      }
      return(sim_results)
    }, cl = cl)
    results_list <- do.call(c, results_list)
  } else if (sim_type == "spin_brain") {
    message("Running spin_brain simulation.")
    if (is.null(perm_id)) {
      stop("Permutation ID must be provided for spin_brain simulation.")
    }
    if (dim(perm_id)[2] < sim_n) {
      stop("The number of simulation must be less than or equal to the number of columns in 'perm_id'.")
    }

    # To boost efficiency we can put this part outside the sim loop

    if (is.null(gsScoreList.null)) {
      message("gsScoreList.null is not provided. Running brainscore to compute null model. This may take long time!")
      gsScoreList.null <- brainscore(
        brain_data = brain_data,
        gene_data = gene_data,
        annoData = annoData,
        cor_method = cor_method,
        aggre_method = aggre_method,
        null_model = "spin_brain",
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        n_cores = n_cores,
        n_perm = n_perm,
        perm_id = perm_id,
        verbose = FALSE
      )
    } else {
      null_model.precomp <- attr(gsScoreList.null, "null_model")
      cor_method.precomp <- attr(gsScoreList.null, "cor_method")
      aggre_method.precomp <- attr(gsScoreList.null, "aggre_method")
      minGSSize.precomp <- attr(gsScoreList.null, "minGSSize")
      maxGSSize.precomp <- attr(gsScoreList.null, "maxGSSize")
      n_perm.precomp <- attr(gsScoreList.null, "n_perm")

      # Check all attributes at once
      if (!((null_model.precomp == sim_type) &&
        (cor_method.precomp == cor_method) &&
        (aggre_method.precomp == aggre_method) &&
        (minGSSize.precomp == minGSSize) &&
        (maxGSSize.precomp == maxGSSize) &&
        (n_perm.precomp == n_perm))) {
        message("Mismatches found between precomputed attributes and input variables.")
        message("Please check the following variables: null_model, cor_method, aggre_method, minGSSize, maxGSSize, n_perm.")
        stop("Please review the mismatches above.")
      } else {
        message("Using precomputed gsScoreList.null.")
      }
    }

    message("Simulation with spin_brain model...")
    if (n_cores == 0) {
      n_cores <- max(detectCores() - 1, 1) # Use all cores minus one, but ensure at least 1 core is used
    } else {
      n_cores <- min(n_cores, detectCores()) # Ensure n_cores does not exceed the number of available cores
    }

    cl <- if (n_cores > 1) makeCluster(n_cores) else NULL
    if (!is.null(cl)) {
      clusterExport(cl, c(
        "sim_n", "subsample_size", "brain_data", "gene_data", "perm_id", "annoData",
        "cor_method", "aggre_method", "minGSSize", "maxGSSize", "gsScoreList.null",
        "pred_df", "cov_df", "simple_lm", "brainscore", "sim_setting"
      ), envir = environment())
    }
    results_list <- pblapply(1:sim_n, function(sim_i) {
      if (sim_setting == "power") {
        sim.brain_data <- brain_data
      } else if (sim_setting == "type1") {
        sim.brain_data <- brain_data[perm_id[, sim_i], , drop = FALSE]
        rownames(sim.brain_data) <- rownames(brain_data)
      }
      gsScore.true <- brainscore(
        brain_data = sim.brain_data,
        gene_data = gene_data,
        annoData = annoData,
        cor_method = cor_method,
        aggre_method = aggre_method,
        null_model = "none",
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        n_cores = 1, # use 1 core for this as it is already within a parallel loop
        verbose = FALSE
      )

      dependent_df.true <- data.frame(gsScore.true, check.names = FALSE)

      sim_results <- list()

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

        stat.tmp <- lapply(1:length(gsScoreList.null), function(i) {
          dependent_df.null <- data.frame(gsScoreList.null[[i]], check.names = FALSE)
          simple_lm(
            dependent_df = dependent_df.null[sampled_idx, , drop = FALSE],
            pred_df = pred_df[sampled_idx, , drop = FALSE],
            cov_df = cov_df[sampled_idx, , drop = FALSE],
            stat2return = "tval_list"
          )
        })

        stat.null <- list_transpose(stat.tmp)
        np_pval <- calculate_pvals(stat.true, stat.null, method = "standard")

        check_names <- all(
          names(stat.true) == names(stat.null),
          names(stat.null) == names(np_pval),
          names(np_pval) == names(np_p.adj),
        )

        if (!check_names) {
          stop("The names of the results are not consistent.")
        } else {
          res$p.adj <- p.adjust(res$pval, method = "fdr")
          res$np_pval <- np_pval
        }

        sampled_res <- res %>%
          dplyr::select(.data$Dependent_vars, .data$p.val, .data$np_p.pval) %>%
          dplyr::rename(
            !!paste0("pval_nofdr_sim_", sim_i, "_subsample_", size2use) := .data$p.val,
            !!paste0("np_pval_nofdr_sim_", sim_i, "_subsample_", size2use) := .data$np_p.val
          )

        sim_results[[paste0("sim_", sim_i, "_subsample_", size2use)]] <- sampled_res
      }

      return(sim_results)
    }, cl = cl)

    # Flatten the results_list
    results_list <- do.call(c, results_list)
  } else if (sim_type == "resample_gene") {
    geneList <- corr_brain_gene(gene_data = gene_data, brain_data = brain_data, method = cor_method)
    geneSetList <- get_geneSetList(annoData)
    selected.gs <- filter_geneSetList(rownames(geneList), geneSetList, minGSSize = minGSSize, maxGSSize = maxGSSize)

    if (is.null(gsScoreList.null)) {
      message("gsScoreList.null is not provided. Running brainscore to compute null model. This may take long time!")
      gsScoreList.null <- brainscore(
        brain_data = brain_data,
        gene_data = gene_data,
        annoData = annoData,
        cor_method = cor_method,
        aggre_method = aggre_method,
        null_model = "resample_gene",
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        n_cores = n_cores,
        n_perm = n_perm,
        verbose = FALSE
      )
    } else {
      null_model.precomp <- attr(gsScoreList.null, "null_model")
      cor_method.precomp <- attr(gsScoreList.null, "cor_method")
      aggre_method.precomp <- attr(gsScoreList.null, "aggre_method")
      minGSSize.precomp <- attr(gsScoreList.null, "minGSSize")
      maxGSSize.precomp <- attr(gsScoreList.null, "maxGSSize")
      n_perm.precomp <- attr(gsScoreList.null, "n_perm")

      # Check all attributes at once
      if (!((null_model.precomp == sim_type) &&
        (cor_method.precomp == cor_method) &&
        (aggre_method.precomp == aggre_method) &&
        (minGSSize.precomp == minGSSize) &&
        (maxGSSize.precomp == maxGSSize) &&
        (n_perm.precomp == n_perm))) {
        message("Mismatches found between precomputed attributes and input variables.")
        message("Please check the following variables: null_model, cor_method, aggre_method, minGSSize, maxGSSize, n_perm.")
        stop("Please review the mismatches above.")
      } else {
        message("Using precomputed gsScoreList.null.")
      }
    }
    message("Simulation with resample_gene model...")
    if (n_cores == 0) {
      n_cores <- max(detectCores() - 1, 1) # Use all cores minus one, but ensure at least 1 core is used
    } else {
      n_cores <- min(n_cores, detectCores()) # Ensure n_cores does not exceed the number of available cores
    }

    cl <- if (n_cores > 1) makeCluster(n_cores) else NULL
    if (!is.null(cl)) {
      clusterExport(cl, c(
        "sim_n", "subsample_size", "geneList", "selected.gs", "gsScoreList.null",
        "aggre_method", "pred_df", "cov_df", "aggregate_geneSetList", "sim_setting"
      ), envir = environment())
    }

    results_list <- pblapply(1:sim_n, function(sim_i) {
      if (sim_setting == "power") {
        sim.geneList <- geneList
      } else if (sim_setting == "type1") {
        sim.geneList <- geneList[sample(1:nrow(geneList), size = nrow(geneList), replace = FALSE), ]
        rownames(sim.geneList) <- rownames(geneList)
      }
      sim.gsScore <- aggregate_geneSetList(sim.geneList, selected.gs, method = aggre_method, n_cores = 1)
      dependent_df.true <- data.frame(sim.gsScore, check.names = FALSE)

      sim_results <- list()

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
        stat.tmp <- lapply(seq_along(gsScoreList.null), function(i) {
          dependent_df.null <- data.frame(gsScoreList.null[[i]], check.names = FALSE)
          simple_lm(
            dependent_df = dependent_df.null[sampled_idx, , drop = FALSE],
            pred_df = pred_df[sampled_idx, , drop = FALSE],
            cov_df = cov_df[sampled_idx, , drop = FALSE],
            stat2return = "tval_list"
          )
        })
        stat.null <- list_transpose(stat.tmp)
        np_pval <- calculate_pvals(stat.true, stat.null, method = "standard")

        check_names <- all(
          names(stat.true) == names(stat.null),
          names(stat.null) == names(np_pval),
          names(np_pval) == names(np_p.adj)
        )
        if (!check_names) {
          stop("The names of the results are not consistent.")
        } else {
          res$p.adj <- p.adjust(res$pval, method = "fdr")
          res$np_pval <- np_pval
        }
        sampled_res <- res %>%
          dplyr::select(.data$Dependent_vars, .data$p.val, .data$np_p.pval) %>%
          dplyr::rename(
            !!paste0("pval_nofdr_sim_", sim_i, "_subsample_", size2use) := .data$p.val,
            !!paste0("np_pval_nofdr_sim_", sim_i, "_subsample_", size2use) := .data$np_p.val
          )
        sim_results[[paste0("sim_", sim_i, "_subsample_", size2use)]] <- sampled_res
      }
      return(sim_results)
    }, cl = cl)

    # Flatten the results_list
    results_list <- do.call(c, results_list)
  }

  if (!is.null(cl)) stopCluster(cl)
  message("Simulation completed.")
  return(results_list)
}
