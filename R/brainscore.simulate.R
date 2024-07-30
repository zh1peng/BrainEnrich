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
                                  "ks_orig", "ks_weighted", "ks_pos_neg_sum", "sign_test", "rank_sum", "custom"
                                ),
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
      message(paste("=========Processing simulation:", sim_i, "/", sim_n, "========="))
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
              TRUE ~ 0
            ),
            fdr_ifsig = case_when(
              p.adj < 0.05 ~ 1,
              TRUE ~ 0
            )
          ) %>%
          dplyr::select(Dependent_vars, nofdr_ifsig, fdr_ifsig) %>%
          rename(
            !!paste0("nofdr_sim_", sim_i, "_subsample_", size2use) := nofdr_ifsig,
            !!paste0("fdr_sim_", sim_i, "_subsample_", size2use) := fdr_ifsig
          )
        results_list[[paste0("sim_", sim_i, "_subsample_", size2use)]] <- sampled_res
      }
    }
  } else if (sim_type == "spin_brain") {
    if (dim(perm_id)[2] < sim_n) {
      stop("The number of simulation must be less than or equal to the number of columns in 'perm_id'.")
    }
    for (sim_i in 1:sim_n) {
      message(paste("=========Processing simulation:", sim_i, "/", sim_n, "========="))
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
        np_p.adj <- p.adjust(np_pval, method = "fdr")

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
      message(paste("=========Processing simulation:", sim_i, "/", sim_n, "========="))
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
        np_p.adj <- p.adjust(np_pval, method = "fdr")

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
  }

  return(results_list)
}
