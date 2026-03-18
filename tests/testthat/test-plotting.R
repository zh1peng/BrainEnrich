skip_if_not_installed("ggplot2")

data(brain_data)
data(perm_id_dk_lh_5000)
gene_data <- get_geneExp(atlas = "desikan", rdonor = "r0.6", hem = "L")
annoData <- get_annoData("SynGO")

test_that("native plotting helpers work for brainenrich results", {
  res <- brainenrich(
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    perm_id = perm_id_dk_lh_5000,
    cor_method = "pearson",
    aggre_method = "mean",
    null_model = "spin_brain",
    n_perm = 3,
    n_cores = 1,
    minGSSize = 20,
    maxGSSize = 200,
    pvalueCutoff = 1,
    threshold_type = "sd",
    threshold_value = 1
  )

  expect_s3_class(plot_terms(res, type = "dot", top_n = 5), "ggplot")
  expect_s3_class(plot_terms(res, type = "ridge", top_n = 5), "ggplot")
  expect_s3_class(plot_heatmap_terms(res, top_n = 5, max_perm = 3), "ggplot")
  expect_s3_class(plot_term_network(res, top_n = 5, min_overlap = 0.01), "ggplot")

  top_term <- res$table$ID[[1]]
  expect_s3_class(plot_core_genes(res, term_id = top_term, mode = "impact"), "ggplot")
})

test_that("native plotting helpers work for brainscore.lm_test results", {
  data(sim_hcp)
  sim_hcp_small <- sim_hcp[seq_len(min(80, nrow(sim_hcp))), , drop = FALSE]
  brain_data_ind <- dplyr::select(sim_hcp_small, starts_with("L_")) %>% t()
  colnames(brain_data_ind) <- paste0("sub-", seq_len(ncol(brain_data_ind)))
  pred_df <- dplyr::select(sim_hcp_small, BMI)
  cov_df <- dplyr::select(sim_hcp_small, Age, Sex)

  res <- brainscore.lm_test(
    pred_df = pred_df,
    cov_df = cov_df,
    brain_data = brain_data_ind,
    gene_data = gene_data,
    annoData = annoData,
    cor_method = "pearson",
    aggre_method = "mean",
    n_cores = 1,
    minGSSize = 50,
    maxGSSize = 200,
    null_model = "spin_brain",
    n_perm = 3,
    perm_id = perm_id_dk_lh_5000,
    pvalueCutoff = 1,
    threshold_type = "sd",
    threshold_value = 1
  )

  expect_s3_class(plot_terms(res, type = "volcano"), "ggplot")
  top_term <- res$table$ID[[1]]
  expect_s3_class(plot_core_genes(res, term_id = top_term, mode = "direction"), "ggplot")
})
