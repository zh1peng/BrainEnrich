data(sim_hcp)
sim_hcp_small <- sim_hcp[seq_len(min(80, nrow(sim_hcp))), , drop = FALSE]
brain_data <- dplyr::select(sim_hcp_small, starts_with("L_")) %>% t()
colnames(brain_data) <- paste0("sub-", seq_len(ncol(brain_data)))
gene_data <- get_geneExp(atlas = "desikan", rdonor = "r0.6", hem = "L")
annoData <- get_annoData(type = "SynGO")
cov_df <- sim_hcp_small %>% dplyr::select(Age, Sex)
pred_df_num <- sim_hcp_small %>% dplyr::select(BMI)

test_that("brainscore.lm_test returns EnrichRes in spin_brain mode", {
  data(perm_id_dk_lh_5000)
  res <- brainscore.lm_test(
    pred_df = pred_df_num,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    cor_method = "pearson",
    aggre_method = "mean",
    n_cores = 1,
    minGSSize = 50,
    maxGSSize = 200,
    null_model = "spin_brain",
    n_perm = 10,
    perm_id = perm_id_dk_lh_5000,
    pvalueCutoff = 1
  )
  expect_s3_class(res, "EnrichRes")
  expect_equal(res$analysis_type, "brainscore_lm_test")
})

test_that("brainscore.lm_test returns data frame when gsea_obj = FALSE", {
  data(perm_id_dk_lh_5000)
  res <- brainscore.lm_test(
    pred_df = pred_df_num,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    cor_method = "pearson",
    aggre_method = "mean",
    n_cores = 1,
    minGSSize = 50,
    maxGSSize = 200,
    null_model = "spin_brain",
    n_perm = 10,
    perm_id = perm_id_dk_lh_5000,
    pvalueCutoff = 1,
    gsea_obj = FALSE,
    threshold_type = "none"
  )
  expect_true(is.data.frame(res))
})

test_that("brainscore.lm_test handles factor predictors", {
  pred_df_fac <- sim_hcp_small %>%
    dplyr::mutate(group = dplyr::case_when(BMI > 25 ~ "case", TRUE ~ "hc")) %>%
    dplyr::select(group)
  pred_df_fac$group <- stats::relevel(factor(pred_df_fac$group), ref = "hc")

  data(perm_id_dk_lh_5000)
  res <- brainscore.lm_test(
    pred_df = pred_df_fac,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    cor_method = "pearson",
    aggre_method = "mean",
    n_cores = 1,
    minGSSize = 50,
    maxGSSize = 200,
    null_model = "spin_brain",
    n_perm = 10,
    perm_id = perm_id_dk_lh_5000,
    pvalueCutoff = 1,
    gsea_obj = FALSE,
    threshold_type = "none"
  )
  expect_true(is.data.frame(res))
})
