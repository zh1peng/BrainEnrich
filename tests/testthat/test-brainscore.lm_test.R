data(sim_hcp)
brain_data <- dplyr::select(sim_hcp, starts_with("L_")) %>% t()
colnames(brain_data) <- paste0("sub-", 1:ncol(brain_data))
gene_data <- get_geneExp(atlas = "desikan", rdonor = "r0.6", hem = "L")
annoData <- get_annoData(type = "SynGO")
cov_df <- sim_hcp %>% dplyr::select(Age, Sex)
pred_df <- sim_hcp %>% dplyr::select(BMI)

# Test the brainenrich function
test_that("brainenrich performs gene set analysis correctly with valid input (spin_brain)", {
  data(perm_id_dk_lh_5000)
  # Perform analysis with valid inputs
  res <- brainscore.lm_test(
    pred_df = pred_df,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    cor_method = "pearson",
    aggre_method = "mean",
    n_cores = 0,
    minGSSize = 50,
    maxGSSize = 200,
    null_model = "spin_brain",
    n_perm = 5,
    perm_id = perm_id_dk_lh_5000,
    pvalueCutoff = 1
  )

  # Test that the result is a gseaResult object
  expect_s4_class(res, "gseaResult")
})



# Test the brainenrich function
test_that("brainenrich performs gene set analysis correctly with valid input (spin_brain)", {
  # Perform analysis with valid inputs
  data(coord_dk_lh)
  res <- brainscore.lm_test(
    pred_df = pred_df,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    cor_method = "pearson",
    aggre_method = "mean",
    n_cores = 0,
    minGSSize = 50,
    maxGSSize = 200,
    null_model = "spin_brain",
    n_perm = 5,
    coord.l = coord_dk_lh,
    pvalueCutoff = 1,
    gsea_obj = FALSE,
    threshold_type = "none"
  )

  # Test that the result is a gseaResult object
  expect_type(res, "list")
})



# Test the brainenrich function
test_that("brainenrich performs gene set analysis correctly with valid input (resample_gene)", {
  # Perform analysis with valid inputs
  res <- brainscore.lm_test(
    pred_df = pred_df,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    cor_method = "pearson",
    aggre_method = "mean",
    n_cores = 0,
    minGSSize = 50,
    maxGSSize = 200,
    null_model = "resample_gene",
    n_perm = 5,
    coord.l = coord_dk_lh,
    pvalueCutoff = 1
  )

  # Test that the result is a gseaResult object
  expect_s4_class(res, "gseaResult")
})



pred_df <- sim_hcp %>%
  dplyr::mutate(group = case_when(BMI > 25 ~ "case", TRUE ~ "hc")) %>%
  dplyr::select(group)
pred_df$group <- relevel(factor(pred_df$group), ref = "hc")
# Test the brainenrich function
test_that("brainenrich performs gene set analysis correctly with valid input (resample_gene)", {
  # Perform analysis with valid inputs
  res <- brainscore.lm_test(
    pred_df = pred_df,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    cor_method = "pearson",
    aggre_method = "mean",
    n_cores = 0,
    minGSSize = 50,
    maxGSSize = 200,
    null_model = "resample_gene",
    n_perm = 5,
    pvalueCutoff = 1,
    gsea_obj = FALSE,
    threshold_type = "none"
  )

  # Test that the result is a list
  expect_type(res, "list")
})


pred_df <- sim_hcp %>%
  dplyr::mutate(group = case_when(BMI > 25 ~ 1, TRUE ~ 0)) %>%
  dplyr::select(group)
pred_df$group <- as.factor(pred_df$group)
# Test the brainenrich function
test_that("brainenrich performs gene set analysis correctly with valid input (resample_gene)", {
  # Perform analysis with valid inputs
  res <- brainscore.lm_test(
    pred_df = pred_df,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    cor_method = "pearson",
    aggre_method = "mean",
    n_cores = 0,
    minGSSize = 50,
    maxGSSize = 200,
    null_model = "resample_gene",
    n_perm = 5,
    pvalueCutoff = 1,
    gsea_obj = FALSE,
    threshold_type = "none"
  )

  # Test that the result is a list
  expect_type(res, "list")
})
