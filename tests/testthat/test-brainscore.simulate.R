data(sample_df)
brain_data <- dplyr::select(sample_df, starts_with("L_")) %>% t()
colnames(brain_data) <- paste0("sub-", 1:ncol(brain_data))

gene_data <- get_geneExp(atlas = "desikan", rdonor = "r0.6", hem = "L")
annoData <- get_annoData(type = "SynGO")

cov_df <- sample_df %>% dplyr::select(Age, Sex)
pred_df <- sample_df %>% dplyr::select(BMI)


# Basic test of brainscore.simulate with randomize_pred
test_that("brainscore.simulate works with randomize_pred", {
  # Run the simulation
  res <- brainscore.simulate(
    pred_df = pred_df,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    sim_n = 5,
    subsample_size = 100,
    sim_type = "randomize_pred",
    cor_method = "pearson",
    aggre_method = "mean",
    minGSSize = 20,
    maxGSSize = 200
  )

  # Test that the result is a list
  expect_type(res, "list")

  # Check that the output structure is valid
  expect_length(res, 5)
})

# Test with spin_brain simulation type
test_that("brainscore.simulate works with spin_brain", {
  data(perm_id_dk_lh_5000)
  # Run the simulation
  res <- brainscore.simulate(
    pred_df = pred_df,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    sim_n = 5,
    subsample_size = 100,
    sim_type = "spin_brain",
    cor_method = "pearson",
    aggre_method = "mean",
    n_perm = 10,
    perm_id = perm_id_dk_lh_5000,
    minGSSize = 20,
    maxGSSize = 200
  )

  # Test that the result is a list
  expect_type(res, "list")

  # Check the result length is equal to number of simulations
  expect_length(res, 5)
})

# Test resample_gene simulation type
test_that("brainscore.simulate works with resample_gene", {
  # Run the simulation
  res <- brainscore.simulate(
    pred_df = pred_df,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    sim_n = 5,
    subsample_size = 100,
    sim_type = "resample_gene",
    cor_method = "pearson",
    aggre_method = "mean",
    n_perm = 10,
    minGSSize = 20,
    maxGSSize = 200
  )

  # Test that the result is a list
  expect_type(res, "list")
})



# Basic test of brainscore.simulate with randomize_pred
test_that("brainscore.simulate works with randomize_pred & power setting", {
  # Run the simulation
  res <- brainscore.simulate(
    pred_df = pred_df,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    sim_n = 5,
    subsample_size = 100,
    sim_setting="power",
    sim_type = "randomize_pred",
    cor_method = "pearson",
    aggre_method = "mean",
    minGSSize = 20,
    maxGSSize = 200
  )

  # Test that the result is a list
  expect_type(res, "list")

  # Check that the output structure is valid
  expect_length(res, 5)
})

# Test with spin_brain simulation type
test_that("brainscore.simulate works with spin_brain & power", {
  data(perm_id_dk_lh_5000)
  # Run the simulation
  res <- brainscore.simulate(
    pred_df = pred_df,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    sim_n = 5,
    subsample_size = 100,
    sim_setting="power",
    sim_type = "spin_brain",
    cor_method = "pearson",
    aggre_method = "mean",
    n_perm = 10,
    perm_id = perm_id_dk_lh_5000,
    minGSSize = 20,
    maxGSSize = 200
  )

  # Test that the result is a list
  expect_type(res, "list")

  # Check the result length is equal to number of simulations
  expect_length(res, 5)
})

# Test resample_gene simulation type
test_that("brainscore.simulate works with resample_gene & power", {
  # Run the simulation
  res <- brainscore.simulate(
    pred_df = pred_df,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    sim_n = 5,
    subsample_size = 100,
     sim_setting="power",
    sim_type = "resample_gene",
    cor_method = "pearson",
    aggre_method = "mean",
    n_perm = 10,
    minGSSize = 20,
    maxGSSize = 200
  )

  # Test that the result is a list
  expect_type(res, "list")
})









# Test with pre-calculated null gsScore
test_that("brainscore.simulate works with pre-calculated null gsScore", {
  # Pre-calculate the null gsScore
  gsScore.null <- brainscore(
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    null_model = "spin_brain",
    n_perm = 10,
    perm_id = perm_id_dk_lh_5000,
    minGSSize = 20,
    maxGSSize = 200
  )

  # Run the simulation using the pre-calculated gsScore
  res <- brainscore.simulate(
    pred_df = pred_df,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    gsScoreList.null = gsScore.null,
    sim_n = 5,
    subsample_size = 100,
    sim_type = "spin_brain",
    cor_method = "pearson",
    aggre_method = "mean",
    n_perm = 10,
    perm_id = perm_id_dk_lh_5000,
    minGSSize = 20,
    maxGSSize = 200
  )

  # Test that the result is a list
  expect_type(res, "list")
})


# Test with pre-calculated null gsScore
test_that("brainscore.simulate works with pre-calculated null gsScore", {
  # Pre-calculate the null gsScore
  gsScore.null <- brainscore(
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    null_model = "resample_gene",
    n_perm = 10,
    perm_id = perm_id_dk_lh_5000,
    minGSSize = 20,
    maxGSSize = 200
  )

  # Run the simulation using the pre-calculated gsScore

  expect_error(res <- brainscore.simulate(
    pred_df = pred_df,
    cov_df = cov_df,
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    gsScoreList.null = gsScore.null,
    sim_n = 5,
    subsample_size = 100,
    sim_type = "spin_brain",
    cor_method = "pearson",
    aggre_method = "mean",
    n_perm = 10,
    perm_id = perm_id_dk_lh_5000,
    minGSSize = 20,
    maxGSSize = 200
  ), "Please review the mismatches above.")
})
