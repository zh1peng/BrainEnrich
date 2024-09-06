data(sample_df)
brain_data <- dplyr::select(sample_df, starts_with("L_")) %>% t()
colnames(brain_data) <- paste0("sub-", 1:ncol(brain_data))

gene_data <- get_geneExp(atlas = "desikan", rdonor = "r0.6", hem = "L")
annoData <- get_annoData(type = "SynGO")

test_that("brainscore with coexp_matched", {
  pkgload::load_all()
  mock_ask_user_continue <- function(msg) {
    return(TRUE)
  }

  local_mocked_bindings(ask_user_continue = mock_ask_user_continue)
  res <- brainscore(
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    null_model = "coexp_matched",
    cor_method = "pearson",
    aggre_method = "mean",
    n_perm = 10,
    minGSSize = 20,
    maxGSSize = 200,
    matchcoexp_tol = 0.8,
    n_cores = 0
  )
  expect_type(res, "list")
  expect_equal(length(res), 10)
})



test_that("brainscore with coexp_matched", {
  pkgload::load_all()
  mock_ask_user_continue <- function(msg) {
    return(FALSE)
  }

  local_mocked_bindings(ask_user_continue = mock_ask_user_continue)
  expect_error(
    res <- brainscore(
      brain_data = brain_data,
      gene_data = gene_data,
      annoData = annoData,
      null_model = "coexp_matched",
      cor_method = "pearson",
      aggre_method = "mean",
      n_perm = 10,
      minGSSize = 20,
      maxGSSize = 200,
      matchcoexp_tol = 0.8,
      n_cores = 0
    ), "Operation cancelled by user."
  )
})
