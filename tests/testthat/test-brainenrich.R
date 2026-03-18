data(brain_data)
gene_data <- get_geneExp(atlas = "desikan", rdonor = "r0.6", hem = "L")
annoData <- get_annoData("SynGO")

test_that("brainenrich returns EnrichRes with spin_brain perm_id", {
  data(perm_id_dk_lh_5000)
  res <- brainenrich(
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    perm_id = perm_id_dk_lh_5000,
    cor_method = "pearson",
    aggre_method = "mean",
    null_model = "spin_brain",
    n_perm = 10,
    n_cores = 1,
    minGSSize = 20,
    maxGSSize = 200,
    pvalueCutoff = 1,
    threshold_type = "sd",
    threshold_value = 1
  )
  expect_s3_class(res, "EnrichRes")
  expect_equal(res$analysis_type, "brainenrich")
})

test_that("brainenrich returns EnrichRes with resample_gene", {
  res <- brainenrich(
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    cor_method = "pearson",
    aggre_method = "mean",
    null_model = "resample_gene",
    n_perm = 10,
    n_cores = 1,
    minGSSize = 20,
    maxGSSize = 200,
    pvalueCutoff = 1,
    threshold_type = "none"
  )
  expect_s3_class(res, "EnrichRes")
})

test_that("brainenrich accepts custom methods", {
  custom_cor <- function(gene_data, brain_data) cor(gene_data, brain_data, method = "spearman")
  custom_aggre <- function(genelist, geneSet) mean(genelist[names(genelist) %in% geneSet])
  res <- brainenrich(
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    cor_method = custom_cor,
    aggre_method = custom_aggre,
    null_model = "resample_gene",
    n_perm = 10,
    n_cores = 1,
    minGSSize = 20,
    maxGSSize = 200,
    pvalueCutoff = 1,
    threshold_type = "none"
  )
  expect_s3_class(res, "EnrichRes")
})
