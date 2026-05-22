make_normality_brainscore_fixture <- function() {
  set.seed(10)
  regions <- paste0("R", seq_len(8))
  genes <- paste0("G", seq_len(20))
  subjects <- paste0("sub-", seq_len(24))

  brain_data <- matrix(rnorm(length(regions) * length(subjects)),
    nrow = length(regions),
    dimnames = list(regions, subjects)
  )
  gene_data <- matrix(rnorm(length(regions) * length(genes)),
    nrow = length(regions),
    dimnames = list(regions, genes)
  )
  build_Anno <- getFromNamespace("build_Anno", "DOSE")
  annoData <- build_Anno(
    data.frame(
      term = rep(c("term_a", "term_b"), c(8, 7)),
      gene = c(genes[1:8], genes[9:15])
    ),
    data.frame(
      term = c("term_a", "term_b"),
      description = c("Term A", "Term B")
    )
  )
  cov_df <- data.frame(Age = rnorm(length(subjects)))
  pred_df <- data.frame(BMI = rnorm(length(subjects)))

  list(
    brain_data = brain_data,
    gene_data = gene_data,
    annoData = annoData,
    cov_df = cov_df,
    pred_df = pred_df
  )
}

test_that("brainscore.normality reports KS diagnostics for each term", {
  set.seed(1)
  gsScore <- data.frame(
    normal_term = rnorm(200),
    skewed_term = rexp(200),
    check.names = FALSE
  )

  res <- brainscore.normality(gsScore, method = "ks", p_adjust_method = "fdr")

  expect_s3_class(res, "data.frame")
  expect_equal(res$ID, colnames(gsScore))
  expect_true(all(c(
    "normality_method",
    "normality_n",
    "normality_ks_D",
    "normality_ks_p",
    "normality_ks_padj",
    "normality_flag",
    "normality_reason"
  ) %in% colnames(res)))
  expect_false(any(c(
    "normality_mean",
    "normality_sd",
    "normality_median",
    "normality_min",
    "normality_max",
    "normality_skewness",
    "normality_excess_kurtosis"
  ) %in% colnames(res)))
  expect_true(res$normality_flag[res$ID == "skewed_term"])
})

test_that("brainscore.normality uses adjusted p-values for flags", {
  set.seed(2)
  gsScore <- data.frame(
    term_a = rexp(100),
    term_b = rnorm(100),
    check.names = FALSE
  )

  res <- brainscore.normality(gsScore, method = "ks", p_adjust_method = "none")

  expect_equal(res$normality_ks_padj, res$normality_ks_p)
  expect_equal(res$normality_flag, res$normality_ks_padj < 0.05)
})

test_that("brainscore.normality handles invalid input columns without failing", {
  gsScore <- data.frame(
    too_small = c(1, 2, NA, NA),
    constant = rep(1, 4),
    all_na = rep(NA_real_, 4),
    check.names = FALSE
  )

  res <- brainscore.normality(gsScore)

  expect_true(all(is.na(res$normality_ks_D)))
  expect_true(all(is.na(res$normality_ks_p)))
  expect_true(all(is.na(res$normality_ks_padj)))
  expect_true(all(is.na(res$normality_flag)))
  expect_true(all(!is.na(res$normality_reason)))
})

test_that("brainscore.normality can include Shapiro diagnostics", {
  set.seed(3)
  gsScore <- data.frame(term_a = rnorm(50), check.names = FALSE)

  res <- brainscore.normality(gsScore, method = "both")

  expect_true(all(c(
    "normality_ks_D",
    "normality_ks_p",
    "normality_ks_padj",
    "normality_shapiro_W",
    "normality_shapiro_p",
    "normality_shapiro_padj"
  ) %in% colnames(res)))
  expect_equal(res$normality_method, "both")
})

test_that("brainscore attaches normality diagnostics to empirical scores only", {
  fixture <- make_normality_brainscore_fixture()

  empirical <- suppressMessages(brainscore(
    brain_data = fixture$brain_data,
    gene_data = fixture$gene_data,
    annoData = fixture$annoData,
    null_model = "none",
    minGSSize = 5,
    maxGSSize = 10,
    n_cores = 1,
    normality_check = TRUE
  ))
  expect_s3_class(attr(empirical, "normality"), "data.frame")

  null_scores <- suppressMessages(brainscore(
    brain_data = fixture$brain_data,
    gene_data = fixture$gene_data,
    annoData = fixture$annoData,
    null_model = "resample_gene",
    n_perm = 2,
    minGSSize = 5,
    maxGSSize = 10,
    n_cores = 1,
    normality_check = TRUE
  ))
  expect_null(attr(null_scores, "normality"))
})

test_that("brainscore.lm_test returns normality diagnostics with table output", {
  fixture <- make_normality_brainscore_fixture()

  res <- suppressWarnings(suppressMessages(brainscore.lm_test(
    pred_df = fixture$pred_df,
    cov_df = fixture$cov_df,
    brain_data = fixture$brain_data,
    gene_data = fixture$gene_data,
    annoData = fixture$annoData,
    null_model = "resample_gene",
    n_perm = 2,
    minGSSize = 5,
    maxGSSize = 10,
    n_cores = 1,
    pvalueCutoff = 1,
    gsea_obj = FALSE,
    threshold_type = "none",
    normality_check = TRUE
  )))

  expect_true(all(c(
    "normality_method",
    "normality_ks_D",
    "normality_ks_p",
    "normality_ks_padj",
    "normality_flag"
  ) %in% colnames(res)))
})

test_that("brainscore.lm_test returns normality diagnostics in gseaResult output", {
  fixture <- make_normality_brainscore_fixture()

  res <- suppressWarnings(suppressMessages(brainscore.lm_test(
    pred_df = fixture$pred_df,
    cov_df = fixture$cov_df,
    brain_data = fixture$brain_data,
    gene_data = fixture$gene_data,
    annoData = fixture$annoData,
    null_model = "resample_gene",
    n_perm = 2,
    minGSSize = 5,
    maxGSSize = 10,
    n_cores = 1,
    pvalueCutoff = 1,
    gsea_obj = TRUE,
    normality_check = TRUE
  )))
  res_df <- as.data.frame(res)

  expect_true(all(c(
    "normality_method",
    "normality_ks_D",
    "normality_ks_p",
    "normality_ks_padj",
    "normality_flag"
  ) %in% colnames(res_df)))
})

test_that("brainscore.lm_test rejects unsupported aggregation methods", {
  fixture <- make_normality_brainscore_fixture()

  expect_error(
    brainscore.lm_test(
      pred_df = fixture$pred_df,
      cov_df = fixture$cov_df,
      brain_data = fixture$brain_data,
      gene_data = fixture$gene_data,
      annoData = fixture$annoData,
      aggre_method = "maxmean",
      null_model = "none",
      minGSSize = 5,
      maxGSSize = 10,
      gsea_obj = FALSE
    ),
    "not supported in brainscore.lm_test"
  )
})

test_that("brainscore.lm_test accepts custom aggregation functions", {
  fixture <- make_normality_brainscore_fixture()
  custom_mean <- function(genelist, geneSet) {
    geneSet <- intersect(geneSet, names(genelist))
    mean(genelist[geneSet])
  }

  res <- suppressWarnings(suppressMessages(brainscore.lm_test(
    pred_df = fixture$pred_df,
    cov_df = fixture$cov_df,
    brain_data = fixture$brain_data,
    gene_data = fixture$gene_data,
    annoData = fixture$annoData,
    aggre_method = custom_mean,
    null_model = "none",
    minGSSize = 5,
    maxGSSize = 10,
    n_cores = 1,
    gsea_obj = FALSE
  )))

  expect_s3_class(res, "data.frame")
  expect_true(all(c("term_a", "term_b") %in% res$Dependent_vars))
})
