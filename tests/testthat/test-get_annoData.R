library(testthat)

test_that("get_annoData returns EnrichAnno for SynGO", {
  anno <- get_annoData("SynGO")
  expect_s3_class(anno, "EnrichAnno")
  expect_true(all(c("term_id", "gene_id") %in% names(anno$term2gene)))
  expect_true(all(c("term_id", "term_name") %in% names(anno$term2name)))
  expect_gt(nrow(anno$term2gene), 0)
})
