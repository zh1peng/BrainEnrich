library(testthat)

test_that("get_geneExp loads expression matrix for SynGO workflows", {
  mx <- get_geneExp(atlas = "desikan", rdonor = "r0.6", hem = "L")
  expect_true(is.matrix(mx))
  expect_gt(nrow(mx), 0)
  expect_gt(ncol(mx), 0)
})

test_that("get_geneExp validates atlas argument", {
  expect_error(get_geneExp(atlas = "invalid", rdonor = "r0.6", hem = "L"))
})
