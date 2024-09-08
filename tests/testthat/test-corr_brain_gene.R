

# Sample data for testing
set.seed(123)
gene_data <- matrix(rnorm(100), nrow = 10)
brain_data <- matrix(rnorm(100), nrow = 10)
rownames(gene_data) <- rownames(brain_data) <- paste0("Sample", 1:10)

test_that("corr_brain_gene works with pearson method", {
  # Test pearson correlation
  res <- corr_brain_gene(gene_data, brain_data, method = "pearson", r2z = FALSE)
  expect_true(is.matrix(res))
  expect_equal(ncol(res), ncol(brain_data))  # Ensure correct dimensions
  expect_equal(nrow(res), ncol(gene_data))   # Ensure correct dimensions
  expect_equal(attr(res, "cor_type"), "pearson")  # Check attribute
  expect_null(attr(res, "is_fisherz"))  # Fisher Z should be NULL
})


test_that("corr_brain_gene throws error when rownames do not match", {
  # Create gene_data with mismatched rownames
  brain_data_mismatch <- matrix(rnorm(100), nrow = 10)
  rownames(brain_data_mismatch) <- paste0("DifferentSample", 1:10)

  expect_error(corr_brain_gene(gene_data, brain_data_mismatch, method = "pearson"), 
               "Row names of brain_data and gene_data must match.")
})

test_that("corr_brain_gene works with custom method", {
  # Define a custom association function
  custom_func <- function(gene_data, brain_data) {
    return(cor(gene_data, brain_data, method = "spearman"))
  }

  # Test custom method
  res <- corr_brain_gene(gene_data, brain_data, method = custom_func)
  expect_true(is.matrix(res))
  expect_equal(attr(res, "cor_type"), "custom")
})
