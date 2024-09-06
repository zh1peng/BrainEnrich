library(testthat)
library(utils)

# Retrieve real annotation data
annoData <- get_annoData() # Use a valid annotation type

# Test the split_Anno function
test_that("split_Anno handles valid and invalid inputs", {
  # Test with valid environment
  result <- split_Anno(annoData)

  # Test that path2gene and path2name are data frames
  expect_type(result$path2gene, "list") # A data frame is technically a list in R
  expect_true(is.data.frame(result$path2gene))

  # If PATHID2NAME is available, check it as well
  if (!is.null(result$path2name)) {
    expect_true(is.data.frame(result$path2name))
  } else {
    expect_null(result$path2name)
  }

  # Test with an empty environment
  empty_env <- new.env()

  # This should also raise an error
  expect_error(split_Anno(empty_env), "The environment does not contain the required objects.")
})
