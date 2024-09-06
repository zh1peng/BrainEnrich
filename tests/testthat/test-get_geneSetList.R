library(testthat)

annoData <- get_annoData() # Example of annotation type
# Example test for get_geneSetList()
test_that("get_geneSetList retrieves gene sets correctly", {
  # Assume that get_annoData is working correctly, so we can retrieve annoData
  # Test the output of get_geneSetList
  geneSetList <- get_geneSetList(annoData)
  # Test that geneSetList is a list
  expect_type(geneSetList, "list")
  # Test that the list is not empty (assuming annoData contains data)
  expect_true(length(geneSetList) > 0)
  # Test that each element in the list has the expected structure (assuming they are vectors or data frames)
  expect_true(all(sapply(geneSetList, is.vector) | sapply(geneSetList, is.data.frame)))
  expect_error(get_geneSetList(NULL), "not supported")
})
