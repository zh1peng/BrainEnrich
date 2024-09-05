
library(testthat)

# Example test for filter_geneSetList()
annoData <- get_annoData()
geneSetList <- get_geneSetList(annoData)
gene_data <- get_geneExp("desikan", "r0.6", "L")
bg_genes <- colnames(gene_data)

test_that("filter_geneSetList filters gene sets correctly", {
  
  filtered_geneSetList <- filter_geneSetList(bg_genes, geneSetList,1,200)
  expect_true(all(sapply(filtered_geneSetList, is.vector) | sapply(filtered_geneSetList, is.data.frame)))
  expect_true(all(sapply(filtered_geneSetList, function(x) all(x %in% bg_genes))))
})
