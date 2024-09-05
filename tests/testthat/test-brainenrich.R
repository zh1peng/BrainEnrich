library(testthat)
library(mockery)

# Prepare mock data for testing

# Example brain data (small data frame with region names as rownames)
data(brain_data)
# Example gene data (matching rownames with brain_data)
gene_data <- get_geneExp(atlas = "desikan", rdonor = "r0.6", hem = "L")
# Example annotation data (using get_annoData function)
annoData <- get_annoData("SynGO")
geneSetList <- get_geneSetList(annoData)
geneSetList.filtered <- filter_geneSetList(colnames(gene_data), geneSetList, minGSSize = 20, maxGSSize = 200)

# Test the brainenrich function
test_that("brainenrich performs gene set analysis correctly with valid input (spin brain + perm id)", {
data(perm_id_dk_lh_5000)
  # Perform analysis with valid inputs
res<-brainenrich(
  brain_data = brain_data,
  gene_data = gene_data,
  annoData = annoData,
  perm_id = perm_id_dk_lh_5000,
  cor_method = "pearson",
  aggre_method = "mean",
  null_model = "spin_brain",
  n_perm = 2,
  n_cores = 2,
  minGSSize = 20,
  maxGSSize = 200,
  pvalueCutoff = 50,
  threshold_type = "sd",
  threshold = 1
)
  # Test that the result is a gseaResult object
  expect_s4_class(res, "gseaResult")
})


# Test the brainenrich function
test_that("brainenrich performs gene set analysis correctly with valid input (spin brain + coord.l)", {
data(coord_dk_lh)
  # Perform analysis with valid inputs
res<-brainenrich(
  brain_data = brain_data,
  gene_data = gene_data,
  annoData = annoData,
  coord.l = coord_dk_lh,
  perm_id = NULL,
  seed = 2024,
  cor_method = "pearson",
  aggre_method = "mean",
  null_model = "spin_brain",
  n_perm = 2,
  n_cores = 2,
  minGSSize = 20,
  maxGSSize = 200,
  pvalueCutoff = 50,
  threshold_type = "sd",
  threshold = 1
)
  # Test that the result is a gseaResult object
  expect_s4_class(res, "gseaResult")
})

# Test the brainenrich function
test_that("brainenrich performs gene set analysis correctly with valid input (resample gene)", {
  # Perform analysis with valid inputs
res<-brainenrich(
  brain_data = brain_data,
  gene_data = gene_data,
  annoData = annoData,
  cor_method = "pearson",
  aggre_method = "mean",
  null_model = "resample_gene",
  n_perm = 2,
  n_cores = 0,
  minGSSize = 20,
  maxGSSize = 200,
  pvalueCutoff = 1,
  threshold_type = "percentile",
  threshold = 80
)
  # Test that the result is a gseaResult object
  expect_s4_class(res, "gseaResult")
})



# Test the brainenrich function
test_that("brainenrich performs gene set analysis correctly with valid input (co-exp matched)", {

mock_readline <- function(prompt) {
    return("Y")  # Simulate user input as "Y"
  }

with_mock(
  readline = mock_readline,
  # Perform analysis with valid inputs
res<-brainenrich(
  brain_data = brain_data,
  gene_data = gene_data,
  annoData = annoData,
  cor_method = "pearson",
  aggre_method = "mean",
  null_model = "coexp_matched",
  n_perm = 2,
  n_cores = 2,
  minGSSize = 20,
  maxGSSize = 200,
  pvalueCutoff = 1,
  threshold_type = "none",
  matchcoexp_tol = 0.8,
  matchcoexp_max_iter = 1000,
)
)
  # Test that the result is a gseaResult object
  expect_s4_class(res, "gseaResult")
})
sampled_gs <- resample_geneSetList_matching_coexp(gene_data, 
geneSetList.filtered, 
tol = 0.8, 
max_iter = 1000, 
n_perm = 2, 
n_cores = 0)
seq_along(geneSetList.filtered)



# Define valid cor_method and aggre_method options
valid_cor_methods <- c("spearman", "pls1c", "pls1w")
valid_aggre_methods <- c("median", "meanabs", "meansqr", "maxmean", 
                         "ks_orig", "ks_weighted", "ks_pos_neg_sum", "sign_test", 
                         "rank_sum")

# Loop through all combinations of cor_method and aggre_method
for (cor_method in valid_cor_methods) {
  for (aggre_method in valid_aggre_methods) {
    test_that(paste0("brainenrich works with cor_method = '", cor_method, 
                     "' and aggre_method = '", aggre_method, "'"), {
      # Perform analysis with the current combination of cor_method and aggre_method
      res <- brainenrich(
        brain_data = brain_data,
        gene_data = gene_data,
        annoData = annoData,
        cor_method = cor_method,
        aggre_method = aggre_method,
        null_model = "resample_gene",
        n_perm = 2,
        n_cores = 1,
        minGSSize = 20,
        maxGSSize = 200,
        pvalueCutoff = 1,
        threshold_type = "none"
      )

      # Test that the result is a gseaResult object
      expect_s4_class(res, "gseaResult")
    })
  }
}
