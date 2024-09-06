library(testthat)
library(mockery)

# Test the get_geneExp function
test_that("get_geneExp works with valid and invalid inputs", {
  # Use a temporary directory for the test
  temp_dir <- tempdir()

  # Create a temporary package directory structure
  package_dir <- file.path(temp_dir, "BrainEnrich")
  extdata_dir <- file.path(package_dir, "extdata")
  gene_exp_dir <- file.path(extdata_dir, "geneExp")

  # Ensure the directory is cleaned up before the test
  if (dir.exists(gene_exp_dir)) {
    unlink(gene_exp_dir, recursive = TRUE)
  }

  # Create the necessary directories
  dir.create(gene_exp_dir, recursive = TRUE)

  # Stub the find.package() function to return the temporary package directory
  stub(get_geneExp, "find.package", function(x) package_dir)


  # Test with valid parameters
  atlas <- "schaefer100"
  rdonor <- "r0.6"
  hem <- "L"
  GeneExpCSV <- file.path(gene_exp_dir, sprintf("%s_%s.csv.bz2", atlas, rdonor))

  # Simulate the file not existing and the function downloading it
  if (file.exists(GeneExpCSV)) {
    file.remove(GeneExpCSV)
  }

  expect_message(get_geneExp(atlas, rdonor, hem), "File not found locally. Downloading from GitHub...")
  expect_true(is.matrix(get_geneExp(atlas, rdonor, hem)))

  # Test with valid parameters
  atlas <- "desikan"
  rdonor <- "r0.6"
  hem <- "L"
  GeneExpCSV <- file.path(gene_exp_dir, sprintf("%s_%s.csv.bz2", atlas, rdonor))

  # Simulate the file not existing and the function downloading it
  if (file.exists(GeneExpCSV)) {
    file.remove(GeneExpCSV)
  }

  expect_message(get_geneExp(atlas, rdonor, hem), "File not found locally. Downloading from GitHub...")
  expect_true(is.matrix(get_geneExp(atlas, rdonor, hem)))
  # Simulate a corrupted file and test error handling
  writeLines("corrupted content", GeneExpCSV)

  # Expect an error when trying to read the corrupted file
  expect_error(get_geneExp(atlas, rdonor, hem), regexp = "Failed to read the CSV file.")

  # Verify that the corrupted file was removed
  expect_false(file.exists(GeneExpCSV))



  # Clean up the created files and directories
  unlink(gene_exp_dir, recursive = TRUE)
})
