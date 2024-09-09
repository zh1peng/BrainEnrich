library(testthat)


# Test the get_annoData function using a real temporary directory
test_that("get_annoData works with actual paths in tempdir", {
  # Use a temporary directory for the test
  temp_dir <- tempdir()

  # Create a temporary package directory structure similar to the actual package
  package_dir <- file.path(temp_dir, "BrainEnrich")
  extdata_dir <- file.path(package_dir, "extdata")
  gene_set_dir <- file.path(extdata_dir, "geneSets")


  if (dir.exists(gene_set_dir)) {
    unlink(gene_set_dir, recursive = TRUE)
  }

  # Simulate the package being installed at the tempdir path
  dir.create(gene_set_dir, recursive = TRUE)

  # Stub the find.package() function to return the temporary package directory
  stub(get_annoData, "find.package", function(x) package_dir)

  # Test scenario where the file doesn't exist and needs to be downloaded
  file_path <- file.path(gene_set_dir, "CellTypes_Lake2018.rds")

  # Remove the file if it exists (to simulate downloading)
  if (file.exists(file_path)) {
    file.remove(file_path)
  }

  # Run the function, expect it to attempt the download
  expect_message(get_annoData("CellTypes_Lake2018"), "File not found locally. Downloading from GitHub...")

  # Check that the file was "downloaded" (since no actual download, we check the file path)
  expect_true(file.exists(file_path))

  # Test scenario where the file is corrupted (or unreadable)
  # Corrupt the file by writing some bad data to it
  writeLines("corrupted content", file_path)

  # Expect an error when trying to read the corrupted file
  expect_error(get_annoData("CellTypes_Lake2018"), "Failed to read the RDS file")

  # Verify that the corrupted file was removed
  expect_false(file.exists(file_path))

  unlink(gene_set_dir, recursive = TRUE)
})
