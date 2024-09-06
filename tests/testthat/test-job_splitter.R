library(testthat)

test_that("job_splitter works with valid inputs", {
  # Create mock data
  set.seed(123)
  perm_id <- matrix(rnorm(100), nrow = 10, ncol = 10) # Example data

  # Define mock function for FUN
  mock_fun <- function(perm_id, n_perm) {
    return(list(result = perm_id, n_perm = n_perm))
  }

  # Define temporary output directory
  temp_output_dir <- tempdir()

  # Call job_splitter with mock data
  job_splitter(
    job_id = 1,
    n_iter_per_job = 5,
    iter_total = 10,
    output_dir = temp_output_dir,
    FUN = mock_fun,
    subset_vars = list(perm_id = perm_id),
    subset_total_var = "n_perm"
  )

  # Check if output file was created
  output_file <- file.path(temp_output_dir, "res_job_1.rds")
  expect_true(file.exists(output_file))

  # Check contents of the file
  res <- readRDS(output_file)
  expect_equal(nrow(res$result), 10)
  expect_equal(ncol(res$result), 5)
  expect_equal(res$n_perm, 5)
})

test_that("job_splitter throws error when output_dir is missing", {
  # Expect an error when output_dir is not provided
  expect_error(
    job_splitter(
      job_id = 1,
      n_iter_per_job = 5,
      iter_total = 10,
      FUN = function() NULL
    ),
    "output_dir must be specified."
  )
})
