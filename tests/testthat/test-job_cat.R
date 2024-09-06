test_that("job_cat combines results from multiple RDS files", {
  pkgload::load_all()
  mock_ask_user_continue <- function(msg) {
    return(TRUE)
  }

  local_mocked_bindings(ask_user_continue = mock_ask_user_continue)

  # Create mock data and save as multiple RDS files
  temp_input_dir <- tempdir()
  files_in_dir <- list.files(temp_input_dir, full.names = TRUE)
  if (length(files_in_dir) > 0) {
    unlink(files_in_dir)
  }

  # Save some example RDS files
  for (i in 1:3) {
    saveRDS(list(result = rnorm(5)), file = file.path(temp_input_dir, sprintf("res_job_%d.rds", i)))
  }

  # Call job_cat to combine the RDS files
  combined_res <- job_cat(
    input_dir = temp_input_dir,
    file_pattern = "res_job_%d.rds",
    save_name = "combined_results.rds"
  )

  # Check if combined results file exists
  combined_file <- file.path(temp_input_dir, "combined_results.rds")
  expect_true(file.exists(combined_file))

  # Check that the combined results have the correct length (3 lists combined)
  expect_equal(length(combined_res), 3)
})



test_that("job_cat combines results from multiple RDS files", {
  pkgload::load_all()
  mock_ask_user_continue <- function(msg) {
    return(FALSE)
  }

  local_mocked_bindings(ask_user_continue = mock_ask_user_continue)

  # Create mock data and save as multiple RDS files
  temp_input_dir <- tempdir()

  files_in_dir <- list.files(temp_input_dir, full.names = TRUE)
  if (length(files_in_dir) > 0) {
    unlink(files_in_dir)
  }

  # Save some example RDS files
  for (i in 1:3) {
    saveRDS(list(result = rnorm(5)), file = file.path(temp_input_dir, sprintf("res_job_%d.rds", i)))
  }

  # Call job_cat to combine the RDS files
  expect_error(
    combined_res <- job_cat(
      input_dir = temp_input_dir,
      file_pattern = "res_job_%d.rds",
      save_name = "combined_results.rds"
    ), "Operation cancelled by user."
  )
})




test_that("job_cat detects missing RDS files", {
  temp_input_dir <- tempdir()
  files_in_dir <- list.files(temp_input_dir, full.names = TRUE)
  if (length(files_in_dir) > 0) {
    unlink(files_in_dir)
  }

  # Save only two RDS files out of three expected
  for (i in 1:2) {
    saveRDS(list(result = rnorm(5)), file = file.path(temp_input_dir, sprintf("res_job_%d.rds", i)))
  }

  # Expect an error due to missing file
  expect_error(
    job_cat(
      input_dir = temp_input_dir,
      n_rds = 3, # Expecting 3 files
      file_pattern = "res_job_%d.rds"
    ),
    "Some RDS files are missing."
  )
})

test_that("job_cat deletes original RDS files after combining", {
  # Create mock data and save as multiple RDS files
  temp_input_dir <- tempdir()
  files_in_dir <- list.files(temp_input_dir, full.names = TRUE)
  if (length(files_in_dir) > 0) {
    unlink(files_in_dir)
  }

  # Save some example RDS files
  for (i in 1:3) {
    saveRDS(list(result = rnorm(5)), file = file.path(temp_input_dir, sprintf("res_job_%d.rds", i)))
  }

  # Call job_cat and check if original files are deleted
  job_cat(
    input_dir = temp_input_dir,
    n_rds = 3,
    file_pattern = "res_job_%d.rds",
    delete_originals = TRUE,
    save_name = "combined_results.rds"
  )

  # Check that the original RDS files have been deleted
  for (i in 1:3) {
    expect_false(file.exists(file.path(temp_input_dir, sprintf("res_job_%d.rds", i))))
  }
})
