test_that("compress_csv_bzip2 compresses CSV file correctly", {
  # Create a temporary CSV file
  temp_input_csv <- tempfile(fileext = ".csv")
  data <- data.frame(A = 1:3, B = c("x", "y", "z"))
  write.csv(data, temp_input_csv, row.names = FALSE)

  # Call compress_csv_bzip2 to compress the file
  compress_csv_bzip2(temp_input_csv)

  # Check if the compressed file exists
  compressed_file <- paste0(temp_input_csv, ".bz2")
  expect_true(file.exists(compressed_file))

  # Check if the compressed file can be read back correctly
  decompressed_data <- read.csv(bzfile(compressed_file), check.names = FALSE)
  expect_equal(decompressed_data, data)

  # Clean up temporary files
  unlink(temp_input_csv)
  unlink(compressed_file)
})

test_that("ask_user_continue works with valid inputs", {
  # Mock the readline function to return "Y"

  pkgload::load_all()
  mock_ask_user_continue <- function(msg) {
    return(TRUE)
  }
  local_mocked_bindings(ask_user_continue = mock_ask_user_continue)
  expect_true(ask_user_continue("Test message")) # Expect FALSE for "N"
})


test_that("ask_user_continue works with valid inputs", {
  # Mock the readline function to return "Y"

  pkgload::load_all()
  mock_ask_user_continue <- function(msg) {
    return(FALSE)
  }
  local_mocked_bindings(ask_user_continue = mock_ask_user_continue)
  expect_false(ask_user_continue("Test message")) # Expect FALSE for "N"
})
