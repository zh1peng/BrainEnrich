test_that("compress_csv_bzip2 compresses CSV file correctly", {
  # Create a temporary CSV file
  temp_input_csv <- tempfile(fileext = ".csv")
  data <- data.frame(A = 1:3, B = c("x", "y", "z"))
  write.csv(data, temp_input_csv, row.names = FALSE)

  # Call compress_csv_bzip2 to compress the file
  BrainEnrich:::compress_csv_bzip2(temp_input_csv)

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
  mock_ask_user_continue <- function(msg) {
    return(TRUE)
  }
  local_mocked_bindings(ask_user_continue = mock_ask_user_continue, .package = "BrainEnrich")
  expect_true(ask_user_continue("Test message")) # Expect FALSE for "N"
})

test_that("ask_user_continue works with valid inputs", {
  # Mock the readline function to return "Y"
  mock_ask_user_continue <- function(msg) {
    return(FALSE)
  }
  local_mocked_bindings(ask_user_continue = mock_ask_user_continue, .package = "BrainEnrich")
  expect_false(ask_user_continue("Test message")) # Expect FALSE for "N"
})

test_that("ask_user_continue returns TRUE for 'Y' input", {
  # Mock the package-local readline wrapper so the installed package sees it.
  local_mocked_bindings(
    be_readline = function(prompt) {
      return("Y")
    },
    .package = "BrainEnrich"
  )

  expect_true(ask_user_continue("Test message"))
})



test_that("ask_user_continue returns FALSE for 'N' input", {
  # Mock the package-local readline wrapper so the installed package sees it.
  local_mocked_bindings(
    be_readline = function(prompt) {
      return("N")
    },
    .package = "BrainEnrich"
  )

  expect_false(ask_user_continue("Test message"))
})

test_that("be_resolve_n_cores respects R CMD check core limits", {
  withr::local_envvar(c("_R_CHECK_LIMIT_CORES_" = "true"))

  available <- parallel::detectCores()
  if (is.na(available) || available < 1L) {
    available <- 1L
  }
  expected_default <- max(min(available - 1L, 2L), 1L)

  resolved_default <- BrainEnrich:::be_resolve_n_cores(0)
  resolved_explicit <- BrainEnrich:::be_resolve_n_cores(4)

  expect_equal(resolved_default, expected_default)
  expect_gte(resolved_explicit, 1L)
  expect_lte(resolved_explicit, 2L)
})
