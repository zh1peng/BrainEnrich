library(testthat)

test_that("rotate_parcellation works as expected", {
  # Load the example data for left hemisphere coordinates
  data(coord_dk_lh)
  
  # Run the rotate_parcellation function with a specific seed and 10 rotations
  perm_id_dk_lh_10 <- rotate_parcellation(coord.l = coord_dk_lh, nrot = 10, seed = 2024)
  
  # Test that the output is a matrix
  expect_true(is.matrix(perm_id_dk_lh_10))
  
  # Test that the number of rotations (columns) is correct
  expect_equal(ncol(perm_id_dk_lh_10), 10)
  
  # Test that the number of regions (rows) matches the input coordinates
  expect_equal(nrow(perm_id_dk_lh_10), nrow(coord_dk_lh))
  
  # Test reproducibility: running the same function with the same seed should give the same result
  perm_id_dk_lh_10_repeat <- rotate_parcellation(coord.l = coord_dk_lh, nrot = 10, seed = 2024)
  
  # The output should be identical since the seed is set
  expect_equal(perm_id_dk_lh_10, perm_id_dk_lh_10_repeat)
})

