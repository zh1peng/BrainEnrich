data(brain_data)

test_that("plot_brain works with valid inputs", {
  df2plot <- brain_data %>% tibble::rownames_to_column("label")
  # Prepare brain_data for plotting
  # Call the plot_brain function and ensure it doesn't throw an error
  expect_no_error({
    plot_brain(df2plot,
      ats = "dk",
      filterby = "none",
      limit2show = c(-0.5, 0.5),
      what2plot = "BD",
      hem = "left",
      low = "#0197b2",
      mid = "white",
      high = "orange",
      legend2show = "Effect Size"
    )
  })
})


# Create a list of regions with left (L) and right (R) labels
regions <- c("thal", "caud", "put", "pal", "hippo", "amyg", "accumb")

# Create labels for left and right hemispheres
labels <- c(paste0("SV_L_", regions), paste0("SV_R_", regions))

# Generate some sample brain data for these regions
set.seed(123)
df2plot <- tibble::tibble(
  label = labels,
  value = runif(length(labels), min = -1, max = 1) # Random values for illustration
)

test_that("plot_brain works with valid inputs", {
  # Prepare brain_data for plotting
  # Call the plot_brain function and ensure it doesn't throw an error
  expect_no_error({
    plot_brain(df2plot,
      ats = "aseg",
      filterby = "none",
      limit2show = c(-0.5, 0.5),
      what2plot = "value",
      hem = "both",
      low = "#0197b2",
      mid = "white",
      high = "orange",
      legend2show = "Effect Size"
    )
  })
})
