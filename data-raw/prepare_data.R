# Process PC1 data, filter for regions starting with 'L_' and set row names
brain_data_PC1 <- read.csv("data-raw/desikan_PC1_data.csv") %>%
  filter(stringr::str_detect(Region, "^L_")) %>%
  tibble::column_to_rownames("Region")
usethis::use_data(brain_data_PC1, overwrite = TRUE, compress = "xz")

# Load coordinate data and save it
coord_dk_lh <- read.csv("data-raw/desikan_centroid.csv", stringsAsFactors = FALSE) %>%
  tibble::column_to_rownames("Row")
usethis::use_data(coord_dk_lh, overwrite = TRUE, compress = "xz")


# var_order <- read.csv("data-raw/desikan_PC1_data.csv") %>%
#   filter(stringr::str_detect(Region, "^L_")) %>%
#   select(Region)
# rows <- 34
# cols <- 100
# # Generate a 34 by 100 matrix with uniform random numbers between 0 and 3
# brain_data_random <- matrix(runif(rows * cols, min = 0, max = 3), nrow = rows, ncol = cols)
# colnames(brain_data_random) <- paste0("sub_", 1:cols)
# rownames(brain_data_random) <- var_order$Region
# usethis::use_data(brain_data_random, overwrite = TRUE, compress = "xz")

perm_id_dk_lh_5000 <- readRDS(file.path("E:/xhmhc/BrainEnrich_ms/data", "perm_id_dk_lh_5000.RDS"))
usethis::use_data(perm_id_dk_lh_5000, overwrite = TRUE, compress = "xz")
