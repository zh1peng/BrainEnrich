# Process PC1 data, filter for regions starting with 'L_' and set row names
# brain_data <- read.csv("data-raw/desikan_PC1_data.csv") %>%
#   filter(stringr::str_detect(Region, "^L_")) %>%
#   tibble::column_to_rownames("Region")
# usethis::use_data(brain_data, overwrite = TRUE, compress = "xz")


df <- read.csv(file.path("E:/xhmhc/BrainEnrich_ms/data", "MP_data.csv")) %>% dplyr::filter(stringr::str_detect(label, "L_"))
brain_data <- df %>%
  dplyr::select(label, BD) %>%
  dplyr::filter(stringr::str_detect(label, "L_")) %>%
  tibble::column_to_rownames("label")
usethis::use_data(brain_data, overwrite = TRUE, compress = "xz")


res <- readRDS(file.path("E:/xhmhc/BrainEnrich_ms/results", "bd_MF_res.RDS"))
usethis::use_data(res, overwrite = TRUE, compress = "xz")

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


# simulate data from hcp data
brain_data <- df.hcp %>%
  tibble::column_to_rownames("SubjID") %>%
  dplyr::select(starts_with("L_") & ends_with("thickavg")) %>%
  dplyr::rename_all(~ stringr::str_replace_all(., "_thickavg$", "")) %>%
  t()

mean_vals <- rowMeans(brain_data)
cov_matrix <- cov(t(brain_data))
set.seed(2024)
sample_df <- MASS::mvrnorm(n = 100, mu = mean_vals, Sigma = cov_matrix) %>% as.data.frame()
colnames(sample_df) <- names(mean_vals)

sample_df$Age <- sample(df.hcp$Age_in_Yrs, 100, replace = FALSE)
sample_df$Sex <- sample(df.hcp$Sex, 100, replace = FALSE)
sample_df$BMI <- sample(df.hcp$BMI, 100, replace = FALSE)
usethis::use_data(sample_df, overwrite = TRUE, compress = "xz")
