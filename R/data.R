#' Brain Data PC1 for Left Hemisphere
#'
#' This dataset contains PC1 data filtered for regions starting with 'L_' in the Desikan atlas.
#'
#' @format A data frame with rows as regions and columns as effect sizes of case-control comparisons on regional cortical thickness between bipolar disorders and healthy controls.
#' @source read.csv('data-raw/desikan_PC1_data.csv')
"brain_data"


#' Desikan Centroid Coordinates for Left Hemisphere
#'
#' This dataset contains the centroid coordinates for the Desikan atlas regions in the left hemisphere.
#'
#' @format A data frame with rows as regions and columns as coordinates (x, y, z).
#' @source read.csv('data-raw/desikan_centroid.csv')
"coord_dk_lh"


#' Permutation index for left Desikan regions (5000 permutations)
#' @format A matrix with rows as regions and columns as permutated indices.
#' @source rotate_parcellation(coord.l = coord_dk_lh, nrot = 5000)
"perm_id_dk_lh_5000"


#' Simulated data from HCP data
#' @format A data frame with rows as subjects and columns as age, sex, BMI and regional cortical thickness values
#' @source mvrnorm(n = 100, mu = mean_vals, Sigma = cov_matrix); sample(df.hcp$Age_in_Yrs, 100, replace = FALSE); sample(df.hcp$Sex, 100, replace = FALSE); sample(df.hcp$BMI, 100, replace = FALSE)
"sample_df"