#' Brain Data PC1 for Left Hemisphere
#'
#' This dataset contains PC1 data filtered for regions starting with 'L_' in the Desikan atlas.
#'
#' @format A data frame with rows as regions and columns as PC1 data.
#' @source read.csv('data-raw/desikan_PC1_data.csv')
"brain_data_PC1"


#' Desikan Centroid Coordinates for Left Hemisphere
#'
#' This dataset contains the centroid coordinates for the Desikan atlas regions in the left hemisphere.
#'
#' @format A data frame with rows as regions and columns as coordinates (x, y, z).
#' @source read.csv('data-raw/desikan_centroid.csv')
"coord_dk_lh"
