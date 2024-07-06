#' Get Gene Expression Data
#'
#' This function retrieves gene expression data based on specified parameters.
#' #' This function loads gene expression data from CSV files located in the `inst/extdata/geneExp` directory of the package.
#' If the specified file does not exist locally, it will be downloaded from the GitHub repository.
#'
#' The data is obtained from the ENIGMA-TOOLBOX. Please cite the ENIGMA-TOOLBOX:
#' Larivière, S., Paquola, C., Park, B. Y., Royer, J., Wang, Y., Benkarim, O., ... & Bernhardt, B. C. (2021).
#' The ENIGMA Toolbox: multiscale neural contextualization of multisite neuroimaging datasets. Nature Methods, 18(7), 698-700.
#'
#'
#' @param atlas A character string specifying the atlas to use. Options are "desikan", "schaefer100", "schaefer200", "schaefer300".
#' @param rdonor A character string specifying the donor resolution to use. Options are "r0.2", "r0.4", "r0.6".
#' @param hem A character string specifying the hemisphere to use. Options are "L" (Left), "R" (Right), "B" (Both).
#'
#' @return A matrix containing the gene expression data.
#' @export
#'
#' @examples
#' \dontrun{
#' geneExpMatrix <- get_geneExp("desikan", "r0.4", "L")
#' }

get_geneExp <- function(atlas = c("desikan", "schaefer100", "schaefer200", "schaefer300"),
                        rdonor = c("r0.2", "r0.4", "r0.6"),
                        hem = c("L", "R", "B")) {
  atlas <- match.arg(atlas)
  rdonor <- match.arg(rdonor)
  hem <- match.arg(hem)

  if (atlas == "desikan") {
    search_pattern <- sprintf("^%s_", hem)
  } else if (atlas %in% c("schaefer100", "schaefer200", "schaefer300")) {
    search_pattern <- sprintf("%sH", hem)
  }

  # Define the path to the gene expression data
  package_dir <- find.package("BrainEnrich")
  # Specify the path for the new 'extdata' directory
  extdata_dir <- file.path(package_dir, "extdata")
  # Create the 'extdata' directory if it does not exist
  if (!dir.exists(extdata_dir)) {
    dir.create(extdata_dir)
  }

  gene_exp_dir <- file.path(extdata_dir, "geneExp")
    # Ensure the extdata/geneExp directory exists
  if (!dir.exists(gene_exp_dir)) {
    dir.create(gene_exp_dir, recursive = TRUE)
  }
  GeneExpCSV <- file.path(gene_exp_dir, sprintf("%s_%s.csv", atlas, rdonor))



  # Define GitHub URL for downloading the file
  url <- paste0("https://github.com/zh1peng/BrainEnrich/raw/master/extdata/geneExp/", atlas, "_", rdonor, ".csv")
  
  if (!file.exists(GeneExpCSV)) {
      options(timeout = 300) 
      download.file(url, GeneExpCSV,method='libcurl')
  }

  # Read the CSV file
  gene.df <- read.csv(GeneExpCSV, stringsAsFactors = FALSE)

  # Filter based on hemisphere
  if (hem %in% c("L", "R")) {
    gene.df <- dplyr::filter(gene.df, grepl(search_pattern, Region))
  }

  # Filter complete cases and convert to matrix
  gene.df <- gene.df %>%
    dplyr::filter(complete.cases(.)) %>%
    tibble::column_to_rownames("Region")

  gene.mx <- as.matrix(gene.df)

  return(gene.mx)
}
