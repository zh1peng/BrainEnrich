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
#' @importFrom utils download.file
#' @importFrom dplyr %>% filter everything across
#' @importFrom tibble column_to_rownames
#' @importFrom stats complete.cases
#' @importFrom rlang .data
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
  GeneExpCSV <- file.path(gene_exp_dir, sprintf("%s_%s.csv.bz2", atlas, rdonor))



  # Define GitHub URL for downloading the file
  url <- paste0("https://github.com/zh1peng/BrainEnrich/raw/master/inst/extdata/geneExp/", atlas, "_", rdonor, ".csv.bz2")


  if (!file.exists(GeneExpCSV)) {
    message(sprintf("File not found locally. Downloading from GitHub... %s", url))
    message("If the download is slow, download manually.")
    message(sprintf("and save files as %s", GeneExpCSV))

    GeneExpDir <- dirname(GeneExpCSV)
    if (!dir.exists(GeneExpDir)) {
      dir.create(GeneExpDir, recursive = TRUE)
    }
    options(timeout = 600)
    tryCatch(
      {
        download.file(url, GeneExpCSV)
        message("File successfully downloaded.")
      },
      error = function(e) {
        message("Download failed. Please copy and paste the following URL into your browser:")
        message(url)
        message(sprintf("and save the file to the following folder: %s", GeneExpDir))
      }
    )
  }


  # Read the CSV file
  message("Reading gene expression data...")

  gene.df <- tryCatch(
    {
      # Attempt to read the CSV file
      read.csv_bzip2(GeneExpCSV, stringsAsFactors = FALSE, check.names = FALSE)
    },
    error = function(e) {
      # If an error occurs (e.g., file is corrupted), handle it here
      message("Error reading the compressed CSV file: ", conditionMessage(e))
      message("The file might be incomplete or corrupted. Removing the file.")

      # Attempt to remove the corrupted file
      if (file.exists(GeneExpCSV)) {
        file.remove(GeneExpCSV)
      }

      # Stop execution and raise an error
      stop("Failed to read the CSV file. The file has been removed. Please try again.")
    }
  )
  # Filter based on hemisphere
  if (hem %in% c("L", "R")) {
    gene.df <- dplyr::filter(gene.df, grepl(search_pattern, .data$Region))
  }

  # Filter complete cases and convert to matrix
  gene.df <- gene.df %>%
    dplyr::filter(complete.cases(across(everything()))) %>%
    tibble::column_to_rownames("Region")

  gene.mx <- as.matrix(gene.df)

  return(gene.mx)
}
