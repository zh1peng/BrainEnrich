#' Load gene sets annotation data
#'
#' This function loads annotation data from RDS files located in the `inst/extdata/geneSets` directory of the package.
#' If the specified file does not exist locally, it will be downloaded from the GitHub repository.
#'
#' @param type A character string specifying the type of gene set to load. Options are:
#'   \describe{
#'     \item{"CellTypes_Lake2018"}{Cell types data from Lake2018}
#'     \item{"CellTypes_Martins2021"}{Cell types data from Martins2021}
#'     \item{"CellTypes_Seidlitz2020"}{Cell types data from Seidlitz2020}
#'     \item{"DGN"}{DisGeNET gene sets}
#'     \item{"GO_BP"}{Gene Ontology Biological Process}
#'     \item{"GO_CC"}{Gene Ontology Cellular Component}
#'     \item{"GO_MF"}{Gene Ontology Molecular Function}
#'     \item{"KEGG"}{KEGG gene sets}
#'     \item{"Reactome"}{Reactome gene sets}
#'     \item{"SynGO"}{SynGO gene sets}
#'     \item{"WikiPathways"}{WikiPathways gene sets}
#'   }
#'
#' @return A data frame containing the annotation data.
#' @export
#'
#' @examples
#' \dontrun{
#' annoData <- get_annoData("GO_BP")
#' }

get_annoData <- function(type = c("CellTypes_Lake2018", "CellTypes_Martins2021", "CellTypes_Seidlitz2020",
                                  "DGN", "GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "SynGO", "WikiPathways")) {
  type <- match.arg(type)
  
  # Define local file path
  local_file <- system.file("extdata", "geneSets", paste0(type, ".rds"), package = "BrainEnrich")

  # Define GitHub URL for downloading the file
  url <- paste0("https://github.com/zh1peng/BrainEnrich/tree/master/inst/extdata/geneSets/", type, ".rds")
  
  if (!file.exists(local_file)) {
    message("File not found locally. Downloading from GitHub...")
    temp_file <- tempfile(fileext = ".rds")
    
    # Clean-up function to ensure temp files are removed
    on.exit({
      if (file.exists(temp_file)) unlink(temp_file)
    }, add = TRUE)
    
    tryCatch({
      download.file(url, temp_file, mode = "wb")
      file.copy(temp_file, local_file)
    }, error = function(e) {
      stop("An error occurred while downloading the file: ", e$message)
    })
  }
  
  message("Loading annotation data...")
  annoData <- readRDS(local_file)
  
  return(annoData)
}


#' Filter Gene Set List
#'
#' This function filters a list of gene sets based on the background genes and specified size constraints.
#'
#' @param bg_genes A vector of background gene symbols to be used for filtering.
#' @param geneSetList A list of gene sets to be filtered.
#' @param minGSSize Minimum gene set size for filtering.
#' @param maxGSSize Maximum gene set size for filtering.
#' @return A filtered list of gene sets that meet the size constraints and background genes criteria.
#' @import DOSE
#' @export
filter_geneSetList <- function(bg_genes, geneSetList, minGSSize, maxGSSize) {
  # Create temporary values to name the background genes
  tmp.val <- seq_along(bg_genes)
  names(tmp.val) <- bg_genes
  # Filter the gene set list
  geneSet_filter <- getFromNamespace("geneSet_filter", "DOSE")
  geneSetList_filtered <- geneSet_filter(geneSetList, tmp.val, minGSSize, maxGSSize)
  return(geneSetList_filtered)
}



