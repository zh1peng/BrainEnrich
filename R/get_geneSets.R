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


  # Determine the package installation directory
  package_dir <- find.package("BrainEnrich")
  # Specify the path for the new 'extdata' directory
  extdata_dir <- file.path(package_dir, "extdata")
  # Create the 'extdata' directory if it does not exist
  if (!dir.exists(extdata_dir)) {
    dir.create(extdata_dir)
  }


  gene_set_dir <- file.path(extdata_dir, "geneSets")

   if (!dir.exists(gene_set_dir)) {
      dir.create(gene_set_dir)
    }
  GeneSetsRDS <- file.path(gene_set_dir, paste0(type, ".rds"))

    # Ensure the gene_set_dir directory exists

  # Define GitHub URL for downloading the file
  url <- paste0("https://github.com/zh1peng/BrainEnrich/raw/master/extdata/geneSets/", type, ".rds")
  
  if (!file.exists(GeneSetsRDS)) {
    message(sprintf("File not found locally. Downloading from GitHub... %s", url))
    message("If the download is slow, download manually.")
    message(sprintf("and save files as %s", GeneSetsRDS))
    options(timeout = 600) 
    download.file(url, GeneSetsRDS, method = 'libcurl')
  }
  
  message("Loading annotation data...")
  annoData <- readRDS(GeneSetsRDS)
  
  return(annoData)
}


#' Get Gene Set List
#'
#' This function retrieves a gene set list from annotation data. It optionally converts
#' gene identifiers to gene symbols.
#'
#' @param annoData Annotation data to retrieve gene sets from.
#' @param convert_to_symbol Logical; if TRUE, converts gene identifiers to gene symbols.
#' @return A list of gene sets.
#' @import DOSE
#' @export 
get_geneSetList <- function(annoData) {
  getGeneSet <- getFromNamespace('getGeneSet', 'DOSE')
  geneSetList <- getGeneSet(annoData)
  return(geneSetList)
}

#' Get Gene Set Descriptions
#'
#' This function retrieves descriptions for gene sets from annotation data.
#' @param term term to search from annoData (can be a vector of terms).
#' @param annoData An environment containing annotation data.
#' @param strip_prefix A character string to remove from the beginning of each term. 
#' @return A character vector of gene set descriptions.
#' @import DOSE
#' @export 
get_termDescription <- function(term, annoData, strip_prefix = '') {
  term = sapply(term, function(x) gsub(strip_prefix, '', x))
  TERM2NAME <- getFromNamespace("TERM2NAME", "DOSE")
  termDescrip <- TERM2NAME(term, annoData)
  return(termDescrip)
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



