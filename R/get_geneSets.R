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
#'     \item{"MeSH_A"}{MeSH category A gene sets}
#'     \item{"MeSH_B"}{MeSH category B gene sets}
#'     \item{"MeSH_C"}{MeSH category C gene sets}
#'     \item{"MeSH_D"}{MeSH category D gene sets}
#'     \item{"MeSH_E"}{MeSH category E gene sets}
#'     \item{"MeSH_F"}{MeSH category F gene sets}
#'     \item{"MeSH_G"}{MeSH category G gene sets}
#'     \item{"MeSH_H"}{MeSH category H gene sets}
#'     \item{"MeSH_I"}{MeSH category I gene sets}
#'     \item{"MeSH_J"}{MeSH category J gene sets}
#'     \item{"MeSH_K"}{MeSH category K gene sets}
#'     \item{"MeSH_L"}{MeSH category L gene sets}
#'     \item{"MeSH_M"}{MeSH category M gene sets}
#'     \item{"MeSH_N"}{MeSH category N gene sets}
#'     \item{"MeSH_Z"}{MeSH category Z gene sets}
#'   }
#' @importFrom utils download.file
#' @return A data frame containing the annotation data.
#' @export
#'
#' @examples
#' \dontrun{
#' annoData <- get_annoData("GO_BP")
#' }
get_annoData <- function(type = c(
                           "CellTypes_Lake2018", "CellTypes_Martins2021", "CellTypes_Seidlitz2020",
                           "DGN", "GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "SynGO", "WikiPathways",
                           "MeSH_A", "MeSH_B", "MeSH_C", "MeSH_D", "MeSH_E", "MeSH_F", "MeSH_G",
                           "MeSH_H", "MeSH_I", "MeSH_J", "MeSH_K", "MeSH_L", "MeSH_M", "MeSH_N", "MeSH_Z"
                         )) {
  type <- match.arg(type)

  # Determine the package installation directory
  package_dir <- find.package("BrainEnrich")
  # Specify the path for the 'extdata' directory
  extdata_dir <- file.path(package_dir, "extdata")

  # Create the 'extdata' and 'geneSets' directories if they do not exist
  if (!dir.exists(extdata_dir)) {
    dir.create(extdata_dir)
  }

  gene_set_dir <- file.path(extdata_dir, "geneSets")
  if (!dir.exists(gene_set_dir)) {
    dir.create(gene_set_dir)
  }

  # Define the local file path for the gene set
  GeneSetsRDS <- file.path(gene_set_dir, paste0(type, ".rds"))

  # Define GitHub URL for downloading the file
  url <- paste0("https://github.com/zh1peng/BrainEnrich/raw/master/inst/extdata/geneSets/", type, ".rds")

  # Download the file if it does not exist locally
  if (!file.exists(GeneSetsRDS)) {
    message(sprintf("File not found locally. Downloading from GitHub... %s", url))
    message("If the download is slow, download manually.")
    message(sprintf("and save files as %s", GeneSetsRDS))

    tryCatch(
      {
        download.file(url, GeneSetsRDS, method = "libcurl")
        message("File successfully downloaded.")
      },
      error = function(e) {
        message("Download failed. Please copy and paste the following URL into your browser:")
        message(url)
        message(sprintf("and save the file to the following folder: %s", gene_set_dir))
      }
    )
  }

  message("Loading annotation data...")
  annoData <- tryCatch(
    {
      # Attempt to read the RDS file
      readRDS(GeneSetsRDS)
    },
    error = function(e) {
      # If an error occurs (e.g., file is corrupted), handle it here
      message("Error reading the RDS file: ", conditionMessage(e))
      message("The file might be incomplete or corrupted. Removing the file.")

      # Attempt to remove the corrupted file
      if (file.exists(GeneSetsRDS)) {
        file.remove(GeneSetsRDS)
      }
      # Stop execution and raise an error
      stop("Failed to read the RDS file. The file has been removed. Please try again.")
    }
  )

  return(annoData)
}


#' Get Gene Set List
#'
#' This function retrieves a gene set list from annotation data. It optionally converts
#' gene identifiers to gene symbols.
#'
#' @param annoData Annotation data to retrieve gene sets from.
#' @return A list of gene sets.
#' @import DOSE
#' @importFrom utils getFromNamespace
#' @export
get_geneSetList <- function(annoData) {
  getGeneSet <- getFromNamespace("getGeneSet", "DOSE")
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
#' @importFrom utils getFromNamespace
#' @export
get_termDescription <- function(term, annoData, strip_prefix = "") {
  term <- sapply(term, function(x) gsub(strip_prefix, "", x))
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
#' @importFrom utils getFromNamespace
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



#' Split Annotation Data
#'
#' This function extracts and reconstructs the `path2gene` and `path2name` data frames
#' from an annotation environment created by the `build_Anno` function. It retrieves
#' and organizes pathway IDs and associated gene IDs, as well as pathway names, into
#' two separate data frames.
#'
#' @param AnnoEnv An environment containing annotation data, including
#'        `PATHID2EXTID` and `PATHID2NAME`, created by the `build_Anno` function.
#' @importFrom utils stack
#' @return A list with two components:
#' \describe{
#'   \item{path2gene}{A data frame with two columns: `pathID` and `geneID`, representing
#'   the relationship between pathway IDs and gene IDs.}
#'   \item{path2name}{A data frame with two columns: `pathID` and `pathName`, representing
#'   the mapping between pathway IDs and their corresponding names. If `PATHID2NAME` is not
#'   present in the environment, this component will be `NULL`.}
#' }
#' @export
split_Anno <- function(AnnoEnv) {
  if (!exists("PATHID2EXTID", envir = AnnoEnv) ||
    !exists("EXTID2PATHID", envir = AnnoEnv)) {
    stop("The environment does not contain the required objects.")
  }
  PATHID2EXTID <- get("PATHID2EXTID", envir = AnnoEnv)
  PATHID2NAME <- get("PATHID2NAME", envir = AnnoEnv, inherits = TRUE)

  # Reconstruct path2gene from PATHID2EXTID
  path2gene <- stack(PATHID2EXTID)
  path2gene <- path2gene[, c(2, 1)]
  colnames(path2gene) <- c("pathID", "geneID")

  # Reconstruct path2name from PATHID2NAME
  if (!is.null(PATHID2NAME)) {
    path2name <- stack(PATHID2NAME)
    path2name <- path2name[, c(2, 1)]
    colnames(path2name) <- c("pathID", "pathName")
  } else {
    path2name <- NULL
  }
  return(list(path2gene = path2gene, path2name = path2name))
}

#' Convert Annotation Data to Pathway Table
#'
#' This function takes an annotation environment, typically from `clusterProfiler`,
#' and converts it into a pathway table. Each pathway will have its associated genes
#' concatenated into a single string, along with the pathway name and the size of the gene set.
#'
#' @param AnnoEnv An environment or list containing annotation data.
#' This environment is expected to contain pathway-to-gene (`path2gene`) and pathway-to-name
#' (`path2name`) mappings.
#' @importFrom stats aggregate
#' @return A data frame with four columns:
#' \describe{
#'   \item{pathID}{Pathway identifier}
#'   \item{pathName}{Pathway name}
#'   \item{geneSetSize}{Size of the gene set for each pathway}
#'   \item{geneID}{Concatenated string of gene IDs for each pathway, separated by ";"}
#' }

#' @export
Anno2Table <- function(AnnoEnv) {
  # Split annotation environment
  tmp.list <- split_Anno(AnnoEnv)
  path2gene <- tmp.list$path2gene
  path2name <- tmp.list$path2name

  # Concatenate genes by pathID
  gene_concat <- aggregate(geneID ~ pathID, path2gene, paste, collapse = ";")

  # Count the number of genes per pathID
  gene_count <- aggregate(geneID ~ pathID, path2gene, length)
  names(gene_count)[2] <- "geneSetSize" # Rename the count column

  # Merge gene concatenation and gene count with path names
  df <- merge(gene_concat, gene_count, by = "pathID")
  df <- merge(df, path2name, by = "pathID")

  # Reorder columns if needed
  df <- df[, c("pathID", "pathName", "geneSetSize", "geneID")]

  return(df)
}




#' Convert a Dataframe to Annotation Environment
#' This function converts a dataframe containing annotation information back to an annotation environment
#' using the `build_Anno` function from the `DOSE` package. The resulting environment can be used for
#' enrichment analysis or other annotation-related tasks.
#' @param df A dataframe containing the annotation information with columns:
#'   - `pathID`: The ID of the pathway/term.
#'   - `pathName`: The descriptive name of the pathway/term.
#'   - `geneID`: A string of gene identifiers associated with the pathway/term, separated by a delimiter.
#' @param sep A character string that defines the delimiter for splitting the `geneID` column. Defaults to `";"`.
#' @return An annotation environment containing the term-to-gene mapping and term descriptions.
#' @importFrom dplyr select mutate distinct %>%
#' @importFrom tidyr unnest
#' @importFrom rlang .data 
#' @export
Table2Anno <- function(df, sep = ";") {
  # Check if required columns are present
  required_cols <- c("pathID", "pathName", "geneID")
  if (!all(required_cols %in% names(df))) {
    stop("The dataframe must contain the columns: pathID, pathName, and geneID.")
  }

  # Split geneID into individual genes using the specified separator
  TERM2GENE <- df %>%
    dplyr::select(term = .data$pathID, gene = .data$geneID) %>%
    dplyr::mutate(gene = strsplit(as.character(.data$gene), sep)) %>%
    tidyr::unnest(cols = .data$gene)

  # Prepare TERM2NAME
  TERM2NAME <- df %>%
    dplyr::select(term = .data$pathID, description = .data$pathName) %>%
    dplyr::distinct()

  # Use build_Anno to reconstruct the annotation object
  build_Anno <- getFromNamespace("build_Anno", "DOSE")
  annoEnv <- build_Anno(TERM2GENE, TERM2NAME)
  return(annoEnv)
}

#' Filter and Intersect Genes in Annotation Data
#'
#' This function filters the genes in the `geneID` column of an annotation dataframe
#' based on a given list of background genes. Only genes present in the background list
#' are retained in the resulting dataframe. It also calculates and includes the `geneSetSize` column,
#' and filters gene sets based on specified size constraints.
#'
#' @param df A dataframe containing the annotation information with columns:
#'   - `pathID`: The ID of the pathway/term.
#'   - `pathName`: The descriptive name of the pathway/term.
#'   - `geneID`: A string of gene identifiers associated with the pathway/term, separated by a delimiter.
#' @param background A character vector of background gene names.
#' @param sep A character string that defines the delimiter for splitting the `geneID` column. Defaults to `";"`.
#' @param minGSSize An integer specifying the minimum gene set size to retain. Defaults to `0`.
#' @param maxGSSize An integer specifying the maximum gene set size to retain. Defaults to `Inf`.
#'
#' @return A dataframe with columns:
#'   - `pathID`: The ID of the pathway/term.
#'   - `pathName`: The descriptive name of the pathway/term.
#'   - `geneSetSize`: The size of the filtered gene set.
#'   - `geneID`: The filtered and concatenated gene identifiers.
#'
#' @importFrom dplyr mutate filter group_by summarise
#' @importFrom tidyr unnest
#' @export
FilterTable <- function(df, bg_genes, sep = ";", minGSSize = 0, maxGSSize = Inf) {
  # Check if required columns are present
  required_cols <- c("pathID", "pathName", "geneID")
  if (!all(required_cols %in% names(df))) {
    stop("The dataframe must contain the columns: pathID, pathName, and geneID.")
  }
  
   # Check if bg_genes is a character vector
  if (!is.character(bg_genes)) {
    stop("bg_genes must be a character vector containing gene identifiers.")
  }
  
  # Check if minGSSize and maxGSSize are numeric
  if (!is.numeric(minGSSize) || !is.numeric(maxGSSize)) {
    stop("minGSSize and maxGSSize must be numeric values.")
  }

  # Split geneID into individual genes, filter against background, and reconstruct geneID
  filtered_df <- df %>%
    dplyr::mutate(
      gene = strsplit(as.character(.data$geneID), sep)  # Split geneID into lists
    ) %>%
    tidyr::unnest(cols = .data$gene) %>%  # Expand into long format
    dplyr::filter(.data$gene %in% background) %>%  # Keep only genes in the background
    dplyr::group_by(.data$pathID, .data$pathName) %>%
    dplyr::summarise(
      geneID = paste(unique(.data$gene), collapse = sep),  # Recombine filtered genes
      geneSetSize = n_distinct(.data$gene),  # Count the number of unique genes
      .groups = "drop"  # Ungroup after summarizing
    ) %>%
    dplyr::filter(
      geneSetSize >= minGSSize & geneSetSize <= maxGSSize  # Filter by gene set size
    ) %>%
    dplyr::select(pathID, pathName, geneSetSize, geneID)  # Reorder columns
  
  return(filtered_df)
}
