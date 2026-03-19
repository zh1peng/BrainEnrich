#' Load gene set annotation data
#'
#' This function loads annotation data from RDS files in package data resources.
#' If data is not found in packaged resources, it attempts to download the file
#' to the user cache directory from the BrainEnrichData companion package.
#'
#' @param type A character string specifying the type of gene set to load.
#' @return An object of class `EnrichAnno`.
#' @importFrom utils download.file
#' @export
#'
#' @examples
#' \dontrun{
#' annoData <- get_annoData("GO_BP")
#' }
get_annoData <- function(type = c(
                           "CellTypes_Lake2018", "CellTypes_Martins2021", "CellTypes_Seidlitz2020",
                           "DGN", "GO_BP", "GO_CC", "GO_MF", "KEGG", "MKEGG", "Reactome", "SynGO", "WikiPathways",
                           "MeSH_A", "MeSH_B", "MeSH_C", "MeSH_D", "MeSH_E", "MeSH_F", "MeSH_G",
                           "MeSH_H", "MeSH_I", "MeSH_J", "MeSH_K", "MeSH_L", "MeSH_M", "MeSH_N", "MeSH_Z"
                         )) {
  type <- match.arg(type)

  relpath <- file.path("geneSets", paste0(type, ".rds"))
  resolved <- be_resolve_data_file(relpath)
  anno_file <- resolved$path

  if (!file.exists(anno_file)) {
    url <- be_data_download_url(relpath)
    message(sprintf("File not found locally. Downloading to user cache... %s", url))
    message(sprintf("Cache target: %s", anno_file))
    options(timeout = max(600, getOption("timeout")))
    tryCatch(
      {
        utils::download.file(url, anno_file, method = "libcurl")
        message("File successfully downloaded.")
      },
      error = function(e) {
        stop(
          "Download failed for annotation file. ",
          "Install BrainEnrichData or download the file manually from: ", url,
          " and place it at: ", anno_file
        )
      }
    )
  }

  raw_obj <- tryCatch(
    readRDS(anno_file),
    error = function(e) {
      stop("Failed to read annotation RDS file: ", conditionMessage(e))
    }
  )

  if (is_EnrichAnno(raw_obj)) {
    return(raw_obj)
  }
  if (is.environment(raw_obj)) {
    return(anno_env_to_EnrichAnno(raw_obj))
  }
  if (is.list(raw_obj) && all(c("term2gene", "term2name") %in% names(raw_obj))) {
    return(new_EnrichAnno(
      term2gene = raw_obj$term2gene,
      term2name = raw_obj$term2name,
      meta = if ("meta" %in% names(raw_obj)) raw_obj$meta else list(source = "list_rds")
    ))
  }

  stop("Unsupported annotation file format. Expected EnrichAnno/list/environment.")
}


#' Get gene set list
#'
#' @param annoData Annotation object of class `EnrichAnno`.
#' @return A named list of gene sets.
#' @export
get_geneSetList <- function(annoData) {
  anno_to_geneSetList(annoData)
}


#' Get gene set descriptions
#'
#' @param term Term IDs.
#' @param annoData An object of class `EnrichAnno`.
#' @param strip_prefix A character string to remove from the beginning of each term.
#' @return A character vector of gene set descriptions.
#' @export
get_termDescription <- function(term, annoData, strip_prefix = "") {
  term <- sapply(term, function(x) gsub(strip_prefix, "", x))
  anno_term_description(term, annoData)
}


#' Filter gene set list
#'
#' This function filters a list of gene sets based on background genes and size.
#'
#' @param bg_genes Character vector of background genes.
#' @param geneSetList Named list of gene sets.
#' @param minGSSize Minimum size threshold.
#' @param maxGSSize Maximum size threshold.
#' @return Filtered named list of gene sets.
#' @export
filter_geneSetList <- function(bg_genes, geneSetList, minGSSize, maxGSSize) {
  if (!is.character(bg_genes)) {
    stop("bg_genes must be a character vector.")
  }
  if (!is.list(geneSetList)) {
    stop("geneSetList must be a list.")
  }
  if (!is.numeric(minGSSize) || !is.numeric(maxGSSize)) {
    stop("minGSSize and maxGSSize must be numeric.")
  }

  filtered <- lapply(geneSetList, function(gs) {
    intersect(unique(as.character(gs)), bg_genes)
  })
  gs_size <- vapply(filtered, length, FUN.VALUE = integer(1))
  keep <- gs_size >= minGSSize & gs_size <= maxGSSize
  filtered[keep]
}


#' Split annotation data
#'
#' Convert an annotation object to path2gene/path2name data frames.
#'
#' @param AnnoEnv An `EnrichAnno` object or legacy annotation environment.
#' @return A list with `path2gene` and `path2name` data frames.
#' @importFrom utils stack
#' @export
split_Anno <- function(AnnoEnv) {
  anno <- NULL
  if (is_EnrichAnno(AnnoEnv)) {
    anno <- AnnoEnv
  } else if (is.environment(AnnoEnv)) {
    anno <- anno_env_to_EnrichAnno(AnnoEnv)
  } else {
    stop("AnnoEnv must be EnrichAnno or a legacy annotation environment.")
  }

  path2gene <- anno$term2gene
  names(path2gene) <- c("pathID", "geneID")
  path2name <- anno$term2name
  names(path2name) <- c("pathID", "pathName")
  list(path2gene = path2gene, path2name = path2name)
}


#' Convert annotation object to pathway table
#'
#' @param AnnoEnv An `EnrichAnno` object or legacy annotation environment.
#' @importFrom stats aggregate
#' @return Data frame with pathway metadata and concatenated genes.
#' @export
Anno2Table <- function(AnnoEnv) {
  tmp.list <- split_Anno(AnnoEnv)
  path2gene <- tmp.list$path2gene
  path2name <- tmp.list$path2name

  gene_concat <- aggregate(geneID ~ pathID, path2gene, paste, collapse = ";")
  gene_count <- aggregate(geneID ~ pathID, path2gene, length)
  names(gene_count)[2] <- "geneSetSize"

  df <- merge(gene_concat, gene_count, by = "pathID")
  df <- merge(df, path2name, by = "pathID")
  df[, c("pathID", "pathName", "geneSetSize", "geneID")]
}


#' Convert a data frame to EnrichAnno
#'
#' @param df A data frame with columns `pathID`, `pathName`, and `geneID`.
#' @param sep Gene ID separator in `geneID`.
#' @return An object of class `EnrichAnno`.
#' @importFrom dplyr select mutate distinct %>%
#' @importFrom tidyr unnest
#' @importFrom rlang .data
#' @export
Table2Anno <- function(df, sep = ";") {
  required_cols <- c("pathID", "pathName", "geneID")
  if (!all(required_cols %in% names(df))) {
    stop("The dataframe must contain the columns: pathID, pathName, geneID.")
  }

  term2gene <- df %>%
    dplyr::select(term_id = pathID, gene_id = geneID) %>%
    dplyr::mutate(gene_id = strsplit(as.character(.data$gene_id), sep)) %>%
    tidyr::unnest(cols = .data$gene_id)

  term2name <- df %>%
    dplyr::select(term_id = pathID, term_name = pathName) %>%
    dplyr::distinct()

  new_EnrichAnno(term2gene = term2gene, term2name = term2name, meta = list(source = "Table2Anno"))
}


#' Filter and intersect genes in annotation table
#'
#' @param df Data frame with `pathID`, `pathName`, `geneID`.
#' @param bg_genes Character vector of background genes.
#' @param sep Gene separator in `geneID`.
#' @param minGSSize Minimum retained gene set size.
#' @param maxGSSize Maximum retained gene set size.
#' @return Filtered annotation table.
#' @importFrom dplyr mutate filter group_by summarise select n_distinct
#' @importFrom tidyr unnest
#' @importFrom rlang .data
#' @export
FilterTable <- function(df, bg_genes, sep = ";", minGSSize = 0, maxGSSize = Inf) {
  required_cols <- c("pathID", "pathName", "geneID")
  if (!all(required_cols %in% names(df))) {
    stop("The dataframe must contain the columns: pathID, pathName, geneID.")
  }
  if (!is.character(bg_genes)) {
    stop("bg_genes must be a character vector containing gene identifiers.")
  }
  if (!is.numeric(minGSSize) || !is.numeric(maxGSSize)) {
    stop("minGSSize and maxGSSize must be numeric values.")
  }

  filtered_df <- df %>%
    dplyr::mutate(gene = strsplit(as.character(.data$geneID), sep)) %>%
    tidyr::unnest(cols = .data$gene) %>%
    dplyr::filter(.data$gene %in% bg_genes) %>%
    dplyr::group_by(.data$pathID, .data$pathName) %>%
    dplyr::summarise(
      geneID = paste(unique(.data$gene), collapse = sep),
      geneSetSize = dplyr::n_distinct(.data$gene),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$geneSetSize >= minGSSize & .data$geneSetSize <= maxGSSize) %>%
    dplyr::select(pathID, pathName, geneSetSize, geneID)

  filtered_df
}
