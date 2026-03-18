#' Get gene expression data
#'
#' This function retrieves expression matrices for supported atlases. It first
#' looks for packaged resources (preferring `BrainEnrichData`), including an
#' installed companion package or a sibling development checkout, and falls
#' back to downloading into the user cache directory.
#'
#' @param atlas Atlas name. One of "desikan", "schaefer100", "schaefer200", "schaefer300".
#' @param rdonor Donor threshold. One of "r0.2", "r0.4", "r0.6".
#' @param hem Hemisphere selection. One of "L", "R", "B".
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
  } else {
    search_pattern <- sprintf("%sH", hem)
  }

  relpath <- file.path("geneExp", sprintf("%s_%s.csv.bz2", atlas, rdonor))
  resolved <- be_resolve_data_file(relpath)
  gene_exp_file <- resolved$path

  if (!file.exists(gene_exp_file)) {
    url <- be_data_download_url(relpath)
    message(sprintf("File not found locally. Downloading to user cache... %s", url))
    message(sprintf("Cache target: %s", gene_exp_file))

    options(timeout = max(600, getOption("timeout")))
    tryCatch(
      {
        utils::download.file(url, gene_exp_file, method = "libcurl")
        message("File successfully downloaded.")
      },
      error = function(e) {
        stop(
          "Download failed for gene expression file. ",
          "Install BrainEnrichData or download the file manually from: ", url,
          " and place it at: ", gene_exp_file
        )
      }
    )
  }

  message("Reading gene expression data...")
  gene.df <- tryCatch(
    read.csv_bzip2(gene_exp_file, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      stop("Failed to read gene expression file: ", conditionMessage(e))
    }
  )

  if (hem %in% c("L", "R")) {
    gene.df <- dplyr::filter(gene.df, grepl(search_pattern, .data$Region))
  }

  gene.df <- gene.df %>%
    dplyr::filter(stats::complete.cases(dplyr::across(dplyr::everything()))) %>%
    tibble::column_to_rownames("Region")

  as.matrix(gene.df)
}
