#' Construct an EnrichAnno object
#'
#' @param term2gene Data frame with columns `term_id` and `gene_id`.
#' @param term2name Data frame with columns `term_id` and `term_name`.
#' @param meta Optional named list of metadata.
#' @return An object of class `EnrichAnno`.
new_EnrichAnno <- function(term2gene, term2name, meta = list()) {
  if (!is.data.frame(term2gene) || !all(c("term_id", "gene_id") %in% names(term2gene))) {
    stop("term2gene must be a data frame with columns: term_id, gene_id.")
  }
  if (!is.data.frame(term2name) || !all(c("term_id", "term_name") %in% names(term2name))) {
    stop("term2name must be a data frame with columns: term_id, term_name.")
  }

  term2gene <- term2gene[, c("term_id", "gene_id"), drop = FALSE]
  term2name <- term2name[, c("term_id", "term_name"), drop = FALSE]
  term2gene <- unique(term2gene)
  term2name <- unique(term2name)

  structure(
    list(
      term2gene = term2gene,
      term2name = term2name,
      meta = meta
    ),
    class = "EnrichAnno"
  )
}

#' Check if object is EnrichAnno
#'
#' @param x Object.
#' @return Logical scalar.
is_EnrichAnno <- function(x) {
  is.list(x) &&
    inherits(x, "EnrichAnno") &&
    all(c("term2gene", "term2name") %in% names(x))
}

#' Construct an EnrichRes object
#'
#' @param analysis_type One of `brainenrich` or `brainscore_lm_test`.
#' @param table Data frame of term-level results.
#' @param gene_sets Named list of gene sets.
#' @param perm_scores Matrix of null/permutation scores.
#' @param params Named list of run parameters.
#' @param diagnostics Named list of diagnostics.
#' @param model_table Optional data frame for model-level results.
#' @param np_table Optional data frame for non-parametric summaries.
#' @param core_genes Optional data frame/list with core genes.
#' @param gene_scores Optional matrix/vector with gene-level scores.
#' @return An object of class `EnrichRes`.
new_EnrichRes <- function(analysis_type,
                          table,
                          gene_sets,
                          perm_scores = NULL,
                          params = list(),
                          diagnostics = list(),
                          model_table = NULL,
                          np_table = NULL,
                          core_genes = NULL,
                          gene_scores = NULL) {
  if (!analysis_type %in% c("brainenrich", "brainscore_lm_test")) {
    stop("analysis_type must be one of: brainenrich, brainscore_lm_test.")
  }
  if (!is.data.frame(table)) {
    stop("table must be a data frame.")
  }
  if (!is.list(gene_sets)) {
    stop("gene_sets must be a list.")
  }

  structure(
    list(
      analysis_type = analysis_type,
      table = table,
      gene_sets = gene_sets,
      perm_scores = perm_scores,
      params = params,
      diagnostics = diagnostics,
      model_table = model_table,
      np_table = np_table,
      core_genes = core_genes,
      gene_scores = gene_scores
    ),
    class = "EnrichRes"
  )
}

#' Check if object is EnrichRes
#'
#' @param x Object.
#' @return Logical scalar.
is_EnrichRes <- function(x) {
  is.list(x) &&
    inherits(x, "EnrichRes") &&
    all(c("analysis_type", "table", "gene_sets") %in% names(x))
}

#' @export
print.EnrichRes <- function(x, ...) {
  n_terms <- if (is.data.frame(x$table)) nrow(x$table) else NA_integer_
  cat("<EnrichRes>\n")
  cat("  analysis_type:", x$analysis_type, "\n")
  cat("  terms:", n_terms, "\n")
  cat("  gene_sets:", length(x$gene_sets), "\n")
  invisible(x)
}

#' @export
as.data.frame.EnrichRes <- function(x, ...) {
  x$table
}

#' Internal q-value helper aligned with the legacy DOSE behavior.
#'
#' @param pvals Numeric vector of p-values.
#' @return Numeric vector of q-values.
calculate_qvalue_local <- function(pvals) {
  if (length(pvals) == 0L) {
    return(numeric(0))
  }

  qobj <- tryCatch(
    qvalue::qvalue(pvals, lambda = 0.05, pi0.method = "bootstrap"),
    error = function(e) NULL
  )

  if (inherits(qobj, "qvalue")) {
    qobj$qvalues
  } else {
    NA_real_
  }
}

#' Internal helper to convert a legacy annotation environment to EnrichAnno.
#'
#' @param anno_env Environment with PATHID2EXTID/PATHID2NAME.
#' @return `EnrichAnno`.
anno_env_to_EnrichAnno <- function(anno_env) {
  if (!is.environment(anno_env)) {
    stop("anno_env must be an environment.")
  }
  if (!exists("PATHID2EXTID", envir = anno_env, inherits = TRUE)) {
    stop("Annotation environment missing PATHID2EXTID.")
  }

  path2extid <- get("PATHID2EXTID", envir = anno_env, inherits = TRUE)
  path2name <- if (exists("PATHID2NAME", envir = anno_env, inherits = TRUE)) {
    get("PATHID2NAME", envir = anno_env, inherits = TRUE)
  } else {
    NULL
  }

  term2gene <- utils::stack(path2extid)
  term2gene <- term2gene[, c("ind", "values")]
  names(term2gene) <- c("term_id", "gene_id")

  if (!is.null(path2name)) {
    term2name <- utils::stack(path2name)
    term2name <- term2name[, c("ind", "values")]
    names(term2name) <- c("term_id", "term_name")
  } else {
    uniq_terms <- unique(term2gene$term_id)
    term2name <- data.frame(term_id = uniq_terms, term_name = uniq_terms, stringsAsFactors = FALSE)
  }

  new_EnrichAnno(term2gene = term2gene, term2name = term2name, meta = list(source = "legacy_env"))
}

#' Internal helper to build a named gene-set list from EnrichAnno.
#'
#' @param annoData `EnrichAnno`.
#' @return Named list of gene vectors.
anno_to_geneSetList <- function(annoData) {
  if (!is_EnrichAnno(annoData)) {
    stop("annoData must be an EnrichAnno object.")
  }
  split(annoData$term2gene$gene_id, annoData$term2gene$term_id)
}

#' Internal helper to map term id to description.
#'
#' @param term Character vector of term IDs.
#' @param annoData `EnrichAnno`.
#' @return Character vector.
anno_term_description <- function(term, annoData) {
  if (!is_EnrichAnno(annoData)) {
    stop("annoData must be an EnrichAnno object.")
  }
  idx <- match(term, annoData$term2name$term_id)
  annoData$term2name$term_name[idx]
}
