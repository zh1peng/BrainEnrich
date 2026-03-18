#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
baseline_path <- if (length(args) >= 1L) args[[1]] else file.path(getwd(), "reports", "parity_baseline.rds")
refactor_path <- if (length(args) >= 2L) args[[2]] else file.path(getwd(), "reports", "parity_refactor.rds")
output_rds <- if (length(args) >= 3L) args[[3]] else file.path(getwd(), "reports", "parity_compare.rds")
output_md <- if (length(args) >= 4L) args[[4]] else file.path(getwd(), "reports", "parity_summary.md")
tolerance <- if (length(args) >= 5L) as.numeric(args[[5]]) else 1e-6

is_gsea_result <- function(x) {
  methods::is(x, "gseaResult")
}

is_enrich_res <- function(x) {
  inherits(x, "EnrichRes")
}

normalize_scalar <- function(x) {
  if (is.factor(x)) {
    return(as.character(x))
  }
  x
}

normalize_core_string <- function(x) {
  x <- normalize_scalar(x)
  vapply(x, function(value) {
    if (is.na(value) || !nzchar(value)) {
      return("")
    }
    genes <- unlist(strsplit(as.character(value), "/", fixed = TRUE), use.names = FALSE)
    genes <- gsub("\\s*\\([^)]*\\)$", "", genes)
    genes <- sort(unique(genes[nzchar(genes)]))
    paste(genes, collapse = "/")
  }, FUN.VALUE = character(1))
}

canonicalize_table <- function(df) {
  if (is.null(df)) {
    return(data.frame(ID = character(), stringsAsFactors = FALSE))
  }

  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  df[] <- lapply(df, normalize_scalar)

  if ("Dependent_vars" %in% names(df) && !"ID" %in% names(df)) {
    names(df)[names(df) == "Dependent_vars"] <- "ID"
  }
  if (!"ID" %in% names(df)) {
    if (nrow(df) > 0L && !is.null(rownames(df))) {
      df$ID <- rownames(df)
    } else {
      df$ID <- character(nrow(df))
    }
  }

  if ("core_enrichment" %in% names(df)) {
    df$core_enrichment <- normalize_core_string(df$core_enrichment)
  }
  if ("core_genes" %in% names(df)) {
    df$core_genes <- normalize_core_string(df$core_genes)
  }

  if (nrow(df) > 0L) {
    df <- df[order(df$ID), , drop = FALSE]
  }
  rownames(df) <- NULL
  df
}

extract_enrich_table <- function(x) {
  if (is_gsea_result(x)) {
    return(methods::slot(x, "result"))
  }
  if (is_enrich_res(x)) {
    return(x$table)
  }
  if (is.data.frame(x)) {
    return(x)
  }
  NULL
}

extract_enrich_gene_sets <- function(x) {
  if (is_gsea_result(x)) {
    return(methods::slot(x, "geneSets"))
  }
  if (is_enrich_res(x)) {
    return(x$gene_sets)
  }
  NULL
}

extract_enrich_perm_scores <- function(x) {
  if (is_gsea_result(x)) {
    return(methods::slot(x, "permScores"))
  }
  if (is_enrich_res(x)) {
    return(x$perm_scores)
  }
  NULL
}

normalize_core_list <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(stats::setNames(character(), character()))
  }
  vals <- vapply(x, function(item) {
    item <- as.character(item)
    item <- gsub("\\s*\\([^)]*\\)$", "", item)
    item <- sort(unique(item[nzchar(item)]))
    paste(item, collapse = "/")
  }, FUN.VALUE = character(1))
  vals[order(names(vals))]
}

extract_core_map <- function(x) {
  table <- canonicalize_table(extract_enrich_table(x))
  if ("core_enrichment" %in% names(table)) {
    out <- stats::setNames(table$core_enrichment, table$ID)
    return(out[order(names(out))])
  }
  if ("core_genes" %in% names(table)) {
    out <- stats::setNames(table$core_genes, table$ID)
    return(out[order(names(out))])
  }
  if (is_enrich_res(x) && is.list(x$core_genes)) {
    return(normalize_core_list(x$core_genes))
  }
  stats::setNames(character(), character())
}

canonicalize_gene_sets <- function(x) {
  if (is.null(x) || !is.list(x)) {
    return(NULL)
  }
  x <- x[order(names(x))]
  lapply(x, function(genes) sort(unique(as.character(genes))))
}

canonicalize_perm_scores <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  m <- as.matrix(x)
  list(
    dim = dim(m),
    colnames = colnames(m),
    values = unname(m)
  )
}

compare_gene_sets <- function(a, b) {
  if (is.null(a) && is.null(b)) {
    return(list(ok = TRUE, details = "OK"))
  }
  if (is.null(a) || is.null(b)) {
    return(list(ok = FALSE, details = "gene sets missing on one side"))
  }
  if (!identical(names(a), names(b))) {
    return(list(ok = FALSE, details = "gene set names differ"))
  }
  mismatch <- names(a)[!mapply(identical, a, b)]
  if (length(mismatch) > 0L) {
    return(list(ok = FALSE, details = paste("gene membership differs for:", paste(mismatch, collapse = ", "))))
  }
  list(ok = TRUE, details = "OK")
}

compare_core_maps <- function(a, b) {
  if (!identical(names(a), names(b))) {
    return(list(ok = FALSE, details = "core gene term IDs differ"))
  }
  if (!identical(unname(a), unname(b))) {
    return(list(ok = FALSE, details = "core gene membership differs"))
  }
  list(ok = TRUE, details = "OK")
}

compare_perm_scores <- function(a, b) {
  if (is.null(a) && is.null(b)) {
    return(list(ok = TRUE, dim_match = TRUE, numeric_match = TRUE, details = "OK"))
  }
  if (is.null(a) || is.null(b)) {
    return(list(ok = FALSE, dim_match = FALSE, numeric_match = FALSE, details = "permScores missing on one side"))
  }
  dim_match <- identical(a$dim, b$dim)
  col_match <- TRUE
  if (!is.null(a$colnames) && !is.null(b$colnames)) {
    col_match <- identical(a$colnames, b$colnames)
  }
  numeric_match <- isTRUE(all.equal(a$values, b$values, tolerance = tolerance, check.attributes = FALSE))
  details <- c()
  if (!dim_match) details <- c(details, "permScores dimensions differ")
  if (!col_match) details <- c(details, "permScores columns differ")
  if (!numeric_match) details <- c(details, "permScores differ beyond tolerance")
  list(
    ok = dim_match && col_match && numeric_match,
    dim_match = dim_match,
    numeric_match = numeric_match,
    details = if (length(details) > 0L) paste(details, collapse = "; ") else "OK"
  )
}

compare_tables <- function(a, b) {
  a <- canonicalize_table(a)
  b <- canonicalize_table(b)

  column_match <- identical(names(a), names(b))
  id_match <- identical(a$ID, b$ID)
  dim_match <- identical(dim(a), dim(b))

  char_cols_a <- names(a)[!vapply(a, is.numeric, FUN.VALUE = logical(1))]
  char_cols_b <- names(b)[!vapply(b, is.numeric, FUN.VALUE = logical(1))]
  numeric_cols_a <- names(a)[vapply(a, is.numeric, FUN.VALUE = logical(1))]
  numeric_cols_b <- names(b)[vapply(b, is.numeric, FUN.VALUE = logical(1))]

  char_col_match <- identical(char_cols_a, char_cols_b)
  numeric_col_match <- identical(numeric_cols_a, numeric_cols_b)

  char_match <- if (length(char_cols_a) > 0L && length(char_cols_b) > 0L) {
    isTRUE(all.equal(a[, char_cols_a, drop = FALSE], b[, char_cols_b, drop = FALSE], check.attributes = FALSE))
  } else {
    TRUE
  }
  numeric_match <- if (length(numeric_cols_a) > 0L && length(numeric_cols_b) > 0L) {
    isTRUE(all.equal(a[, numeric_cols_a, drop = FALSE], b[, numeric_cols_b, drop = FALSE], tolerance = tolerance, check.attributes = FALSE))
  } else {
    TRUE
  }

  details <- c()
  if (!column_match) details <- c(details, "table columns differ")
  if (!id_match) details <- c(details, "table IDs differ")
  if (!dim_match) details <- c(details, "table dimensions differ")
  if (!char_col_match) details <- c(details, "table character columns differ")
  if (!numeric_col_match) details <- c(details, "table numeric columns differ")
  if (!char_match) details <- c(details, "table text payload differs")
  if (!numeric_match) details <- c(details, "table numeric payload differs beyond tolerance")

  list(
    pass = column_match && id_match && dim_match && char_col_match && numeric_col_match && char_match && numeric_match,
    key_match = column_match && id_match && char_col_match && char_match,
    dim_match = dim_match,
    numeric_match = numeric_col_match && numeric_match,
    details = if (length(details) > 0L) paste(details, collapse = "; ") else "OK"
  )
}

normalize_nested <- function(x) {
  if (is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
    x[] <- lapply(x, normalize_scalar)
    return(x)
  }
  if (is.matrix(x) || is.array(x)) {
    return(x)
  }
  if (is.list(x)) {
    out <- lapply(x, normalize_nested)
    names(out) <- names(x)
    return(out)
  }
  normalize_scalar(x)
}

collect_relevant_attrs <- function(x) {
  attrs <- attributes(x)
  keep <- c("null_model", "cor_method", "aggre_method", "minGSSize", "maxGSSize", "n_perm")
  out <- attrs[intersect(names(attrs), keep)]
  if (length(out) == 0L) {
    return(list())
  }
  out
}

extract_key_signature <- function(x) {
  if (is.data.frame(x)) {
    return(list(type = "data.frame", names = names(x), rownames = rownames(x)))
  }
  if (is.matrix(x)) {
    return(list(type = "matrix", rownames = rownames(x), colnames = colnames(x)))
  }
  if (is.list(x)) {
    out <- lapply(x, extract_key_signature)
    names(out) <- names(x)
    return(out)
  }
  class(x)
}

extract_dim_signature <- function(x) {
  if (is.data.frame(x) || is.matrix(x) || is.array(x)) {
    return(dim(x))
  }
  if (is.list(x)) {
    out <- lapply(x, extract_dim_signature)
    names(out) <- names(x)
    return(out)
  }
  length(x)
}

compare_generic <- function(a, b) {
  na <- list(attrs = collect_relevant_attrs(a), payload = normalize_nested(a))
  nb <- list(attrs = collect_relevant_attrs(b), payload = normalize_nested(b))

  key_match <- identical(extract_key_signature(na), extract_key_signature(nb))
  dim_match <- identical(extract_dim_signature(na), extract_dim_signature(nb))
  numeric_match <- isTRUE(all.equal(na, nb, tolerance = tolerance, check.attributes = FALSE))

  details <- c()
  if (!key_match) details <- c(details, "keys differ")
  if (!dim_match) details <- c(details, "dimensions differ")
  if (!numeric_match) details <- c(details, "payload differs beyond tolerance")

  list(
    key_match = key_match,
    dim_match = dim_match,
    gene_set_match = NA,
    core_gene_match = NA,
    numeric_match = numeric_match,
    pass = key_match && dim_match && numeric_match,
    details = if (length(details) > 0L) paste(details, collapse = "; ") else "OK"
  )
}

compare_enrich_like <- function(a, b) {
  table_cmp <- compare_tables(extract_enrich_table(a), extract_enrich_table(b))
  gene_set_cmp <- compare_gene_sets(
    canonicalize_gene_sets(extract_enrich_gene_sets(a)),
    canonicalize_gene_sets(extract_enrich_gene_sets(b))
  )
  core_cmp <- compare_core_maps(extract_core_map(a), extract_core_map(b))
  perm_cmp <- compare_perm_scores(
    canonicalize_perm_scores(extract_enrich_perm_scores(a)),
    canonicalize_perm_scores(extract_enrich_perm_scores(b))
  )

  details <- c(table_cmp$details)
  if (!gene_set_cmp$ok) details <- c(details, gene_set_cmp$details)
  if (!core_cmp$ok) details <- c(details, core_cmp$details)
  if (!perm_cmp$ok) details <- c(details, perm_cmp$details)
  details <- unique(details[details != "OK"])

  list(
    key_match = table_cmp$key_match,
    dim_match = table_cmp$dim_match && perm_cmp$dim_match,
    gene_set_match = gene_set_cmp$ok,
    core_gene_match = core_cmp$ok,
    numeric_match = table_cmp$numeric_match && perm_cmp$numeric_match,
    pass = table_cmp$pass && gene_set_cmp$ok && core_cmp$ok && perm_cmp$ok,
    details = if (length(details) > 0L) paste(details, collapse = "; ") else "OK"
  )
}

compare_table_only <- function(a, b) {
  table_cmp <- compare_tables(a, b)
  core_cmp <- compare_core_maps(extract_core_map(a), extract_core_map(b))
  details <- c(table_cmp$details)
  if (!core_cmp$ok) details <- c(details, core_cmp$details)
  details <- unique(details[details != "OK"])

  list(
    key_match = table_cmp$key_match,
    dim_match = table_cmp$dim_match,
    gene_set_match = NA,
    core_gene_match = core_cmp$ok,
    numeric_match = table_cmp$numeric_match,
    pass = table_cmp$pass && core_cmp$ok,
    details = if (length(details) > 0L) paste(details, collapse = "; ") else "OK"
  )
}

read_payload <- function(path) {
  payload <- readRDS(path)
  if (!"scenarios" %in% names(payload) && "results" %in% names(payload)) {
    payload$scenarios <- lapply(names(payload$results), function(name) {
      list(workflow = name, tier = "legacy", result = payload$results[[name]])
    })
    names(payload$scenarios) <- names(payload$results)
  }
  payload
}

baseline <- read_payload(baseline_path)
refactor <- read_payload(refactor_path)
scenario_names <- union(names(baseline$scenarios), names(refactor$scenarios))

comparisons <- lapply(scenario_names, function(scenario_name) {
  a <- baseline$scenarios[[scenario_name]]
  b <- refactor$scenarios[[scenario_name]]

  if (is.null(a) || is.null(b)) {
    return(data.frame(
      scenario = scenario_name,
      workflow = if (!is.null(a)) a$workflow else if (!is.null(b)) b$workflow else NA_character_,
      tier = if (!is.null(a)) a$tier else if (!is.null(b)) b$tier else NA_character_,
      status = "FAIL",
      key_match = FALSE,
      dim_match = FALSE,
      gene_set_match = NA,
      core_gene_match = NA,
      numeric_match = FALSE,
      details = "scenario missing on one side",
      stringsAsFactors = FALSE
    ))
  }

  result_cmp <- if ((is_gsea_result(a$result) || is_enrich_res(a$result)) ||
    (is_gsea_result(b$result) || is_enrich_res(b$result))) {
    compare_enrich_like(a$result, b$result)
  } else if (is.data.frame(a$result) && is.data.frame(b$result)) {
    compare_table_only(a$result, b$result)
  } else {
    compare_generic(a$result, b$result)
  }

  data.frame(
    scenario = scenario_name,
    workflow = b$workflow,
    tier = b$tier,
    status = if (result_cmp$pass) "PASS" else "FAIL",
    key_match = result_cmp$key_match,
    dim_match = result_cmp$dim_match,
    gene_set_match = ifelse(is.na(result_cmp$gene_set_match), NA, result_cmp$gene_set_match),
    core_gene_match = ifelse(is.na(result_cmp$core_gene_match), NA, result_cmp$core_gene_match),
    numeric_match = result_cmp$numeric_match,
    details = result_cmp$details,
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, comparisons)
summary_df <- summary_df[, c("scenario", "workflow", "tier", "status", "key_match", "dim_match", "gene_set_match", "core_gene_match", "numeric_match", "details")]

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
payload <- list(
  baseline = baseline_path,
  refactor = refactor_path,
  tolerance = tolerance,
  summary = summary_df
)
saveRDS(payload, output_rds)

md_lines <- c(
  "# Parity Summary",
  "",
  paste("- Baseline:", baseline_path),
  paste("- Refactor:", refactor_path),
  paste("- Tolerance:", tolerance),
  "",
  "| Scenario | Workflow | Tier | Status | Keys | Dims | Gene Sets | Core Genes | Numeric | Details |",
  "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
)

for (i in seq_len(nrow(summary_df))) {
  md_lines <- c(
    md_lines,
    sprintf(
      "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
      summary_df$scenario[[i]],
      summary_df$workflow[[i]],
      summary_df$tier[[i]],
      summary_df$status[[i]],
      summary_df$key_match[[i]],
      summary_df$dim_match[[i]],
      ifelse(is.na(summary_df$gene_set_match[[i]]), "NA", summary_df$gene_set_match[[i]]),
      ifelse(is.na(summary_df$core_gene_match[[i]]), "NA", summary_df$core_gene_match[[i]]),
      summary_df$numeric_match[[i]],
      summary_df$details[[i]]
    )
  )
}

writeLines(md_lines, output_md)
message("Saved comparison report to ", output_rds)
message("Saved markdown summary to ", output_md)
