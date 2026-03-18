#' Plot term-level enrichment results
#'
#' @param x An `EnrichRes` object.
#' @param type Plot type: `"dot"`, `"bar"`, `"volcano"`, or `"ridge"`.
#' @param top_n Number of top terms to display for non-volcano views.
#' @return A `ggplot2` object.
#' @export
plot_terms <- function(x, type = c("dot", "bar", "volcano", "ridge"), top_n = 20) {
  be_require_ggplot2()
  be_validate_enrichres(x)

  type <- match.arg(type)
  term_df <- be_term_plot_data(x)

  if (type == "volcano") {
    pvals <- if ("pvalue" %in% names(term_df)) term_df$pvalue else rep(NA_real_, nrow(term_df))
    term_df$neg_log10_p <- -log10(pmax(pvals, .Machine$double.xmin))
    term_df$significant <- if ("p.adjust" %in% names(term_df)) {
      !is.na(term_df$`p.adjust`) & term_df$`p.adjust` <= 0.05
    } else if ("pvalue" %in% names(term_df)) {
      !is.na(term_df$pvalue) & term_df$pvalue <= 0.05
    } else {
      FALSE
    }

    return(
      ggplot2::ggplot(term_df, ggplot2::aes(x = .data$effect, y = .data$neg_log10_p, color = .data$significant)) +
        ggplot2::geom_point(alpha = 0.8, size = 2.5) +
        ggplot2::scale_color_manual(values = c(`TRUE` = "#d95f02", `FALSE` = "#7570b3")) +
        ggplot2::labs(
          x = term_df$effect_label[[1]],
          y = expression(-log[10](pvalue)),
          color = "Significant",
          title = "Term Volcano Plot"
        ) +
        ggplot2::theme_minimal(base_size = 11)
    )
  }

  term_df <- be_top_terms(term_df, top_n = top_n)
  term_df$label <- stats::reorder(term_df$label, term_df$effect)

  if (type == "dot") {
    size_col <- if ("setSize" %in% names(term_df)) term_df$setSize else rep(3, nrow(term_df))
    color_col <- if ("pvalue" %in% names(term_df)) {
      -log10(pmax(term_df$pvalue, .Machine$double.xmin))
    } else {
      term_df$effect
    }
    term_df$plot_size <- size_col
    term_df$plot_color <- color_col

    return(
      ggplot2::ggplot(term_df, ggplot2::aes(x = .data$effect, y = .data$label)) +
        ggplot2::geom_point(ggplot2::aes(size = .data$plot_size, color = .data$plot_color), alpha = 0.9) +
        ggplot2::labs(
          x = term_df$effect_label[[1]],
          y = NULL,
          size = "Set size",
          color = if ("pvalue" %in% names(term_df)) expression(-log[10](pvalue)) else term_df$effect_label[[1]],
          title = "Top Enriched Terms"
        ) +
        ggplot2::theme_minimal(base_size = 11)
    )
  }

  if (type == "bar") {
    return(
      ggplot2::ggplot(term_df, ggplot2::aes(x = .data$effect, y = .data$label, fill = .data$effect)) +
        ggplot2::geom_col(show.legend = FALSE) +
        ggplot2::labs(
          x = term_df$effect_label[[1]],
          y = NULL,
          title = "Top Enriched Terms"
        ) +
        ggplot2::theme_minimal(base_size = 11)
    )
  }

  ridge_df <- be_ridge_data(x, top_n = top_n)
  ggplot2::ggplot(ridge_df$density_df, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$baseline,
        ymax = .data$baseline + .data$scaled_density,
        group = .data$term
      ),
      fill = "#9ecae1",
      alpha = 0.8
    ) +
    ggplot2::geom_segment(
      data = ridge_df$obs_df,
      ggplot2::aes(
        x = .data$observed,
        xend = .data$observed,
        y = .data$baseline,
        yend = .data$baseline + 0.85
      ),
      color = "#d95f02",
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_y_continuous(
      breaks = ridge_df$obs_df$baseline + 0.4,
      labels = ridge_df$obs_df$label
    ) +
    ggplot2::labs(
      x = ridge_df$effect_label,
      y = NULL,
      title = "Observed Scores Against Null Densities"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Plot core genes for a term
#'
#' @param x An `EnrichRes` object.
#' @param term_id Term identifier.
#' @param mode Plot either core-gene `"impact"` or role `"direction"`.
#' @return A `ggplot2` object.
#' @export
plot_core_genes <- function(x, term_id, mode = c("impact", "direction")) {
  be_require_ggplot2()
  be_validate_enrichres(x)

  mode <- match.arg(mode)
  core_df <- be_core_gene_details(x, term_id)

  if (nrow(core_df) == 0L) {
    stop("No core genes available for term: ", term_id)
  }

  if (mode == "impact") {
    core_df$gene_id <- stats::reorder(core_df$gene_id, core_df$impact)
    return(
      ggplot2::ggplot(core_df, ggplot2::aes(x = .data$impact, y = .data$gene_id, fill = .data$role)) +
        ggplot2::geom_col() +
        ggplot2::labs(
          x = "Leave-one-out impact",
          y = NULL,
          fill = "Role",
          title = paste("Core Genes for", term_id)
        ) +
        ggplot2::theme_minimal(base_size = 11)
    )
  }

  role_df <- as.data.frame(table(core_df$role), stringsAsFactors = FALSE)
  names(role_df) <- c("role", "count")
  ggplot2::ggplot(role_df, ggplot2::aes(x = .data$role, y = .data$count, fill = .data$role)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::labs(
      x = NULL,
      y = "Gene count",
      title = paste("Core Gene Directionality for", term_id)
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Plot a permutation heatmap for top terms
#'
#' @param x An `EnrichRes` object.
#' @param top_n Number of terms to display.
#' @param max_perm Maximum number of permutations to display.
#' @return A `ggplot2` object.
#' @export
plot_heatmap_terms <- function(x, top_n = 20, max_perm = 100) {
  be_require_ggplot2()
  be_validate_enrichres(x)

  perm_scores <- be_perm_score_matrix(x)
  term_df <- be_top_terms(be_term_plot_data(x), top_n = top_n)
  term_ids <- intersect(term_df$ID, rownames(perm_scores))

  if (length(term_ids) == 0L) {
    stop("No overlapping term IDs found between `x$table` and `x$perm_scores`.")
  }

  perm_scores <- perm_scores[term_ids, seq_len(min(ncol(perm_scores), max_perm)), drop = FALSE]
  heat_df <- as.data.frame(as.table(perm_scores), stringsAsFactors = FALSE)
  names(heat_df) <- c("term_id", "perm_idx", "score")
  heat_df$term_id <- factor(heat_df$term_id, levels = rev(term_ids))

  ggplot2::ggplot(heat_df, ggplot2::aes(x = .data$perm_idx, y = .data$term_id, fill = .data$score)) +
    ggplot2::geom_tile() +
    ggplot2::labs(
      x = "Permutation",
      y = NULL,
      fill = "Score",
      title = "Null Score Heatmap"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

#' Plot a lightweight term overlap network
#'
#' @param x An `EnrichRes` object.
#' @param top_n Number of terms to include.
#' @param min_overlap Minimum Jaccard overlap required for an edge.
#' @return A `ggplot2` object.
#' @export
plot_term_network <- function(x, top_n = 15, min_overlap = 0.1) {
  be_require_ggplot2()
  be_validate_enrichres(x)

  term_df <- be_top_terms(be_term_plot_data(x), top_n = top_n)
  gene_sets <- x$gene_sets[term_df$ID]
  gene_sets <- gene_sets[lengths(gene_sets) > 0L]

  if (length(gene_sets) == 0L) {
    stop("`x$gene_sets` does not contain non-empty gene sets for the selected terms.")
  }

  effect_label <- term_df$effect_label[[1]]
  node_ids <- names(gene_sets)
  angles <- seq(0, 2 * pi, length.out = length(node_ids) + 1L)[-1L]
  node_df <- data.frame(
    ID = node_ids,
    x = cos(angles),
    y = sin(angles),
    stringsAsFactors = FALSE
  )
  node_df <- merge(node_df, term_df[, c("ID", "label", "effect")], by = "ID", all.x = TRUE, sort = FALSE)

  edge_list <- vector("list", 0L)
  edge_idx <- 0L
  for (i in seq_len(length(gene_sets) - 1L)) {
    for (j in seq((i + 1L), length(gene_sets))) {
      gs_i <- unique(gene_sets[[i]])
      gs_j <- unique(gene_sets[[j]])
      overlap <- length(intersect(gs_i, gs_j)) / length(unique(c(gs_i, gs_j)))
      if (is.finite(overlap) && overlap >= min_overlap) {
        edge_idx <- edge_idx + 1L
        edge_list[[edge_idx]] <- data.frame(
          from = node_ids[[i]],
          to = node_ids[[j]],
          overlap = overlap,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  edge_df <- if (length(edge_list) > 0L) do.call(rbind, edge_list) else data.frame(from = character(), to = character(), overlap = numeric())
  if (nrow(edge_df) > 0L) {
    edge_df <- merge(edge_df, node_df[, c("ID", "x", "y")], by.x = "from", by.y = "ID", sort = FALSE)
    names(edge_df)[names(edge_df) %in% c("x", "y")] <- c("x_from", "y_from")
    edge_df <- merge(edge_df, node_df[, c("ID", "x", "y")], by.x = "to", by.y = "ID", sort = FALSE)
    names(edge_df)[names(edge_df) %in% c("x", "y")] <- c("x_to", "y_to")
  }

  p <- ggplot2::ggplot()
  if (nrow(edge_df) > 0L) {
    p <- p + ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(
        x = .data$x_from,
        y = .data$y_from,
        xend = .data$x_to,
        yend = .data$y_to,
        linewidth = .data$overlap,
        alpha = .data$overlap
      ),
      color = "#636363",
      show.legend = FALSE
    )
  }

  p +
    ggplot2::geom_point(
      data = node_df,
      ggplot2::aes(x = .data$x, y = .data$y, size = abs(.data$effect), color = .data$effect)
    ) +
    ggplot2::geom_text(
      data = node_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      vjust = -1,
      size = 3
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      color = effect_label,
      size = paste0("|", effect_label, "|"),
      title = "Term Overlap Network"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right")
}

be_validate_enrichres <- function(x) {
  if (!is_EnrichRes(x)) {
    stop("`x` must be an EnrichRes object.")
  }
  if (!is.data.frame(x$table) || nrow(x$table) == 0L) {
    stop("`x$table` must be a non-empty data frame.")
  }
}

be_require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plotting. Please install it first.")
  }
}

be_term_plot_data <- function(x) {
  tbl <- x$table
  effect_candidates <- c("gsScore", "Standardized_Coefficient", "Estimate", "t_Value", "Statistic")
  effect_col <- effect_candidates[effect_candidates %in% names(tbl)][1]
  if (is.na(effect_col) || is.null(effect_col)) {
    stop("Could not infer an effect column from `x$table`.")
  }

  out <- tbl
  out$ID <- as.character(out$ID)
  out$label <- if ("Description" %in% names(out)) {
    ifelse(is.na(out$Description) | out$Description == "", out$ID, out$Description)
  } else {
    out$ID
  }
  out$effect <- as.numeric(out[[effect_col]])
  out$effect_label <- effect_col
  out
}

be_top_terms <- function(term_df, top_n = 20) {
  top_n <- max(as.integer(top_n), 1L)
  order_idx <- if ("pvalue" %in% names(term_df)) {
    order(term_df$pvalue, -abs(term_df$effect), na.last = TRUE)
  } else {
    order(-abs(term_df$effect), na.last = TRUE)
  }
  term_df[utils::head(order_idx, top_n), , drop = FALSE]
}

be_perm_score_matrix <- function(x) {
  if (is.null(x$perm_scores)) {
    stop("`x$perm_scores` is required for this plot.")
  }
  perm_scores <- as.matrix(x$perm_scores)
  if (is.null(rownames(perm_scores))) {
    if ("ID" %in% names(x$table) && nrow(x$table) == nrow(perm_scores)) {
      rownames(perm_scores) <- x$table$ID
    } else {
      stop("`x$perm_scores` must have row names or match `x$table$ID` row count.")
    }
  }
  perm_scores
}

be_ridge_data <- function(x, top_n = 10) {
  perm_scores <- be_perm_score_matrix(x)
  term_df <- be_top_terms(be_term_plot_data(x), top_n = top_n)
  term_ids <- intersect(term_df$ID, rownames(perm_scores))

  if (length(term_ids) == 0L) {
    stop("No overlapping term IDs found between `x$table` and `x$perm_scores`.")
  }

  density_list <- vector("list", length(term_ids))
  obs_df <- data.frame(
    term = term_ids,
    label = term_df$label[match(term_ids, term_df$ID)],
    observed = term_df$effect[match(term_ids, term_df$ID)],
    baseline = seq_along(term_ids) - 1,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(term_ids)) {
    term_id <- term_ids[[i]]
    density_obj <- stats::density(as.numeric(perm_scores[term_id, ]), na.rm = TRUE)
    scaled_y <- density_obj$y / max(density_obj$y, na.rm = TRUE) * 0.8
    density_list[[i]] <- data.frame(
      term = term_id,
      x = density_obj$x,
      scaled_density = scaled_y,
      baseline = obs_df$baseline[[i]],
      stringsAsFactors = FALSE
    )
  }

  list(
    density_df = do.call(rbind, density_list),
    obs_df = obs_df,
    effect_label = term_df$effect_label[[1]]
  )
}

be_core_gene_details <- function(x, term_id) {
  labels <- NULL
  if (is.list(x$core_genes) && term_id %in% names(x$core_genes)) {
    labels <- x$core_genes[[term_id]]
  } else if ("core_enrichment" %in% names(x$table)) {
    row_idx <- match(term_id, x$table$ID)
    if (!is.na(row_idx)) {
      labels <- strsplit(as.character(x$table$core_enrichment[[row_idx]]), "/", fixed = TRUE)[[1]]
    }
  } else if ("core_genes" %in% names(x$table)) {
    row_idx <- match(term_id, x$table$ID)
    if (!is.na(row_idx)) {
      labels <- strsplit(as.character(x$table$core_genes[[row_idx]]), "/", fixed = TRUE)[[1]]
    }
  }

  if (is.null(labels) || length(labels) == 0L || identical(labels, "NOT_FOUND")) {
    return(data.frame(gene_id = character(), role = character(), impact = numeric()))
  }

  gene_id <- sub("\\(.*$", "", labels)
  role <- sub("^.*\\((.*)\\)$", "\\1", labels)
  role[role == labels] <- "unknown"
  impact <- attr(labels, "impact")
  if (is.null(impact) || length(impact) != length(labels)) {
    impact <- rev(seq_along(labels))
  }

  data.frame(
    gene_id = gene_id,
    role = role,
    impact = as.numeric(impact),
    stringsAsFactors = FALSE
  )
}
