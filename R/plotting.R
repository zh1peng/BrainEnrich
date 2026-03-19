#' Plot term-level enrichment results
#'
#' @param x An `EnrichRes` object.
#' @param type Plot type: `"dot"`, `"bar"`, `"volcano"`, or `"ridge"`.
#' @param top_n Number of top terms to display for non-volcano views.
#' @return A `ggplot2` object.
#' @export
plot_terms <- function(x, type = c("dot", "bar", "volcano", "ridge"), top_n = 20) {
  require_ggplot2()
  validate_enrichres(x)

  type <- match.arg(type)
  term_df <- term_plot_data(x)

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
      ggplot2::ggplot(term_df, ggplot2::aes(x = .data$effect, y = .data$neg_log10_p)) +
        ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#8c8c8c", linewidth = 0.35) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "#8c8c8c", linewidth = 0.35) +
        ggplot2::geom_point(
          ggplot2::aes(fill = .data$significant),
          shape = 21,
          color = "white",
          stroke = 0.35,
          alpha = 0.95,
          size = 3
        ) +
        ggplot2::scale_fill_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "#2c7fb8"), name = "FDR <= 0.05") +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.08))) +
        ggplot2::labs(
          x = term_df$effect_label[[1]],
          y = expression(-log[10](pvalue)),
          title = "Term Volcano Plot",
          subtitle = "Observed effect against permutation p-values"
        ) +
        plot_theme() +
        ggplot2::theme(legend.position = "top")
    )
  }

  term_df <- top_terms(term_df, top_n = top_n)
  term_df$label_display <- wrap_label(term_df$label, width = 34)
  term_df$label <- stats::reorder(term_df$label, term_df$effect)

  if (type == "dot") {
    size_col <- if ("setSize" %in% names(term_df)) term_df$setSize else rep(3, nrow(term_df))
    fill_col <- if ("pvalue" %in% names(term_df)) {
      -log10(pmax(term_df$pvalue, .Machine$double.xmin))
    } else {
      term_df$effect
    }
    term_df$plot_size <- size_col
    term_df$plot_fill <- fill_col
    fill_lab <- if ("pvalue" %in% names(term_df)) {
      expression(-log[10](pvalue))
    } else {
      term_df$effect_label[[1]]
    }
    fill_scale <- if ("pvalue" %in% names(term_df)) {
      ggplot2::scale_fill_gradient(low = "#eaf2fb", high = "#08306b")
    } else {
      ggplot2::scale_fill_gradient2(low = "#2c7fb8", mid = "#f7f7f7", high = "#c0392b", midpoint = 0)
    }

    return(
      ggplot2::ggplot(term_df, ggplot2::aes(x = .data$effect, y = .data$label_display)) +
        ggplot2::geom_point(
          ggplot2::aes(size = .data$plot_size, fill = .data$plot_fill),
          shape = 21,
          color = "white",
          stroke = 0.35,
          alpha = 0.95
        ) +
        fill_scale +
        ggplot2::scale_size_continuous(range = c(2.8, 8.5)) +
        ggplot2::labs(
          x = term_df$effect_label[[1]],
          y = NULL,
          size = "Set size",
          fill = fill_lab,
          title = "Top Enriched Terms",
          subtitle = "Larger points indicate larger gene sets"
        ) +
        plot_theme() +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
    )
  }

  if (type == "bar") {
    return(
      ggplot2::ggplot(term_df, ggplot2::aes(x = .data$effect, y = .data$label_display, fill = .data$effect)) +
        ggplot2::geom_col(width = 0.72, color = "white", linewidth = 0.25, show.legend = FALSE) +
        ggplot2::scale_fill_gradient2(low = "#2c7fb8", mid = "#f7f7f7", high = "#c0392b", midpoint = 0) +
        ggplot2::labs(
          x = term_df$effect_label[[1]],
          y = NULL,
          title = "Top Enriched Terms",
          subtitle = "Ranked by effect size"
        ) +
        plot_theme() +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
    )
  }

  ridge_df <- ridge_data(x, top_n = top_n)
  ggplot2::ggplot(ridge_df$density_df, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$baseline,
        ymax = .data$baseline + .data$scaled_density,
        group = .data$term
      ),
      fill = "#2c7fb8",
      alpha = 0.28
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = .data$baseline + .data$scaled_density,
        group = .data$term
      ),
      color = "#08306b",
      linewidth = 0.35
    ) +
    ggplot2::geom_segment(
      data = ridge_df$obs_df,
      ggplot2::aes(
        x = .data$observed,
        xend = .data$observed,
        y = .data$baseline,
        yend = .data$baseline + 0.85
      ),
      color = "#c0392b",
      linewidth = 0.65,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_y_continuous(
      breaks = ridge_df$obs_df$baseline + 0.4,
      labels = wrap_label(ridge_df$obs_df$label, width = 30)
    ) +
    ggplot2::labs(
      x = ridge_df$effect_label,
      y = NULL,
      title = "Observed Scores Against Null Densities",
      subtitle = paste0("Top ", nrow(ridge_df$obs_df), " terms and the first ", ncol(perm_score_matrix(x)), " permutations")
    ) +
    plot_theme() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.y = ggplot2::element_text(size = 9)
    )
}

#' Plot core genes for a term
#'
#' @param x An `EnrichRes` object.
#' @param term_id Term identifier.
#' @param mode Plot either core-gene `"impact"` or role `"direction"`.
#' @return A `ggplot2` object.
#' @export
plot_core_genes <- function(x, term_id, mode = c("impact", "direction")) {
  require_ggplot2()
  validate_enrichres(x)

  mode <- match.arg(mode)
  core_df <- core_gene_details(x, term_id)

  if (nrow(core_df) == 0L) {
    stop("No core genes available for term: ", term_id)
  }

  if (mode == "impact") {
    core_df$gene_id <- stats::reorder(core_df$gene_id, core_df$impact)
    return(
      ggplot2::ggplot(core_df, ggplot2::aes(x = .data$impact, y = .data$gene_id, fill = .data$role)) +
        ggplot2::geom_col(width = 0.72, color = "white", linewidth = 0.25) +
        ggplot2::scale_fill_manual(values = role_palette(), drop = FALSE) +
        ggplot2::labs(
          x = "Leave-one-out impact",
          y = NULL,
          fill = "Role",
          title = paste("Core Genes for", term_id),
          subtitle = "Driver genes reinforce the observed statistic; buffer genes oppose it"
        ) +
        plot_theme(legend_position = "bottom") +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
    )
  }

  role_df <- as.data.frame(table(core_df$role), stringsAsFactors = FALSE)
  names(role_df) <- c("role", "count")
  role_df$role <- factor(role_df$role, levels = names(role_palette()))

  ggplot2::ggplot(role_df, ggplot2::aes(x = .data$role, y = .data$count, fill = .data$role)) +
    ggplot2::geom_col(width = 0.68, color = "white", linewidth = 0.25, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = .data$count), vjust = -0.45, fontface = "bold", size = 3.5) +
    ggplot2::scale_fill_manual(values = role_palette(), drop = FALSE) +
    ggplot2::labs(
      x = NULL,
      y = "Gene count",
      title = paste("Core Gene Directionality for", term_id),
      subtitle = "Counts of directional core-gene roles"
    ) +
    plot_theme(legend_position = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
}

#' Plot a permutation heatmap for top terms
#'
#' @param x An `EnrichRes` object.
#' @param top_n Number of terms to display.
#' @param max_perm Maximum number of permutations to display.
#' @return A `ggplot2` object.
#' @export
plot_heatmap_terms <- function(x, top_n = 20, max_perm = 100) {
  require_ggplot2()
  validate_enrichres(x)

  perm_scores <- perm_score_matrix(x)
  term_df <- top_terms(term_plot_data(x), top_n = top_n)
  term_ids <- intersect(term_df$ID, rownames(perm_scores))

  if (length(term_ids) == 0L) {
    stop("No overlapping term IDs found between `x$table` and `x$perm_scores`.")
  }

  perm_scores <- perm_scores[term_ids, seq_len(min(ncol(perm_scores), max_perm)), drop = FALSE]
  heat_df <- as.data.frame(as.table(perm_scores), stringsAsFactors = FALSE)
  names(heat_df) <- c("term_id", "perm_idx", "score")
  heat_df$perm_idx <- as.integer(as.character(heat_df$perm_idx))
  term_labels <- wrap_label(term_df$label[match(term_ids, term_df$ID)], width = 28)
  heat_df$term_label <- factor(term_labels[match(heat_df$term_id, term_ids)], levels = rev(unique(term_labels)))

  ggplot2::ggplot(heat_df, ggplot2::aes(x = .data$perm_idx, y = .data$term_label, fill = .data$score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.12) +
    ggplot2::scale_fill_gradient2(low = "#2c7fb8", mid = "#f7f7f7", high = "#c0392b", midpoint = 0, name = "Score") +
    ggplot2::labs(
      x = "Permutation",
      y = NULL,
      title = "Null Score Heatmap",
      subtitle = paste0("Top ", length(term_ids), " terms across ", ncol(perm_scores), " permutations")
    ) +
    plot_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 9)
    )
}

#' Plot a lightweight term overlap network
#'
#' @param x An `EnrichRes` object.
#' @param top_n Number of terms to include.
#' @param min_overlap Minimum Jaccard overlap required for an edge.
#' @return A `ggplot2` object.
#' @export
plot_term_network <- function(x, top_n = 15, min_overlap = 0.1) {
  require_ggplot2()
  validate_enrichres(x)

  term_df <- top_terms(term_plot_data(x), top_n = top_n)
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
  node_df$label_display <- wrap_label(node_df$label, width = 18)
  node_df$label_x <- node_df$x * 1.12
  node_df$label_y <- node_df$y * 1.12

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
      color = "#7a7a7a",
      show.legend = FALSE
    )
  }

  p +
    ggplot2::geom_point(
      data = node_df,
      ggplot2::aes(x = .data$x, y = .data$y, size = abs(.data$effect), fill = .data$effect),
      shape = 21,
      color = "white",
      stroke = 0.65
    ) +
    ggplot2::geom_label(
      data = node_df,
      ggplot2::aes(x = .data$label_x, y = .data$label_y, label = .data$label_display),
      size = 3,
      label.size = 0.2,
      label.padding = grid::unit(0.12, "lines"),
      fill = "#ffffffe6",
      color = "#2b2b2b",
      fontface = "bold"
    ) +
    ggplot2::scale_fill_gradient2(low = "#2c7fb8", mid = "#f7f7f7", high = "#c0392b", midpoint = 0) +
    ggplot2::scale_size_continuous(range = c(4, 11)) +
    ggplot2::coord_equal(xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3), clip = "off") +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      size = paste0("|", effect_label, "|"),
      fill = effect_label,
      title = "Term Overlap Network",
      subtitle = "Node size reflects effect magnitude; edge thickness reflects Jaccard overlap"
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, margin = ggplot2::margin(b = 8)),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

plot_theme <- function(base_size = 12, legend_position = "right") {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "#1f1f1f"),
      plot.title = ggplot2::element_text(face = "bold", margin = ggplot2::margin(b = 6)),
      plot.subtitle = ggplot2::element_text(color = "#4f4f4f", margin = ggplot2::margin(b = 8)),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "#2b2b2b"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = legend_position,
      panel.grid.major = ggplot2::element_line(color = "#e6e6e6", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 12, 10, 12)
    )
}

wrap_label <- function(x, width = 32) {
  vapply(
    x,
    function(label) {
      if (is.na(label) || !nzchar(label)) {
        return("")
      }
      paste(strwrap(as.character(label), width = width), collapse = "\n")
    },
    character(1),
    USE.NAMES = FALSE
  )
}

role_palette <- function() {
  c(
    driver = "#c0392b",
    buffer = "#2c7fb8",
    neutral = "#7f7f7f",
    unknown = "#7f7f7f"
  )
}

validate_enrichres <- function(x) {
  if (!is_EnrichRes(x)) {
    stop("`x` must be an EnrichRes object.")
  }
  if (!is.data.frame(x$table) || nrow(x$table) == 0L) {
    stop("`x$table` must be a non-empty data frame.")
  }
}

require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plotting. Please install it first.")
  }
}

term_plot_data <- function(x) {
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

top_terms <- function(term_df, top_n = 20) {
  top_n <- max(as.integer(top_n), 1L)
  order_idx <- if ("pvalue" %in% names(term_df)) {
    order(term_df$pvalue, -abs(term_df$effect), na.last = TRUE)
  } else {
    order(-abs(term_df$effect), na.last = TRUE)
  }
  term_df[utils::head(order_idx, top_n), , drop = FALSE]
}

perm_score_matrix <- function(x) {
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

ridge_data <- function(x, top_n = 10) {
  perm_scores <- perm_score_matrix(x)
  term_df <- top_terms(term_plot_data(x), top_n = top_n)
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

core_gene_details <- function(x, term_id) {
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
