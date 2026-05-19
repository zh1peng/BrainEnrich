#' Assess Normality of Individual Gene Set Scores
#'
#' This function reports normality diagnostics for individual-level gene set
#' scores. The diagnostics are intended to guide downstream model choice and are
#' not used as automatic exclusion criteria.
#'
#' @param gsScore A data frame, matrix, or empirical `brainscore()` output list
#'   containing individual gene set scores, with subjects/participants as rows
#'   and gene sets as columns.
#' @param method Character string specifying the normality diagnostic method.
#'   Default is `"ks"`. Other options are `"shapiro"` and `"both"`.
#' @param alpha Numeric significance threshold used to flag non-normal score
#'   distributions after p-value adjustment. Default is 0.05.
#' @param p_adjust_method Character string specifying the method for p-value
#'   adjustment. Default is `"fdr"`. See [stats::p.adjust()] for details.
#' @param shapiro_max_n Maximum sample size used for Shapiro-Wilk tests.
#'   `stats::shapiro.test` supports at most 5000 observations. Default is 5000.
#' @param seed Optional random seed used when subsampling observations for
#'   Shapiro-Wilk tests with more than `shapiro_max_n` observations.
#'
#' @return A data frame with one row per gene set and normality diagnostic
#'   columns.
#' @export
brainscore.normality <- function(gsScore,
                                 method = c("ks", "shapiro", "both"),
                                 alpha = 0.05,
                                 p_adjust_method = "fdr",
                                 shapiro_max_n = 5000,
                                 seed = NULL) {
  method <- match.arg(method)

  if (is.list(gsScore) && !is.data.frame(gsScore)) {
    if (!all(vapply(gsScore, is.numeric, logical(1)))) {
      stop("gsScore list columns must be numeric vectors.")
    }
    if (length(unique(vapply(gsScore, length, integer(1)))) != 1) {
      stop("gsScore list columns must have the same length.")
    }
    gsScore <- data.frame(gsScore, check.names = FALSE)
  }

  if (!is.data.frame(gsScore) && !is.matrix(gsScore)) {
    stop("gsScore must be a data frame, matrix, or empirical brainscore output list.")
  }
  if (is.null(colnames(gsScore))) {
    stop("gsScore must have column names for gene set IDs.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(shapiro_max_n) || length(shapiro_max_n) != 1 || is.na(shapiro_max_n) || shapiro_max_n < 3 || shapiro_max_n > 5000) {
    stop("shapiro_max_n must be a single numeric value between 3 and 5000.")
  }

  shapiro_max_n <- as.integer(shapiro_max_n)
  gsScore <- as.data.frame(gsScore, check.names = FALSE)
  if (!all(vapply(gsScore, is.numeric, logical(1)))) {
    stop("All gsScore columns must be numeric.")
  }

  if (!is.null(seed)) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (has_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  diagnostics <- lapply(colnames(gsScore), function(id) {
    x <- gsScore[[id]]
    x <- x[is.finite(x)]
    n <- length(x)
    reason <- NA_character_

    ks_D <- NA_real_
    ks_p <- NA_real_
    shapiro_W <- NA_real_
    shapiro_p <- NA_real_

    if (n < 3) {
      reason <- "fewer than 3 finite observations"
    } else {
      x_sd <- stats::sd(x)
      if (is.na(x_sd) || x_sd == 0) {
        reason <- "zero variance"
      } else {
        if (method %in% c("ks", "both")) {
          z <- (x - mean(x)) / x_sd
          ks_res <- suppressWarnings(stats::ks.test(z, "pnorm"))
          ks_D <- unname(ks_res$statistic)
          ks_p <- ks_res$p.value
        }

        if (method %in% c("shapiro", "both")) {
          x_shapiro <- x
          if (n > shapiro_max_n) {
            x_shapiro <- sample(x_shapiro, shapiro_max_n)
          }
          shapiro_res <- stats::shapiro.test(x_shapiro)
          shapiro_W <- unname(shapiro_res$statistic)
          shapiro_p <- shapiro_res$p.value
        }
      }
    }

    data.frame(
      ID = id,
      normality_method = method,
      normality_n = n,
      normality_ks_D = ks_D,
      normality_ks_p = ks_p,
      normality_shapiro_W = shapiro_W,
      normality_shapiro_p = shapiro_p,
      normality_reason = reason,
      check.names = FALSE
    )
  })

  res <- do.call(rbind, diagnostics)
  rownames(res) <- NULL

  if (method %in% c("ks", "both")) {
    res$normality_ks_padj <- stats::p.adjust(res$normality_ks_p, method = p_adjust_method)
  }
  if (method %in% c("shapiro", "both")) {
    res$normality_shapiro_padj <- stats::p.adjust(res$normality_shapiro_p, method = p_adjust_method)
  }

  if (method == "shapiro") {
    res$normality_flag <- ifelse(is.na(res$normality_shapiro_padj), NA, res$normality_shapiro_padj < alpha)
  } else {
    res$normality_flag <- ifelse(is.na(res$normality_ks_padj), NA, res$normality_ks_padj < alpha)
  }

  if (method == "ks") {
    res <- res[, c(
      "ID",
      "normality_method",
      "normality_n",
      "normality_ks_D",
      "normality_ks_p",
      "normality_ks_padj",
      "normality_flag",
      "normality_reason"
    )]
  } else if (method == "shapiro") {
    res <- res[, c(
      "ID",
      "normality_method",
      "normality_n",
      "normality_shapiro_W",
      "normality_shapiro_p",
      "normality_shapiro_padj",
      "normality_flag",
      "normality_reason"
    )]
  } else {
    res <- res[, c(
      "ID",
      "normality_method",
      "normality_n",
      "normality_ks_D",
      "normality_ks_p",
      "normality_ks_padj",
      "normality_shapiro_W",
      "normality_shapiro_p",
      "normality_shapiro_padj",
      "normality_flag",
      "normality_reason"
    )]
  }

  return(res)
}
