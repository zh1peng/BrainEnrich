#' Find core genes influencing aggregated score or LM coefficients
#'
#' Leave-one-out analysis to identify genes with the largest influence on
#' enrichment scores or standardized LM coefficients.
#'
#' @param geneList A matrix of genes by subs/models.
#' @param geneSetList A named list of gene sets.
#' @param pred_df Optional predictor data frame. If NULL with `cov_df`, group-level mode is used.
#' @param cov_df Optional covariate data frame.
#' @param aggre_method Aggregation method name or custom function.
#' @param n_cores Number of cores.
#' @param threshold_type Threshold method (`sd` or `percentile`).
#' @param threshold_value Threshold value.
#' @import pbapply parallel
#' @return A named list of character vectors with directional labels (`driver`/`buffer`).
#' The label is assigned by comparing the sign of the leave-one-out change with the sign of the
#' full statistic: `driver` means the gene reinforces the observed effect, while `buffer` means
#' the gene attenuates or opposes it.
#' @export
find_core_genes <- function(geneList, geneSetList, pred_df = NULL, cov_df = NULL, aggre_method, n_cores = 1, threshold_type = c("sd", "percentile"), threshold_value = 1) {
  type1_analysis <- is.null(pred_df) && is.null(cov_df)
  threshold_type <- match.arg(threshold_type)

  if (threshold_type == "sd" && (threshold_value < 0 || threshold_value > 3)) {
    stop("For 'sd', threshold_value should be between 0 and 3.")
  }
  if (threshold_type == "percentile" && (threshold_value < 1 || threshold_value > 99)) {
    stop("For 'percentile', threshold_value should be between 1 and 99.")
  }

  n_cores <- be_resolve_n_cores(n_cores)

  if (type1_analysis) {
    loo_stats <- be_parallel_lapply(seq_along(geneSetList), function(i) {
      gs <- geneSetList[[i]]
      full_score <- as.numeric(aggregate_geneSet(geneList, gs, method = aggre_method))
      changes <- sapply(seq_along(gs), function(gs_i) {
        modified_geneSet <- setdiff(gs, gs[gs_i])
        loo_score <- as.numeric(aggregate_geneSet(geneList, modified_geneSet, method = aggre_method))
        full_score - loo_score
      })
      names(changes) <- gs
      list(changes = changes, reference = full_score)
    },
    n_cores = n_cores,
    export = c("geneList", "geneSetList", "aggregate_geneSet", "aggre_method"),
    envir = environment()
    )
  } else if (is.data.frame(pred_df) && is.data.frame(cov_df)) {
    loo_stats <- be_parallel_lapply(seq_along(geneSetList), function(i) {
      gs <- geneSetList[[i]]
      gs_name <- names(geneSetList)[i]
      full_score <- aggregate_geneSet(geneList, gs, method = aggre_method)
      dependent_df.full <- data.frame(full_score)
      colnames(dependent_df.full) <- gs_name
      lm_res.full <- simple_lm(dependent_df = dependent_df.full, pred_df = pred_df, cov_df = cov_df, stat2return = "all")
      ref_coef <- lm_res.full$Standardized_Coefficient

      changes <- sapply(seq_along(gs), function(gs_i) {
        modified_geneSet <- setdiff(gs, gs[gs_i])
        loo_score <- aggregate_geneSet(geneList, modified_geneSet, method = aggre_method)
        dependent_df.loo <- data.frame(loo_score)
        colnames(dependent_df.loo) <- gs_name
        lm_res.loo <- simple_lm(dependent_df = dependent_df.loo, pred_df = pred_df, cov_df = cov_df, stat2return = "all")
        ref_coef - lm_res.loo$Standardized_Coefficient
      })
      names(changes) <- gs
      list(changes = changes, reference = ref_coef)
    },
    n_cores = n_cores,
    export = c("geneList", "geneSetList", "aggregate_geneSet", "aggre_method", "pred_df", "cov_df", "simple_lm"),
    envir = environment()
    )
  } else {
    stop("pred_df and cov_df must both be NULL or both be data.frames.")
  }

  core_genes <- lapply(loo_stats, function(x) {
    identify_core_genes(
      changes = x$changes,
      reference_stat = x$reference,
      threshold_type = threshold_type,
      threshold_value = threshold_value
    )
  })
  names(core_genes) <- names(geneSetList)
  core_genes
}


#' Identify core genes with directional labels
#'
#' @param changes Named vector of signed LOO changes.
#' @param reference_stat Numeric scalar of full model statistic (score or coef).
#' @param threshold_type Threshold method.
#' @param threshold_value Threshold value.
#' @importFrom stats quantile sd
#' @return Character vector such as `GENE(driver)`/`GENE(buffer)`.
#' The role is `driver` when the leave-one-out change has the same sign as the reference statistic
#' and `buffer` when it has the opposite sign. `neutral` is used when either sign is zero.
identify_core_genes <- function(changes, reference_stat, threshold_type = c("sd", "percentile"), threshold_value = 1) {
  threshold_type <- match.arg(threshold_type)
  impact <- abs(changes)

  if (threshold_type == "percentile") {
    if (threshold_value < 1 || threshold_value > 99) {
      stop("For 'percentile', threshold_value should be between 1 and 99.")
    }
    threshold <- quantile(impact, probs = threshold_value / 100)
  } else {
    if (threshold_value < 0 || threshold_value > 4) {
      stop("For 'sd', threshold_value should be between 0 and 4.")
    }
    threshold <- mean(impact) + threshold_value * sd(impact)
  }

  keep <- names(impact[impact > threshold & impact != 0])
  if (length(keep) == 0) {
    out <- "NOT_FOUND"
    attr(out, "impact") <- numeric(0)
    attr(out, "role") <- character(0)
    return(out)
  }

  ref_sign <- sign(reference_stat)
  roles <- vapply(keep, function(g) {
    delta_sign <- sign(changes[g])
    if (ref_sign == 0 || delta_sign == 0) {
      "neutral"
    } else if (delta_sign == ref_sign) {
      "driver"
    } else {
      "buffer"
    }
  }, FUN.VALUE = character(1))
  labels <- paste0(keep, "(", roles, ")")
  names(labels) <- keep
  attr(labels, "impact") <- unname(impact[keep])
  attr(labels, "role") <- unname(roles)
  labels
}
