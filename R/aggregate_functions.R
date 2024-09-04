#' Aggregate Gene Set Scores
#'
#' Function to aggregate geneList based on geneSet, evaluating one geneSet at a time.
#' The function supports multiple aggregation methods as specified by the user.
#'
#' @param geneList A matrix of genes by models, with each column representing a true or null model.
#' @param geneSet A vector containing names of genes in the gene set of interest.
#' @param method A character string specifying the method to use for aggregation.
#'               Options include 'mean', 'median', 'meanabs', 'meansqr', 'maxmean', 'sign_test', 'rank_sum', 'ks_orig', 'ks_weighted', 'ks_sum'. Default is 'mean'.
#' @importFrom stats median p.adjust
#' @importFrom methods new
#' @return Returns a numeric score based on the specified aggregation method.
#' @export
#'
#'
aggregate_geneSet <- function(geneList, # named correlation/coefficient matrix
                              geneSet, # one geneSet of interest
                              method = c("mean", "median", "meanabs", "meansqr", "maxmean", "ks_orig", "ks_weighted", "ks_pos_neg_sum", "sign_test", "rank_sum", "custom")) {
  if (is.function(method)) { # if method is a custom function
    aggre_func <- method
    method <- "custom"
  } else {
    method <- match.arg(method)
  }

  if (!is.matrix(geneList)) {
    stop("geneList should be a matrix of genes*models (true model or null models)")
  }

  aggre_func <- switch(method,
    mean = {
      function(genelist, geneSet) {
        geneSet <- intersect(geneSet, names(genelist))
        hits <- names(genelist) %in% geneSet
        res <- mean(genelist[hits])
        return(res)
      }
    },
    median = {
      function(genelist, geneSet) {
        geneSet <- intersect(geneSet, names(genelist))
        hits <- names(genelist) %in% geneSet
        res <- median(genelist[hits])
        return(res)
      }
    },
    meanabs = {
      function(genelist, geneSet) {
        geneSet <- intersect(geneSet, names(genelist))
        hits <- names(genelist) %in% geneSet
        res <- mean(abs(genelist[hits]))
      }
    },
    meansqr = {
      function(genelist, geneSet) {
        geneSet <- intersect(geneSet, names(genelist))
        hits <- names(genelist) %in% geneSet
        res <- mean(genelist[hits]^2)
        return(res)
      }
    },
    maxmean = {
      function(genelist, geneSet) {
        geneSet <- intersect(geneSet, names(genelist))
        hits <- names(genelist) %in% geneSet
        X <- genelist[hits]
        # Mean of positive numbers
        pos.mean <- mean(X[X > 0])
        # Mean of negative numbers
        neg.mean <- mean(X[X < 0])
        if (is.na(pos.mean)) {
          pos.mean <- 0
        } # when all value are neg
        if (is.na(neg.mean)) {
          neg.mean <- 0
        } # when all value are pos
        res <- ifelse(pos.mean > abs(neg.mean), pos.mean, neg.mean)
        return(res)
      }
    },
    ks_orig = {
      function(genelist, geneSet) {
        # code obtained from DOSE: https://rdrr.io/bioc/DOSE/src/R/gsea.R
        genelist <- sort(genelist, decreasing = TRUE)
        geneSet <- intersect(geneSet, names(genelist))
        N <- length(genelist)
        Nh <- length(geneSet)
        Phit <- Pmiss <- numeric(N)
        hits <- names(genelist) %in% geneSet
        Phit[hits] <- abs(genelist[hits])^0 # raw rank
        NR <- sum(Phit)
        Phit <- cumsum(Phit / NR)
        Pmiss[!hits] <- 1 / (N - Nh)
        Pmiss <- cumsum(Pmiss)
        runningES <- Phit - Pmiss
        max.ES <- max(runningES)
        min.ES <- min(runningES)
        if (abs(max.ES) > abs(min.ES)) {
          res <- max.ES
        } else {
          res <- min.ES
        }
        return(res)
      }
    },
    ks_weighted = {
      function(genelist, geneSet) {
        # code obtained from DOSE: https://rdrr.io/bioc/DOSE/src/R/gsea.R
        genelist <- sort(genelist, decreasing = TRUE)
        geneSet <- intersect(geneSet, names(genelist))
        N <- length(genelist)
        Nh <- length(geneSet)
        Phit <- Pmiss <- numeric(N)
        hits <- names(genelist) %in% geneSet
        Phit[hits] <- abs(genelist[hits])^1 # weighted rank
        NR <- sum(Phit)
        Phit <- cumsum(Phit / NR)
        Pmiss[!hits] <- 1 / (N - Nh)
        Pmiss <- cumsum(Pmiss)
        runningES <- Phit - Pmiss
        max.ES <- max(runningES)
        min.ES <- min(runningES)
        if (abs(max.ES) > abs(min.ES)) {
          res <- max.ES
        } else {
          res <- min.ES
        }
        return(res)
      }
    },
    ks_pos_neg_sum = {
      function(genelist, geneSet) {
        # function borrowed from clusterprofiler
        genelist <- sort(genelist, decreasing = TRUE)
        geneSet <- intersect(geneSet, names(genelist))
        N <- length(genelist)
        Nh <- length(geneSet)
        Phit <- Pmiss <- numeric(N)
        hits <- names(genelist) %in% geneSet
        Phit[hits] <- abs(genelist[hits])^1 # weighted rank
        NR <- sum(Phit)
        Phit <- cumsum(Phit / NR)
        Pmiss[!hits] <- 1 / (N - Nh)
        Pmiss <- cumsum(Pmiss)
        runningES <- Phit - Pmiss
        max.ES <- max(runningES)
        min.ES <- min(runningES)
        res <- max.ES + min.ES
        return(res)
      }
    },
    sign_test = {
      function(genelist, geneSet) {
        geneSet <- intersect(geneSet, names(genelist))
        hits <- names(genelist) %in% geneSet
        X <- genelist[hits]
        n.pos <- sum(X > 0)
        n.neg <- sum(X < 0)
        n.smaller <- ifelse(n.pos < n.neg, n.pos, n.neg)
        n.total <- n.pos + n.neg
        res <- ifelse(n.total <= 25, n.smaller, ((n.smaller + 0.5) - (n.total / 2)) / (sqrt(n.total) / 2))
        return(res)
      }
    },
    # local_fdr = {
    #   function(genelist, geneSet) {
    #     geneSet <- intersect(geneSet, names(genelist))
    #     hits <- names(genelist) %in% geneSet ## logical
    #     df <- data.frame(vals = as.numeric(genelist), hits = as.numeric(hits))
    #     fit.pos <- glm(hits ~ vals, family = binomial(link = "logit"), data = df[df$vals >= 0, ])
    #     fit.neg <- glm(hits ~ vals, family = binomial(link = "logit"), data = df[df$vals < 0, ])
    #     S.pos <- coef(summary(fit.pos))["vals", "z value"]
    #     S.neg <- coef(summary(fit.neg))["vals", "z value"]
    #     res <- ifelse(abs(S.pos) > abs(S.neg), S.pos, S.neg)
    #     return(res)
    #   }
    # },
    rank_sum = {
      function(genelist, geneSet) {
        abslist.sorted <- sort(abs(genelist), decreasing = F) # sort by abs value
        ranklist <- c(1:length(genelist)) # create rank
        names(ranklist) <- names(abslist.sorted) # set the gene names
        geneSet <- intersect(geneSet, names(genelist))
        hits <- names(genelist) %in% geneSet
        pos <- genelist > 0
        neg <- genelist < 0
        pos.name <- names(genelist)[hits & pos]
        neg.name <- names(genelist)[hits & neg]
        pos.ranksum <- sum(ranklist[pos.name])
        neg.ranksum <- sum(ranklist[neg.name])
        res <- ifelse(pos.ranksum < neg.ranksum, pos.ranksum, neg.ranksum)
        return(res)
      }
    },
    custom = {
      aggre_func
    },
    stop("Invalid aggregation method name")
  )

  if (!is.function(aggre_func)) {
    stop("Invalid method name")
  }

  gs_score <- sapply(1:ncol(geneList), function(x) {
    aggre_func(genelist = geneList[, x], geneSet = geneSet)
  })
  return(gs_score)
}



#' Aggregate Gene Set List in Parallel
#'
#' This function aggregates gene sets in parallel using the `parLapply` function
#' from the `parallel` package, ensuring cross-platform compatibility.
#'
#' @param geneList A list of genes.
#' @param geneSetList A list of gene sets to be aggregated.
#' @param n_cores Number of cores to use for parallel processing. Default is 1.
#' If set to 0, it uses all available cores minus one.
#' @param method aggregation method used.
#' @param n_cores Number of cores to use for parallel processing. Default is 1.
#' @return A list of aggregated gene set scores.
#' @import pbapply
#' @import parallel
#' @export
aggregate_geneSetList <- function(geneList, geneSetList, method, n_cores = 1) {
  # Determine the number of cores to use
  if (n_cores == 0) {
    n_cores <- max(detectCores() - 1, 1) # Use all cores minus one, but ensure at least 1 core is used
  } else {
    n_cores <- min(n_cores, detectCores()) # Ensure n_cores does not exceed the number of available cores
  }

  # Initialize a cluster of workers
  if (n_cores == 1) {
    cl <- NULL
  } else {
    cl <- makeCluster(n_cores)
    # Export necessary variables to the cluster
    clusterExport(cl, c("geneList", "aggregate_geneSet", "geneSetList", "method"),
      envir = environment()
    )
  }
  # Parallelize the processing using pblapply for progress bar
  allgs.scores <- pblapply(seq_along(geneSetList), function(i) {
    gs <- geneSetList[[i]]
    aggregate_geneSet(geneList = geneList, geneSet = gs, method = method)
  }, cl = cl)

  # Stop the cluster after processing
  if (!is.null(cl)) {
    stopCluster(cl)
  }
  names(allgs.scores) <- names(geneSetList)
  return(allgs.scores)
}




#' Aggregate Gene Set List with Matching Coexpression in Parallel
#'
#' This function aggregates gene set scores in parallel using `pblapply` from the `pbapply` package.
#'
#' @param geneList.true A m x 1 matrix of true gene sets. Ensure to include `drop=FALSE` when subsetting.
#' @param geneSetList A list of gene sets.
#' @param sampled_geneSetList A list of sampled gene sets.
#' @param method The method to be used for aggregation.
#' @param n_cores Number of cores to use for parallel processing. Default is 1.
#' If set to 0, it uses all available cores minus one.
#' @return A list of aggregated gene set scores.
#' @import parallel
#' @import pbapply
#' @export
aggregate_geneSetList_matching_coexp <- function(geneList.true,
                                                 geneSetList,
                                                 sampled_geneSetList,
                                                 method,
                                                 n_cores = 1) {
  # Ensure geneList.true is a matrix with one column
  if (!is.matrix(geneList.true) || ncol(geneList.true) != 1) {
    stop("geneList.true should be a m x 1 matrix. Please include drop=FALSE when subsetting.")
  }

  # Ensure the order in geneSetList and sampled_geneSetList are the same
  if (!identical(names(geneSetList), names(sampled_geneSetList))) {
    stop("geneSetList and sampled_geneSetList are not matched.")
  }


  # Determine the number of cores to use
  if (n_cores == 0) {
    n_cores <- max(detectCores() - 1, 1) # Use all cores minus one, but ensure at least 1 core is used
  } else {
    n_cores <- min(n_cores, detectCores()) # Ensure n_cores does not exceed the number of available cores
  }


  # Initialize a cluster of workers
  cl <- makeCluster(n_cores)

  # Export necessary variables and functions to the cluster
  clusterExport(cl,
    varlist = c("geneList.true", "swap_geneList", "aggregate_geneSet", "method", "geneSetList", "sampled_geneSetList"),
    envir = environment()
  )

  # Parallelize the processing using pblapply for progress bar
  allgs.scores <- pblapply(seq_along(geneSetList), function(i) {
    gs <- geneSetList[[i]]
    sampled_gs <- sampled_geneSetList[[i]]
    geneList.null <- swap_geneList(
      geneList.true = geneList.true,
      orig_gs = gs,
      sampled_gs = sampled_gs
    )
    gs.score <- aggregate_geneSet(
      geneList = geneList.null,
      geneSet = gs,
      method = method
    )
    return(gs.score)
  }, cl = cl)

  # Stop the cluster after processing
  stopCluster(cl)
  names(allgs.scores) <- names(geneSetList)
  return(allgs.scores)
}
