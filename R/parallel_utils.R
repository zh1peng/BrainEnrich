#' Internal helper: resolve requested cores
#'
#' @param n_cores Requested number of cores. `0` means all detected cores minus one.
#'   When `_R_CHECK_LIMIT_CORES_` is set, the resolved value is capped at 2 to
#'   keep `R CMD check` compatible with CRAN-style core limits.
#' @return Integer scalar number of cores to use.
be_resolve_n_cores <- function(n_cores) {
  if (!is.numeric(n_cores) || length(n_cores) != 1L || is.na(n_cores)) {
    stop("n_cores must be a single numeric value.")
  }

  available <- parallel::detectCores()
  if (is.na(available) || available < 1L) {
    available <- 1L
  }

  max_allowed <- available
  check_limit <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  if (nzchar(check_limit) && check_limit != "false") {
    max_allowed <- min(max_allowed, 2L)
  }

  if (n_cores == 0) {
    return(max(min(available - 1L, max_allowed), 1L))
  }

  max(min(as.integer(n_cores), max_allowed), 1L)
}

#' Internal helper: apply a function with a safe parallel backend
#'
#' @param X Vector or list to iterate over.
#' @param FUN Function to apply.
#' @param n_cores Number of cores to use.
#' @param seed Optional seed for reproducibility.
#' @param export Character vector of objects to export for PSOCK workers.
#' @param envir Environment used when exporting objects.
#' @return A list of results.
be_parallel_lapply <- function(X,
                               FUN,
                               n_cores = 1,
                               seed = NULL,
                               export = character(),
                               envir = parent.frame()) {
  n_cores <- be_resolve_n_cores(n_cores)

  if (n_cores <= 1L) {
    return(pbapply::pblapply(X, FUN))
  }

  if (.Platform$OS.type != "windows") {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    return(parallel::mclapply(
      X,
      FUN,
      mc.cores = n_cores,
      mc.set.seed = !is.null(seed)
    ))
  }

  cl <- tryCatch(
    parallel::makeCluster(n_cores),
    error = function(e) {
      warning(
        "Parallel cluster startup failed; falling back to serial execution: ",
        conditionMessage(e),
        call. = FALSE
      )
      NULL
    }
  )

  if (is.null(cl)) {
    return(pbapply::pblapply(X, FUN))
  }

  on.exit(parallel::stopCluster(cl), add = TRUE)

  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, seed = seed)
  }
  if (length(export) > 0L) {
    parallel::clusterExport(cl, varlist = export, envir = envir)
  }

  pbapply::pblapply(X, FUN, cl = cl)
}
