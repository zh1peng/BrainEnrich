BE_DATA_PACKAGE <- "BrainEnrichData"
BE_DATA_VERSION <- "2.0.0"

#' Internal helper: user cache directory for BrainEnrich assets.
#'
#' @param ... Optional path components under cache root.
#' @return Character scalar path.
be_cache_dir <- function(...) {
  cache_root <- tryCatch(
    tools::R_user_dir("BrainEnrich", which = "cache"),
    error = function(e) file.path(path.expand("~"), ".cache", "BrainEnrich")
  )
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(cache_root, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

be_find_package_root <- function(package_name = "BrainEnrich", start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = FALSE)

  repeat {
    desc_file <- file.path(current, "DESCRIPTION")
    if (file.exists(desc_file)) {
      desc <- tryCatch(read.dcf(desc_file), error = function(e) NULL)
      if (!is.null(desc) &&
          "Package" %in% colnames(desc) &&
          identical(unname(desc[1, "Package"]), package_name)) {
        return(current)
      }
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }

  ""
}

be_description_version <- function(package_dir) {
  if (!nzchar(package_dir)) {
    return(NA_character_)
  }

  desc_file <- file.path(package_dir, "DESCRIPTION")
  if (!file.exists(desc_file)) {
    return(NA_character_)
  }

  desc <- tryCatch(read.dcf(desc_file), error = function(e) NULL)
  if (is.null(desc) || !"Version" %in% colnames(desc)) {
    return(NA_character_)
  }

  unname(desc[1, "Version"])
}

be_current_package_version <- function() {
  repo_root <- be_find_package_root("BrainEnrich")
  repo_version <- be_description_version(repo_root)
  if (nzchar(repo_version)) {
    return(repo_version)
  }

  tryCatch(
    as.character(utils::packageVersion("BrainEnrich")),
    error = function(e) NA_character_
  )
}

be_version_major <- function(version) {
  if (length(version) != 1 || is.na(version) || !nzchar(version)) {
    return(NA_integer_)
  }

  parts <- strsplit(as.character(version), "\\.", fixed = FALSE)[[1]]
  suppressWarnings(as.integer(parts[[1]]))
}

be_add_extdata_candidate <- function(candidates, root, label, version = NA_character_, enforce_major = FALSE) {
  if (!nzchar(root) || !dir.exists(root)) {
    return(candidates)
  }

  roots <- vapply(candidates, `[[`, character(1), "root", USE.NAMES = FALSE)
  if (root %in% roots) {
    return(candidates)
  }

  c(
    candidates,
    list(list(
      root = root,
      label = label,
      version = version,
      enforce_major = enforce_major
    ))
  )
}

be_extdata_candidates <- function() {
  candidates <- list()

  installed_data_root <- system.file("extdata", package = BE_DATA_PACKAGE)
  installed_data_version <- tryCatch(
    as.character(utils::packageVersion(BE_DATA_PACKAGE)),
    error = function(e) NA_character_
  )
  candidates <- be_add_extdata_candidate(
    candidates,
    installed_data_root,
    "installed BrainEnrichData package",
    version = installed_data_version,
    enforce_major = TRUE
  )

  repo_root <- be_find_package_root("BrainEnrich")
  sibling_pkg_root <- if (nzchar(repo_root)) file.path(dirname(repo_root), BE_DATA_PACKAGE) else ""
  sibling_data_root <- if (nzchar(sibling_pkg_root)) file.path(sibling_pkg_root, "inst", "extdata") else ""
  candidates <- be_add_extdata_candidate(
    candidates,
    sibling_data_root,
    "sibling BrainEnrichData checkout",
    version = be_description_version(sibling_pkg_root),
    enforce_major = TRUE
  )

  installed_pkg_root <- system.file("extdata", package = "BrainEnrich")
  installed_pkg_version <- tryCatch(
    as.character(utils::packageVersion("BrainEnrich")),
    error = function(e) NA_character_
  )
  candidates <- be_add_extdata_candidate(
    candidates,
    installed_pkg_root,
    "installed BrainEnrich package",
    version = installed_pkg_version
  )

  dev_pkg_root <- if (nzchar(repo_root)) file.path(repo_root, "inst", "extdata") else ""
  candidates <- be_add_extdata_candidate(
    candidates,
    dev_pkg_root,
    "BrainEnrich development checkout",
    version = be_description_version(repo_root)
  )

  candidates
}

be_validate_extdata_candidate <- function(candidate) {
  if (!isTRUE(candidate$enforce_major)) {
    return(invisible(candidate))
  }

  expected_major <- be_version_major(BE_DATA_VERSION)
  candidate_major <- be_version_major(candidate$version)
  if (!is.na(candidate_major) && !is.na(expected_major) && candidate_major != expected_major) {
    stop(
      "Found ", candidate$label, " at ", candidate$root,
      " but its version (", candidate$version, ") is incompatible with BrainEnrich ",
      be_current_package_version(), ". ",
      "Install BrainEnrichData ", expected_major, ".x to continue."
    )
  }

  invisible(candidate)
}

#' Internal helper: locate packaged extdata root.
#'
#' Lookup order:
#' 1) BrainEnrichData package extdata (if installed)
#' 2) Sibling BrainEnrichData checkout (`../BrainEnrichData/inst/extdata`)
#' 3) BrainEnrich package extdata (installed package)
#' 4) Development repo `inst/extdata` under the BrainEnrich checkout
#'
#' @return Character scalar directory path or empty string if not found.
be_extdata_root <- function() {
  candidates <- be_extdata_candidates()
  if (length(candidates) == 0) {
    return("")
  }

  be_validate_extdata_candidate(candidates[[1]])
  candidates[[1]]$root
}

#' Internal helper: resolve data file from package extdata or cache.
#'
#' @param relpath Relative path under extdata (e.g., "geneExp/desikan_r0.6.csv.bz2").
#' @return Named list with `path` and `from_pkg`.
be_resolve_data_file <- function(relpath) {
  for (candidate in be_extdata_candidates()) {
    pkg_path <- file.path(candidate$root, relpath)
    if (file.exists(pkg_path)) {
      be_validate_extdata_candidate(candidate)
      return(list(
        path = pkg_path,
        from_pkg = TRUE,
        source = candidate$label,
        version = candidate$version
      ))
    }
  }

  cache_path <- be_cache_dir(relpath)
  list(path = cache_path, from_pkg = FALSE, source = "user cache", version = NA_character_)
}

be_data_download_url <- function(relpath) {
  sprintf(
    "https://raw.githubusercontent.com/zh1peng/%s/v%s/inst/extdata/%s",
    BE_DATA_PACKAGE,
    BE_DATA_VERSION,
    relpath
  )
}
