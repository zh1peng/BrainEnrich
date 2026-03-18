#!/usr/bin/env Rscript

trim_arg <- function(x) {
  trimws(as.character(x))
}

split_arg <- function(x) {
  if (is.null(x) || !nzchar(trim_arg(x))) {
    return(NULL)
  }
  parts <- strsplit(trim_arg(x), ",", fixed = TRUE)[[1]]
  parts <- trim_arg(parts)
  parts[nzchar(parts)]
}

this_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
script_dir <- dirname(normalizePath(this_file))

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1L && nzchar(trim_arg(args[[1]]))) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
output_dir <- if (length(args) >= 2L && nzchar(trim_arg(args[[2]]))) {
  args[[2]]
} else {
  file.path(repo, "reports", "parity")
}
tiers_arg <- if (length(args) >= 3L) trim_arg(args[[3]]) else "tier1,tier2"
scenarios_arg <- if (length(args) >= 4L) trim_arg(args[[4]]) else ""
tolerance_arg <- if (length(args) >= 5L && nzchar(trim_arg(args[[5]]))) {
  trim_arg(args[[5]])
} else {
  "1e-6"
}

run_rscript <- function(script_name, script_args) {
  cmd <- file.path(R.home("bin"), "Rscript")
  status <- system2(cmd, c(file.path(script_dir, script_name), script_args))
  if (!identical(status, 0L)) {
    stop("Command failed for ", script_name, " with status ", status)
  }
}

export_head_snapshot <- function(repo_path, target_dir) {
  tar_path <- tempfile("brainenrich-head-", fileext = ".tar")
  on.exit(unlink(tar_path), add = TRUE)

  git_status <- system2(
    "git",
    c("-c", paste0("safe.directory=", repo_path), "-C", repo_path, "archive", "--format=tar", "HEAD", "-o", tar_path)
  )
  if (!identical(git_status, 0L)) {
    stop("Failed to export HEAD snapshot from ", repo_path)
  }
  utils::untar(tar_path, exdir = target_dir)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

baseline_snapshot <- tempfile("BrainEnrich-baseline-")
dir.create(baseline_snapshot, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(baseline_snapshot, recursive = TRUE, force = TRUE), add = TRUE)
export_head_snapshot(repo, baseline_snapshot)

baseline_rds <- file.path(output_dir, "parity_baseline.rds")
refactor_rds <- file.path(output_dir, "parity_refactor.rds")
compare_rds <- file.path(output_dir, "parity_compare.rds")
summary_md <- file.path(output_dir, "parity_summary.md")

message("Running baseline parity bundle from exported HEAD snapshot...")
run_rscript("run_parity_baseline.R", c(baseline_snapshot, baseline_rds, tiers_arg, scenarios_arg))

message("Running refactor parity bundle from current worktree...")
run_rscript("run_parity_refactor.R", c(repo, refactor_rds, tiers_arg, scenarios_arg))

message("Comparing baseline and refactor bundles...")
run_rscript("compare_parity.R", c(baseline_rds, refactor_rds, compare_rds, summary_md, tolerance_arg))

comparison <- readRDS(compare_rds)
failed <- comparison$summary$scenario[comparison$summary$status == "FAIL"]

if (length(failed) > 0L) {
  failed_csv <- paste(failed, collapse = ",")
  baseline_rerun_rds <- file.path(output_dir, "parity_baseline_rerun.rds")
  refactor_rerun_rds <- file.path(output_dir, "parity_refactor_rerun.rds")
  compare_rerun_rds <- file.path(output_dir, "parity_compare_rerun.rds")
  summary_rerun_md <- file.path(output_dir, "parity_summary_rerun.md")

  message("Initial parity failures detected. Rerunning failed scenarios once: ", failed_csv)
  run_rscript("run_parity_baseline.R", c(baseline_snapshot, baseline_rerun_rds, tiers_arg, failed_csv))
  run_rscript("run_parity_refactor.R", c(repo, refactor_rerun_rds, tiers_arg, failed_csv))
  run_rscript("compare_parity.R", c(baseline_rerun_rds, refactor_rerun_rds, compare_rerun_rds, summary_rerun_md, tolerance_arg))

  rerun <- readRDS(compare_rerun_rds)
  remaining <- rerun$summary$scenario[rerun$summary$status == "FAIL"]
  if (length(remaining) > 0L) {
    stop("Parity failed after rerun for scenarios: ", paste(remaining, collapse = ", "))
  }

  message("Parity passed after rerun. Rerun summary: ", summary_rerun_md)
} else {
  message("Parity passed without rerun.")
}

message("Baseline bundle: ", baseline_rds)
message("Refactor bundle: ", refactor_rds)
message("Comparison summary: ", summary_md)
