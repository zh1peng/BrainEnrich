#!/usr/bin/env Rscript

this_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
script_dir <- dirname(normalizePath(this_file))
source(file.path(script_dir, "common.R"))

args <- parse_run_args("parity_refactor.rds")
load_repo_package(args$repo)
save_parity_run("refactor", args$repo, args$output, tiers = args$tiers, scenarios = args$scenarios)
