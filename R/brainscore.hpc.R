#' Calculate Brain Scores for Gene Sets (HPC version)
#'
#' This function calculates scores for gene sets based on brain data, designed for high-performance computing (HPC) environments where permutations are distributed across multiple jobs.
#'
#' @param job_id An integer specifying the job ID in the Slurm array (e.g., 1 to 1000).
#' @param n_perm_per_job An integer specifying the number of permutations to handle per job.
#' @param perm_total An integer specifying the total number of jobs (equal to n_perm in brainscore function).
#' @param perm_id A matrix of permutation indices. This will be subset according to the job's assigned permutations.
#' @param output_dir A character string specifying the directory where the results should be saved. If the directory does not exist, it will be created.
#' @param ... Additional arguments passed to the brainscore function, such as brain_data, gene_data, annoData, etc.
#' @return A list containing the gene set scores for the specific permutations handled by this job.
#' @export
brainscore.hpc <- function(job_id, n_perm_per_job = 1, perm_total, perm_id = NULL, output_dir = NULL, ...) {
  # Check if output_dir is provided and create if it doesn't exist
  if (is.null(output_dir)) {
    stop("output_dir must be specified.")
  }

  if (!dir.exists(output_dir)) {
    message(sprintf("Output directory %s does not exist. Creating it...", output_dir))
    dir.create(output_dir, recursive = TRUE)
  }

  # Calculate the number of permutations this job should handle
  start_perm <- (job_id - 1) * n_perm_per_job + 1
  end_perm <- min(job_id * n_perm_per_job, perm_total)

  # Subset the perm_id for this job
  if (!is.null(perm_id)) {
    perm_id_subset <- perm_id[, start_perm:end_perm, drop = FALSE]
  } else {
    perm_id_subset <- NULL
  }

  message(sprintf("Job %d processing permutations %d to %d...", job_id, start_perm, end_perm))

  # Call the brainscore function with the relevant permutations
  gs.score <- brainscore(
    n_perm = end_perm - start_perm + 1,
    perm_id = perm_id_subset,
    ...
  )
  # Save results to a file in the specified output directory
  saveRDS(gs.score, file = file.path(output_dir, sprintf("res_job_%d.rds", job_id)))
}

#' Combine Results from Saved RDS Files
#'
#' This function combines results from multiple saved RDS files into a single data frame or list,
#' saves the combined results as an RDS file (optional), and checks for missing files.
#'
#' @param output_dir A character string specifying the directory where the RDS files are stored.
#' @param n_rds An integer specifying the expected number of RDS files. If provided, the function will check if any files are missing.
#' @param save_name A character string specifying the name of the output RDS file for the combined results. Default is "combined_results.rds".
#' @param delete_originals A logical indicating whether to delete the original RDS files after combining. Default is TRUE.
#' @return The combined results, either as a saved RDS file or returned directly if `save_combined` is FALSE.
#' @export

combine.rds <- function(output_dir,
                        n_rds = NULL,
                        save_name = NULL,
                        delete_originals = TRUE) {
  # Get a list of all RDS files in the output directory
  rds_files <- list.files(output_dir, pattern = "\\.rds$", full.names = TRUE)

  if (length(rds_files) == 0) {
    stop("No RDS files found in the specified directory.")
  }

  # Check for missing files if n_rds is provided
  if (!is.null(n_rds)) {
    expected_files <- sprintf(file.path(output_dir, "res_job_%d.rds"), 1:n_rds)
    missing_files <- setdiff(expected_files, rds_files)

    if (length(missing_files) > 0) {
      message("The following files are missing:")
      print(missing_files)
      stop("Some RDS files are missing. Please ensure all files are present before combining.")
    } else {
      message("All expected RDS files are present.")
    }
  } else {
    message("Number of RDS files not specified.")
    msg <- paste0("There are ", length(rds_files), " RDS files that will be combined.")
    ask_user_continue(msg)
    delete_originals <- FALSE
  }

  # Initialize a list to store the loaded results and attributes
  message("Loading RDS files...")
  all_results <- lapply(rds_files, readRDS)

  message("Combining results...")
  combined_results <- do.call(c, all_results)



  # Preserve and update attributes
  if (length(all_results) > 0 && !is.null(attributes(all_results[[1]]))) {
    attributes_to_keep <- attributes(all_results[[1]])
    attributes_to_keep$n_perm <- length(combined_results) # Update n_perm to match the length of the combined list
    attributes(combined_results) <- attributes_to_keep
  }

  # Update the names of the combined results
  names(combined_results) <- paste0("null_", 1:length(combined_results))
  # Optionally save the combined results as an RDS file
  if (!is.null(save_name)) {
    # Ensure save_name ends with .rds
    if (!grepl("\\.rds$", save_name)) {
      save_name <- paste0(save_name, ".rds")
    }
    combined_file_path <- file.path(output_dir, save_name)
    saveRDS(combined_results, combined_file_path)
    message(sprintf("Combined results saved to %s.", combined_file_path))

    # Optionally delete the original files
    if (delete_originals) {
      file.remove(rds_files)
      message("Original RDS files deleted.")
    }
  }

  return(combined_results)
}
