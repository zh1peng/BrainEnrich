#' Distribute and Process Iterations Across HPC Jobs
#'
#' This function distributes a set of iterations across multiple jobs in a high-performance computing (HPC) environment.
#' It subsets the necessary variables, calls a specified function with the relevant subsets, and saves the results to a specified directory.
#'
#' @param job_id An integer specifying the job ID in the HPC job array (e.g., 1, 2, 3, ...). Determines which subset of iterations this job will process.
#' @param n_iter_per_job An integer specifying the number of iterations each job should process. Default is 1.
#' @param iter_total An integer specifying the total number of iterations to be processed across all jobs.
#' @param subset_vars A named list of variables (typically matrices) that need to be subset according to the job's assigned iterations.
#'        The names should correspond to the argument names in `FUN`.
#' @param subset_total_var A character string specifying the name of the argument in `FUN` that corresponds to the total number of iterations.
#'        If provided, the total number of iterations for the current job will be assigned to this argument.
#' @param output_dir A character string specifying the directory where the results should be saved. If the directory does not exist, it will be created.
#' @param prefix A character string specifying the base name for the saved RDS files. Default is "res_job_".
#' @param FUN A function that will be called to process the data for the specified iterations. The function should accept the arguments specified in `...`.
#' @param ... Additional arguments passed to `FUN`.
#' @return This function does not return a value but saves the results of `FUN` to an RDS file in the specified `output_dir`.
#' @examples
#' \dontrun{
#' # Call job_splitter, which will subset 'perm_id' and pass it to 'FUN'
#' job_splitter(
#'   job_id = 1,
#'   n_iter_per_job = 10,
#'   iter_total = 100,
#'   output_dir = "/path/to/output",
#'   FUN = brainscore,
#'   subset_vars = list(perm_id = perm_id),
#'   subset_total_var = "n_perm",
#'   brain_data = brain_data,
#'   gene_data = gene_data,
#'   annoData = annoData,
#'   cor_method = "pearson",
#'   aggre_method = "mean",
#'   null_model = "spin_brain",
#'   minGSSize = 20,
#'   maxGSSize = 200
#' )
#' }
#' @export
job_splitter <- function(job_id,
                         n_iter_per_job = 1,
                         iter_total,
                         prefix = "res_job_",
                         output_dir = NULL,
                         FUN,
                         subset_vars = list(),
                         subset_total_var = NULL,
                         ...) {
  # Check if output_dir is provided and create if it doesn't exist
  if (is.null(output_dir)) {
    stop("output_dir must be specified.")
  }

  if (!dir.exists(output_dir)) {
    message(sprintf("Output directory %s does not exist. Creating it...", output_dir))
    dir.create(output_dir, recursive = TRUE)
  }

  # Calculate the number of iterations this job should handle
  start_iter <- (job_id - 1) * n_iter_per_job + 1
  end_iter <- min(job_id * n_iter_per_job, iter_total)
  total_iter <- end_iter - start_iter + 1

  # Subset the variables in subset_vars for this job
  if (!length(subset_vars) == 0) {
    subset_vars_subset <- lapply(subset_vars, function(x) x[, start_iter:end_iter, drop = FALSE])

    # Replace the original variables with their subsets
    names(subset_vars_subset) <- names(subset_vars)
  } else {
    subset_vars_subset <- list()
  }
  # Prepare the argument list for FUN
  args <- c(subset_vars_subset, list(...))

  # Set the total_iter_var in args if specified
  if (!is.null(subset_total_var)) {
    args[[subset_total_var]] <- total_iter
  }
  message(sprintf("Job %d processing iterations %d to %d...", job_id, start_iter, end_iter))
  # Call the specified function with the arguments
  result <- do.call(FUN, args)
  # Save results to a file in the specified output directory
  saveRDS(result, file = file.path(output_dir, sprintf("%s%d.rds", prefix, job_id)))
}


#' Combine Results from Saved RDS Files
#'
#' This function combines results from multiple saved RDS files into a single data frame or list,
#' saves the combined results as an RDS file (optional), and checks for missing files.
#' @param input_dir A character string specifying the directory containing the RDS files to be combined.
#' @param output_dir A character string specifying the directory where the RDS files are stored.
#' @param n_rds An integer specifying the expected number of RDS files. If provided, the function will check if any files are missing.
#' @param save_name A character string specifying the name of the output RDS file for the combined results.
#'                  If not provided, it will default to the name of the output_dir. Must end with .rds.
#' @param file_pattern A character string specifying the pattern of the RDS files to be combined.
#'                     Default is "res_job_%d.rds".
#' @param delete_originals A logical indicating whether to delete the original RDS files after combining. Default is TRUE.
#' @param preserve_attributes A logical indicating whether to preserve and update attributes specific to brainscore output. Default is FALSE.
#' @param result_prefix A character string specifying the prefix for naming the combined results. Default is "null_".
#' @return The combined results, either as a saved RDS file or returned directly if `save_combined` is FALSE.
#' @export

job_cat <- function(input_dir,
                    output_dir = NULL,
                    n_rds = NULL,
                    save_name = NULL,
                    file_pattern = "res_job_%d.rds",
                    delete_originals = TRUE,
                    preserve_attributes = FALSE,
                    result_prefix = NULL) {
  # Set default save_name if not provided
  if (is.null(save_name)) {
    save_name <- basename(input_dir)
  }

  # Ensure save_name ends with .rds
  if (!grepl("\\.rds$", save_name)) {
    save_name <- paste0(save_name, ".rds")
  }

  # Get a list of all RDS files in the output directory matching the pattern
  rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

  if (length(rds_files) == 0) {
    stop("No RDS files found in the specified directory.")
  }

  # Check for missing files if n_rds is provided
  if (!is.null(n_rds)) {
    expected_files <- sprintf(file.path(input_dir, file_pattern), 1:n_rds)
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
    if (!ask_user_continue(msg)) {
      stop("Operation cancelled by user.")
    } else {
      delete_originals <- FALSE
    }
  }

  # Initialize a list to store the loaded results and attributes
  message("Loading RDS files...")
  all_results <- lapply(rds_files, readRDS)

  message("Combining results...")
  combined_results <- do.call(c, all_results)

  # Optionally preserve and update attributes for brainscore output
  if (preserve_attributes && length(all_results) > 0 && !is.null(attributes(all_results[[1]]))) {
    attributes_to_keep <- attributes(all_results[[1]])
    attributes_to_keep$n_perm <- length(combined_results) # Update n_perm to match the length of the combined list
    attributes(combined_results) <- attributes_to_keep
  }

  # Update the names of the combined results with the user-defined prefix
  if (!is.null(result_prefix) && length(combined_results) > 0 && is.list(combined_results) %% is.character(result_prefix)) {
    names(combined_results) <- paste0(result_prefix, 1:length(combined_results))
  }
  # Save the combined results as an RDS file
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      message(sprintf("Output directory %s does not exist. Creating it...", output_dir))
      dir.create(output_dir, recursive = TRUE)
    }
  } else {
    output_dir <- input_dir
  }
  combined_file_path <- file.path(output_dir, save_name)
  saveRDS(combined_results, combined_file_path)
  message(sprintf("Combined results saved to %s.", combined_file_path))

  # Optionally delete the original files
  if (delete_originals) {
    file.remove(rds_files)
    message("Original RDS files deleted.")
  }

  return(combined_results)
}
