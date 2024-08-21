#' Calculate Brain Scores for Gene Sets (HPC version)
#'
#' This function calculates scores for gene sets based on brain data, designed for high-performance computing (HPC) environments where permutations are distributed across multiple jobs.
#'
#' @param job_id An integer specifying the job ID in the Slurm array (e.g., 1 to 1000).
#' @param perm_total An integer specifying the total number of jobs (equal to n_perm in brainscore function).
#' @param perm_id A matrix of permutation indices. This will be subset according to the job's assigned permutations.
#' @param output_dir A character string specifying the directory where the results should be saved. If the directory does not exist, it will be created.
#' @param ... Additional arguments passed to the brainscore function, such as brain_data, gene_data, annoData, etc.
#' @return A list containing the gene set scores for the specific permutations handled by this job.
#' @export
brainscore.hpc <- function(job_id, perm_total, perm_id = NULL, output_dir = NULL, ...) {
  # Check if output_dir is provided and create if it doesn't exist
  if (is.null(output_dir)) {
    stop("output_dir must be specified.")
  }

  if (!dir.exists(output_dir)) {
    message(sprintf("Output directory %s does not exist. Creating it...", output_dir))
    dir.create(output_dir, recursive = TRUE)
  }

  # Calculate the number of permutations this job should handle
  n_perm_per_job <- 1
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
  saveRDS(gs.score, file = file.path(output_dir, sprintf("brainscore_results_job_%d.rds", job_id)))
}
