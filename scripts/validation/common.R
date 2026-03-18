if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Package `pkgload` is required to run validation scripts.")
}

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

parity_settings <- function() {
  list(
    seed = 20260317L,
    n_cores = 1L,
    atlas = "desikan",
    rdonor = "r0.6",
    hem = "L",
    sim_rows = 80L,
    gene_set_source = "SynGO",
    minGSSize = 20L,
    maxGSSize = 200L,
    lm_minGSSize = 50L,
    pvalueCutoff = 1,
    threshold_type = "sd",
    threshold_value = 1,
    matchcoexp_tol = 0.8,
    matchcoexp_max_iter = 1000000L,
    tiers = list(
      tier1 = list(
        n_perm = 10L,
        sim_n = 3L,
        subsample_size = 50L
      ),
      tier2 = list(
        n_perm = 100L,
        sim_n = 10L,
        subsample_size = 50L
      )
    )
  )
}

scenario_catalog <- function() {
  list(
    tier1_brainenrich_spin_brain = list(workflow = "brainenrich", tier = "tier1"),
    tier1_brainenrich_resample_gene = list(workflow = "brainenrich", tier = "tier1"),
    tier1_brainscore_none = list(workflow = "brainscore", tier = "tier1"),
    tier1_brainscore_spin_brain = list(workflow = "brainscore", tier = "tier1"),
    tier1_brainscore_resample_gene = list(workflow = "brainscore", tier = "tier1"),
    tier1_brainscore_coexp_matched = list(workflow = "brainscore", tier = "tier1"),
    tier1_brainscore_lm_test_object_spin_brain = list(workflow = "brainscore.lm_test", tier = "tier1"),
    tier1_brainscore_lm_test_table_spin_brain = list(workflow = "brainscore.lm_test", tier = "tier1"),
    tier1_brainscore_simulate_randomize_pred_type1 = list(workflow = "brainscore.simulate", tier = "tier1"),
    tier1_brainscore_simulate_randomize_pred_power = list(workflow = "brainscore.simulate", tier = "tier1"),
    tier1_brainscore_simulate_spin_brain_type1 = list(workflow = "brainscore.simulate", tier = "tier1"),
    tier1_brainscore_simulate_spin_brain_power = list(workflow = "brainscore.simulate", tier = "tier1"),
    tier1_brainscore_simulate_resample_gene_type1 = list(workflow = "brainscore.simulate", tier = "tier1"),
    tier1_brainscore_simulate_resample_gene_power = list(workflow = "brainscore.simulate", tier = "tier1"),
    tier2_brainenrich_spin_brain = list(workflow = "brainenrich", tier = "tier2"),
    tier2_brainscore_coexp_matched = list(workflow = "brainscore", tier = "tier2"),
    tier2_brainscore_lm_test_object_spin_brain = list(workflow = "brainscore.lm_test", tier = "tier2"),
    tier2_brainscore_simulate_spin_brain_type1 = list(workflow = "brainscore.simulate", tier = "tier2")
  )
}

parse_run_args <- function(default_output_name) {
  args <- commandArgs(trailingOnly = TRUE)
  repo <- if (length(args) >= 1L && nzchar(trim_arg(args[[1]]))) {
    args[[1]]
  } else {
    getwd()
  }
  output <- if (length(args) >= 2L && nzchar(trim_arg(args[[2]]))) {
    args[[2]]
  } else {
    file.path(repo, "reports", default_output_name)
  }
  tiers <- if (length(args) >= 3L) split_arg(args[[3]]) else c("tier1", "tier2")
  scenarios <- if (length(args) >= 4L) split_arg(args[[4]]) else NULL

  if (is.null(tiers) || length(tiers) == 0L) {
    tiers <- c("tier1", "tier2")
  }

  list(
    repo = normalizePath(repo, mustWork = TRUE),
    output = output,
    tiers = tiers,
    scenarios = scenarios
  )
}

validate_tiers <- function(tiers) {
  valid <- names(parity_settings()$tiers)
  if (!all(tiers %in% valid)) {
    stop("Unknown tiers: ", paste(setdiff(tiers, valid), collapse = ", "))
  }
  tiers
}

resolve_scenarios <- function(tiers = c("tier1", "tier2"), scenarios = NULL) {
  tiers <- validate_tiers(tiers)
  catalog <- scenario_catalog()
  selected <- names(catalog)[vapply(catalog, function(x) x$tier %in% tiers, FUN.VALUE = logical(1))]

  if (is.null(scenarios)) {
    return(selected)
  }

  unknown <- setdiff(scenarios, names(catalog))
  if (length(unknown) > 0L) {
    stop("Unknown scenarios: ", paste(unknown, collapse = ", "))
  }
  scenarios
}

load_repo_package <- function(repo) {
  pkgload::load_all(repo, quiet = TRUE, export_all = FALSE)
}

build_shared_inputs <- function(settings) {
  data("brain_data", package = "BrainEnrich")
  data("sim_hcp", package = "BrainEnrich")
  data("perm_id_dk_lh_5000", package = "BrainEnrich")

  gene_data <- get_geneExp(
    atlas = settings$atlas,
    rdonor = settings$rdonor,
    hem = settings$hem
  )
  annoData <- get_annoData(type = settings$gene_set_source)

  sim_hcp_small <- sim_hcp[seq_len(min(settings$sim_rows, nrow(sim_hcp))), , drop = FALSE]
  brain_data_individual <- t(dplyr::select(sim_hcp_small, dplyr::starts_with("L_")))
  colnames(brain_data_individual) <- paste0("sub-", seq_len(ncol(brain_data_individual)))

  list(
    group_brain_data = brain_data,
    individual_brain_data = brain_data_individual,
    gene_data = gene_data,
    annoData = annoData,
    perm_id = perm_id_dk_lh_5000,
    pred_df = dplyr::select(sim_hcp_small, BMI),
    cov_df = dplyr::select(sim_hcp_small, Age, Sex)
  )
}

scenario_seed <- function(scenario_name) {
  parity_settings()$seed + match(scenario_name, names(scenario_catalog())) - 1L
}

run_with_seed <- function(seed, expr) {
  set.seed(seed)
  force(expr)
}

run_brainenrich_scenario <- function(shared, settings, tier_settings, scenario_name, null_model, perm_id = NULL) {
  list(
    workflow = "brainenrich",
    tier = sub("_.*$", "", scenario_name),
    seed = scenario_seed(scenario_name),
    parameters = list(null_model = null_model, n_perm = tier_settings$n_perm),
    result = run_with_seed(
      scenario_seed(scenario_name),
      brainenrich(
        brain_data = shared$group_brain_data,
        gene_data = shared$gene_data,
        annoData = shared$annoData,
        cor_method = "pearson",
        aggre_method = "mean",
        null_model = null_model,
        n_perm = tier_settings$n_perm,
        perm_id = perm_id,
        n_cores = settings$n_cores,
        minGSSize = settings$minGSSize,
        maxGSSize = settings$maxGSSize,
        pvalueCutoff = settings$pvalueCutoff,
        threshold_type = settings$threshold_type,
        threshold_value = settings$threshold_value
      )
    )
  )
}

run_brainscore_scenario <- function(shared, settings, tier_settings, scenario_name, null_model, perm_id = NULL) {
  list(
    workflow = "brainscore",
    tier = sub("_.*$", "", scenario_name),
    seed = scenario_seed(scenario_name),
    parameters = list(
      null_model = null_model,
      n_perm = if (identical(null_model, "none")) NULL else tier_settings$n_perm
    ),
    result = run_with_seed(
      scenario_seed(scenario_name),
      brainscore(
        brain_data = shared$individual_brain_data,
        gene_data = shared$gene_data,
        annoData = shared$annoData,
        cor_method = "pearson",
        aggre_method = "mean",
        null_model = null_model,
        minGSSize = settings$minGSSize,
        maxGSSize = settings$maxGSSize,
        n_cores = settings$n_cores,
        n_perm = if (identical(null_model, "none")) NULL else tier_settings$n_perm,
        perm_id = perm_id,
        matchcoexp_tol = settings$matchcoexp_tol,
        matchcoexp_max_iter = settings$matchcoexp_max_iter,
        verbose = FALSE
      )
    )
  )
}

run_lm_test_scenario <- function(shared, settings, tier_settings, scenario_name, gsea_obj, perm_id) {
  list(
    workflow = "brainscore.lm_test",
    tier = sub("_.*$", "", scenario_name),
    seed = scenario_seed(scenario_name),
    parameters = list(null_model = "spin_brain", n_perm = tier_settings$n_perm, gsea_obj = gsea_obj),
    result = run_with_seed(
      scenario_seed(scenario_name),
      brainscore.lm_test(
        pred_df = shared$pred_df,
        cov_df = shared$cov_df,
        brain_data = shared$individual_brain_data,
        gene_data = shared$gene_data,
        annoData = shared$annoData,
        cor_method = "pearson",
        aggre_method = "mean",
        n_cores = settings$n_cores,
        minGSSize = settings$lm_minGSSize,
        maxGSSize = settings$maxGSSize,
        null_model = "spin_brain",
        n_perm = tier_settings$n_perm,
        perm_id = perm_id,
        pvalueCutoff = settings$pvalueCutoff,
        threshold_type = settings$threshold_type,
        threshold_value = settings$threshold_value,
        gsea_obj = gsea_obj
      )
    )
  )
}

run_simulate_scenario <- function(shared, settings, tier_settings, scenario_name, sim_type, sim_setting, perm_id = NULL) {
  list(
    workflow = "brainscore.simulate",
    tier = sub("_.*$", "", scenario_name),
    seed = scenario_seed(scenario_name),
    parameters = list(
      sim_type = sim_type,
      sim_setting = sim_setting,
      sim_n = tier_settings$sim_n,
      n_perm = tier_settings$n_perm,
      subsample_size = tier_settings$subsample_size
    ),
    result = run_with_seed(
      scenario_seed(scenario_name),
      brainscore.simulate(
        pred_df = shared$pred_df,
        cov_df = shared$cov_df,
        brain_data = shared$individual_brain_data,
        gene_data = shared$gene_data,
        annoData = shared$annoData,
        sim_n = tier_settings$sim_n,
        subsample_size = tier_settings$subsample_size,
        sim_type = sim_type,
        sim_setting = sim_setting,
        cor_method = "pearson",
        aggre_method = "mean",
        minGSSize = settings$minGSSize,
        maxGSSize = settings$maxGSSize,
        n_cores = settings$n_cores,
        perm_id = perm_id,
        n_perm = tier_settings$n_perm
      )
    )
  )
}

run_scenario <- function(scenario_name, shared, settings) {
  spec <- scenario_catalog()[[scenario_name]]
  tier_settings <- settings$tiers[[spec$tier]]
  perm_id <- shared$perm_id[, seq_len(tier_settings$n_perm), drop = FALSE]

  switch(scenario_name,
    tier1_brainenrich_spin_brain = run_brainenrich_scenario(shared, settings, tier_settings, scenario_name, "spin_brain", perm_id = perm_id),
    tier1_brainenrich_resample_gene = run_brainenrich_scenario(shared, settings, tier_settings, scenario_name, "resample_gene"),
    tier1_brainscore_none = run_brainscore_scenario(shared, settings, tier_settings, scenario_name, "none"),
    tier1_brainscore_spin_brain = run_brainscore_scenario(shared, settings, tier_settings, scenario_name, "spin_brain", perm_id = perm_id),
    tier1_brainscore_resample_gene = run_brainscore_scenario(shared, settings, tier_settings, scenario_name, "resample_gene"),
    tier1_brainscore_coexp_matched = run_brainscore_scenario(shared, settings, tier_settings, scenario_name, "coexp_matched"),
    tier1_brainscore_lm_test_object_spin_brain = run_lm_test_scenario(shared, settings, tier_settings, scenario_name, gsea_obj = TRUE, perm_id = perm_id),
    tier1_brainscore_lm_test_table_spin_brain = run_lm_test_scenario(shared, settings, tier_settings, scenario_name, gsea_obj = FALSE, perm_id = perm_id),
    tier1_brainscore_simulate_randomize_pred_type1 = run_simulate_scenario(shared, settings, tier_settings, scenario_name, sim_type = "randomize_pred", sim_setting = "type1"),
    tier1_brainscore_simulate_randomize_pred_power = run_simulate_scenario(shared, settings, tier_settings, scenario_name, sim_type = "randomize_pred", sim_setting = "power"),
    tier1_brainscore_simulate_spin_brain_type1 = run_simulate_scenario(shared, settings, tier_settings, scenario_name, sim_type = "spin_brain", sim_setting = "type1", perm_id = perm_id),
    tier1_brainscore_simulate_spin_brain_power = run_simulate_scenario(shared, settings, tier_settings, scenario_name, sim_type = "spin_brain", sim_setting = "power", perm_id = perm_id),
    tier1_brainscore_simulate_resample_gene_type1 = run_simulate_scenario(shared, settings, tier_settings, scenario_name, sim_type = "resample_gene", sim_setting = "type1"),
    tier1_brainscore_simulate_resample_gene_power = run_simulate_scenario(shared, settings, tier_settings, scenario_name, sim_type = "resample_gene", sim_setting = "power"),
    tier2_brainenrich_spin_brain = run_brainenrich_scenario(shared, settings, tier_settings, scenario_name, "spin_brain", perm_id = perm_id),
    tier2_brainscore_coexp_matched = run_brainscore_scenario(shared, settings, tier_settings, scenario_name, "coexp_matched"),
    tier2_brainscore_lm_test_object_spin_brain = run_lm_test_scenario(shared, settings, tier_settings, scenario_name, gsea_obj = TRUE, perm_id = perm_id),
    tier2_brainscore_simulate_spin_brain_type1 = run_simulate_scenario(shared, settings, tier_settings, scenario_name, sim_type = "spin_brain", sim_setting = "type1", perm_id = perm_id),
    stop("Unhandled scenario: ", scenario_name)
  )
}

save_parity_run <- function(label, repo, output, tiers = c("tier1", "tier2"), scenarios = NULL) {
  settings <- parity_settings()
  shared <- build_shared_inputs(settings)
  scenario_names <- resolve_scenarios(tiers = tiers, scenarios = scenarios)

  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  results <- lapply(scenario_names, function(scenario_name) {
    message("[", label, "] Running ", scenario_name)
    run_scenario(scenario_name, shared, settings)
  })
  names(results) <- scenario_names

  payload <- list(
    label = label,
    repo = normalizePath(repo, mustWork = TRUE),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    settings = settings,
    scenarios = results
  )
  saveRDS(payload, output)
  message("Saved parity run to ", output)
}
