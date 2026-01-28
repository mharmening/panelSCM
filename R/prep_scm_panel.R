#' Prepare Synthetic Control Models for Panel Data
#'
#' Fits unit-specific synthetic control models for multiple treated units
#' and determines the optimal donor pool size by minimizing the pre-treatment
#' normalized mean squared error.
#'
#' @param df A data frame containing panel data.
#' @param outcome_var Name of the outcome variable.
#' @param predictor_vars Character vector of predictor variable names.
#' @param id_var Name of the unit identifier.
#' @param name_var Name of the unit name variable.
#' @param time_var Name of the time variable.
#' @param treatment_indicator Name of the treatment indicator (binary).
#' @param donor_range Integer vector of donor pool sizes to evaluate.
#' @param parallel Logical. Whether to use parallel processing.
#'
#' @return A list containing fitted synthetic control models, data objects,
#' diagnostic fit measures, and treatment timing information.
#'
#' @export





prep_scm_panel <- function(df,
                           outcome_var,
                           predictor_vars,
                           id_var,
                           name_var,
                           time_var,
                           treatment_indicator,
                           donor_range = 3:5,
                           parallel = TRUE) {


  # --- 0. COERCE UNIT ID TO CHARACTER ---
  df[[id_var]] <- as.character(df[[id_var]])

  # --- 1. INITIAL IDENTIFICATION ---
  treated_unit_ids <- unique(df[[id_var]][df[[treatment_indicator]] == 1])
  n_treated <- length(treated_unit_ids)
  all_unit_ids <- unique(df[[id_var]])
  potential_donor_ids <- all_unit_ids[!(all_unit_ids %in% treated_unit_ids)]

  # --- 2. PHASE 1: PRE CALCULATE CORRELATIONS ---
  # Ranking donors by pre-treatment correlation with the target unit
  message("\n[Phase 1/3] Calculating donor similarities (Correlation ranking)...")

  donor_rankings <- list()
  pb1 <- txtProgressBar(min = 0, max = n_treated, style = 3)

  for (i in 1:n_treated) {
    target_id <- treated_unit_ids[i]
    t_treat <- min(df[[time_var]][df[[id_var]] == target_id & df[[treatment_indicator]] == 1])
    pre_period <- sort(unique(df[[time_var]][df[[time_var]] < t_treat]))
    y_target_pre <- df[[outcome_var]][df[[id_var]] == target_id & df[[time_var]] %in% pre_period]

    cor_vals <- sapply(potential_donor_ids, function(cid) {
      y_control_pre <- df[[outcome_var]][df[[id_var]] == cid & df[[time_var]] %in% pre_period]
      return(cor(y_target_pre, y_control_pre, use = "complete.obs"))
    })

    donor_rankings[[as.character(target_id)]] <- potential_donor_ids[order(cor_vals, decreasing = TRUE)]
    setTxtProgressBar(pb1, i)
  }
  close(pb1)

  # --- 3. PHASE 2: nMSE OPTIMIZATION ---
  # Searching for the donor pool size that minimizes pre-treatment prediction error
  message(sprintf("\n[Phase 2/3] Searching for optimal donor pool size in range %d to %d...",
                  min(donor_range), max(donor_range)))

  get_summed_nmse <- function(n_donors) {
    total_nmse <- 0
    for (target_id in treated_unit_ids) {
      t_treat <- min(df[[time_var]][df[[id_var]] == target_id & df[[treatment_indicator]] == 1])
      pre_period <- sort(unique(df[[time_var]][df[[time_var]] < t_treat]))
      selected_donors <- donor_rankings[[as.character(target_id)]][1:n_donors]

      dp <- Synth::dataprep(
        foo = df, predictors = predictor_vars, predictors.op = "mean",
        time.predictors.prior = pre_period, dependent = outcome_var,
        unit.variable = id_var, unit.names.variable = name_var,
        time.variable = time_var, treatment.identifier = target_id,
        controls.identifier = selected_donors, time.optimize.ssr = pre_period,
        time.plot = unique(df[[time_var]])
      )
      synth_res <- Synth::synth(data.prep.obj = dp, verbose = FALSE)

      y_actual <- dp$Y1plot
      y_synth <- dp$Y0plot %*% synth_res$solution.w

      # Calculate Normalized Mean Square Error (nMSE)
      unit_nmse <- mean((y_actual - y_synth)^2) / var(y_actual)
      total_nmse <- total_nmse + unit_nmse
    }
    return(total_nmse)
  }

  # Perform optimization across the specified range
  if (length(donor_range) > 1) {
    if (parallel) {
      message("Parallel processing enabled via future.apply. That might take a while...")
      nmse_results <- future.apply::future_lapply(donor_range, get_summed_nmse)
    } else {
      pb2 <- txtProgressBar(min = 0, max = length(donor_range), style = 3)
      nmse_results <- list()
      for(j in 1:length(donor_range)) {
        nmse_results[[j]] <- get_summed_nmse(donor_range[j])
        setTxtProgressBar(pb2, j)
      }
      close(pb2)
    }

    nmse_vec <- unlist(nmse_results)
    optimal_n <- donor_range[which.min(nmse_vec)]

    # Create the diagnostic fit table for the elbow plot
    fit_table <- data.frame(
      n_donors = donor_range,
      summed_nmse = nmse_vec
    )

    message(sprintf("Optimal donor pool size identified: n = %d", optimal_n))
  } else {
    optimal_n <- donor_range
    fit_table <- data.frame(n_donors = optimal_n, summed_nmse = get_summed_nmse(optimal_n))
  }

  # --- 4. PHASE 3: FINAL SCM COMPUTATION ---
  # Compute final counterfactuals using the optimal donor pool size
  message(sprintf("\n[Phase 3/3] Computing final SCM models for %d units...", n_treated))

  final_synth_models <- list()
  final_dataprep_objects <- list()
  final_treatment_times <- numeric(n_treated)

  pb3 <- txtProgressBar(min = 0, max = n_treated, style = 3)

  for (k in 1:n_treated) {
    target_id <- treated_unit_ids[k]
    t_treat <- min(df[[time_var]][df[[id_var]] == target_id & df[[treatment_indicator]] == 1])
    final_treatment_times[k] <- t_treat
    pre_period <- sort(unique(df[[time_var]][df[[time_var]] < t_treat]))
    selected_donors <- donor_rankings[[as.character(target_id)]][1:optimal_n]

    dp_final <- Synth::dataprep(
      foo = df, predictors = predictor_vars, predictors.op = "mean",
      time.predictors.prior = pre_period, dependent = outcome_var,
      unit.variable = id_var, unit.names.variable = name_var,
      time.variable = time_var, treatment.identifier = target_id,
      controls.identifier = selected_donors, time.optimize.ssr = pre_period,
      time.plot = unique(df[[time_var]])
    )

    final_dataprep_objects[[k]] <- dp_final
    final_synth_models[[k]] <- Synth::synth(data.prep.obj = dp_final, verbose = FALSE)

    setTxtProgressBar(pb3, k)
  }
  close(pb3)

  message("\nSuccess: SCM panel preparation complete.\n")

  return(list(
    optimal_n = optimal_n,
    fit_table = fit_table,           # New: Fit diagnostics table
    models = final_synth_models,
    data_objects = final_dataprep_objects,
    times = final_treatment_times,
    treated_ids = treated_unit_ids
  ))
}

