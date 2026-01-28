#' Estimate Heterogeneous Average Treatment Effects
#'
#' Estimates heterogeneous average treatment effects (HATTs) from individual
#' treatment effects using kernel density estimation and highest density regions.
#'
#' @param scm_results Output from \code{prep_scm_panel()}.
#' @param pre_period_length Number of pre-treatment periods.
#' @param post_period_length Number of post-treatment periods.
#' @param hdr_conf_level Confidence level for the HDR visualization.
#' @param max_hatt_conf_search Maximum confidence level searched for HATT estimation.
#' @param bw_hdr_visual Bandwidth parameter for HDR visualization.
#' @param bw_hatt_estimate Bandwidth parameter for HATT estimation.
#'
#' @return A list containing individual treatment effects, estimated HATTs,
#' HDR bounds, and the optimized confidence level.
#'
#' @export




estimate_hatt <- function(scm_results,
                          pre_period_length = 10,
                          post_period_length = 11,
                          hdr_conf_level = 0.90,
                          max_hatt_conf_search = 0.35,
                          bw_hdr_visual = 0.8,
                          bw_hatt_estimate = 6) {


  # Initialize dimensions based on SCM input
  n_units <- length(scm_results$treated_ids)
  total_periods <- pre_period_length + post_period_length

  # --- STEP A: ALIGN INDIVIDUAL TREATMENT EFFECTS (ITE) ---
  # Creating a matrix for internal math and a tidy list for output data
  ite_matrix <- matrix(NA, nrow = n_units, ncol = total_periods)
  tidy_ite_list <- list()

  for (i in 1:n_units) {
    # Extract data objects for the specific treated unit
    data_obj <- scm_results$data_objects[[i]]
    model_obj <- scm_results$models[[i]]
    t_treat <- scm_results$times[i]

    # Gap calculation: Actual Y minus Synthetic Y
    gap <- as.numeric(data_obj$Y1plot) - (as.numeric(data_obj$Y0plot %*% model_obj$solution.w))

    # Align by relative time centered on the individual treatment date
    start_idx <- t_treat - pre_period_length
    end_idx <- t_treat + post_period_length - 1
    ite_matrix[i, ] <- gap[start_idx:end_idx]

    # Append to list for the final results_df (the "spaghetti" lines)
    tidy_ite_list[[i]] <- data.frame(
      unit = scm_results$treated_ids[i],
      relative_time = (-pre_period_length + 1):post_period_length,
      ite_gap = ite_matrix[i, ]
    )
  }
  results_df <- do.call(rbind, tidy_ite_list)

  # --- STEP B: OPTIMIZE HATT LEVEL (The "j-loop") ---
  # Iteratively searching for the HDR confidence level that maximizes
  # the correlation between raw ITEs and their assigned cluster means.
  test_levels <- seq(0.01, 0.99, 0.01)
  search_limit <- round(max_hatt_conf_search * 100)
  store_correlation <- numeric(search_limit)

  for (j in 1:search_limit) {
    conf_j <- test_levels[j]
    assigned_means <- matrix(NA, nrow = n_units, ncol = total_periods)

    for (t in 1:total_periods) {
      # KDE for current time point
      dens <- density(ite_matrix[, t], bw = bw_hatt_estimate)
      sorted_y <- sort(dens$y, decreasing = TRUE)
      threshold <- sorted_y[which(cumsum(sorted_y)/sum(sorted_y) >= conf_j)[1]]

      # Identify HDR coordinates
      hpd_x <- dens$x[dens$y >= threshold]

      # Determine cluster separation based on the modus of x-differences
      modus <- as.numeric(names(sort(-table(diff(hpd_x)))[1])) + 0.001
      starts <- hpd_x[c(1, which(diff(hpd_x) > modus) + 1)]
      ends <- hpd_x[c(which(diff(hpd_x) > modus), length(hpd_x))]

      # Calculate representative cluster means (HATTs)
      means_t <- sapply(seq_along(starts), function(k) {
        idx <- dens$x >= starts[k] & dens$x <= ends[k]
        # Weighted mean calculation: (y * x) / sum(y)
        sum(dens$y[idx] * dens$x[idx]) / sum(dens$y[idx])
      })

      # Assign each unit's ITE to the closest identified cluster mean
      for (u in 1:n_units) {
        assigned_means[u, t] <- means_t[which.min(abs(ite_matrix[u, t] - means_t))]
      }
    }
    # Correlate assigned group means with original individual effects
    store_correlation[j] <- cor(as.vector(assigned_means), as.vector(ite_matrix))
  }

  # Select the level with the highest correlation score
  opt_conf_hatt <- test_levels[which.max(store_correlation)]

  # --- STEP C: FINAL CALCULATION ---
  hatt_list <- list()
  hdr_list <- list()

  for (t in 1:total_periods) {
    t_rel <- (-pre_period_length + t)

    # 1. HATT (Red Lines) - Uses optimized confidence level
    dens_hatt <- density(ite_matrix[, t], bw = bw_hatt_estimate)
    s_y <- sort(dens_hatt$y, decreasing = TRUE)
    t_h <- s_y[which(cumsum(s_y)/sum(s_y) >= opt_conf_hatt)[1]]
    h_x <- dens_hatt$x[dens_hatt$y >= t_h]
    mod_h <- as.numeric(names(sort(-table(diff(h_x)))[1])) + 0.001
    st_h <- h_x[c(1, which(diff(h_x) > mod_h) + 1)]
    en_h <- h_x[c(which(diff(h_x) > mod_h), length(h_x))]

    means_final <- sapply(seq_along(st_h), function(k) {
      idx <- dens_hatt$x >= st_h[k] & dens_hatt$x <= en_h[k]
      sum(dens_hatt$y[idx] * dens_hatt$x[idx]) / sum(dens_hatt$y[idx])
    })
    hatt_list[[t]] <- data.frame(relative_time = t_rel, hatt_value = means_final, cluster_id = as.factor(1:length(means_final)))

    # 2. HDR (Gray Ribbons) - Uses fixed visual level and scaled bandwidth
    dens_hdr <- density(ite_matrix[, t], bw = 10^bw_hdr_visual)
    s_y_r <- sort(dens_hdr$y, decreasing = TRUE)
    t_h_r <- s_y_r[which(cumsum(s_y_r)/sum(s_y_r) >= hdr_conf_level)[1]]
    h_x_r <- dens_hdr$x[dens_hdr$y >= t_h_r]
    mod_r <- as.numeric(names(sort(-table(diff(h_x_r)))[1])) + 0.001
    st_r <- h_x_r[c(1, which(diff(h_x_r) > mod_r) + 1)]
    en_r <- h_x_r[c(which(diff(h_x_r) > mod_r), length(h_x_r))]

    hdr_list[[t]] <- data.frame(relative_time = t_rel, lower_bound = st_r, upper_bound = en_r, cluster_id = as.factor(1:length(st_r)))
  }

  # Clean up row names and merge lists into master dataframes
  hatt_df <- do.call(rbind, hatt_list)
  hdr_df <- do.call(rbind, hdr_list)
  rownames(hatt_df) <- NULL
  rownames(hdr_df) <- NULL

  # Return everything in a structured list for Step 3 (Plotting)
  return(list(
    results_df = results_df,
    hatt_df = hatt_df,
    hdr_df = hdr_df,
    optimized_hatt_level = opt_conf_hatt
  ))
}
