#' Plot Heterogeneous Treatment Effects from Synthetic Control Models
#'
#' Visualizes individual treatment effects, heterogeneous average treatment
#' effects, and highest density regions over relative time.
#'
#' @param hatt_output Output from \code{estimate_hatt()}.
#' @param title Plot title.
#' @param hdr_label Label for the HDR region.
#' @param hatt_label Label for the HATT lines.
#' @param ite_label Label for individual treatment effect lines.
#' @param color_hdr Color of HDR ribbons.
#' @param color_ite Color of individual treatment effect lines.
#' @param color_hatt Color of HATT lines.
#' @param color_ref Color of reference line at zero.
#' @param color_treat Color of treatment time line.
#'
#' @return A ggplot object.
#'
#' @export



plot_scm_hatt <- function(hatt_output,
                          title = "Heterogeneous Average Treatment Effects",
                          hdr_label = "HDR Region",
                          hatt_label = "Estimated heterogeneous ATT",
                          ite_label = "",
                          color_hdr   = "black",
                          color_ite   = "#66C2A5",
                          color_hatt  = "#E78AC3",
                          color_ref   = "black",
                          color_treat = "#FC8D62") {


  # --- 1. DATA PREPARATION (Connection Logic) ---
  ite_data  <- hatt_output$results_df
  hdr_raw   <- hatt_output$hdr_df
  hatt_raw  <- hatt_output$hatt_df

  all_times <- unique(hatt_raw$relative_time)
  all_clusters <- unique(hatt_raw$cluster_id)
  full_grid <- expand.grid(relative_time = all_times, cluster_id = all_clusters)

  hatt_data <- full_grid %>%
    left_join(hatt_raw, by = c("relative_time", "cluster_id")) %>%
    group_by(relative_time) %>%
    fill(hatt_value, .direction = "downup") %>%
    ungroup()

  hdr_data <- full_grid %>%
    left_join(hdr_raw, by = c("relative_time", "cluster_id")) %>%
    group_by(relative_time) %>%
    fill(lower_bound, upper_bound, .direction = "downup") %>%
    ungroup()

  # --- 2. PLOTTING WITH LEGEND MAPPING ---
  p <- ggplot() +
    # Layer 1: HDR Regions (Ribbons) mapped to 'fill' for legend
    geom_ribbon(data = hdr_data,
                aes(x = relative_time, ymin = lower_bound, ymax = upper_bound,
                    group = cluster_id, fill = "HDR"),
                alpha = 0.3) +

    # Layer 2: ITE Lines mapped to 'color' for legend
    geom_line(data = ite_data,
              aes(x = relative_time, y = ite_gap, group = unit, color = "ITE"),
              size = 0.4, alpha = 0.4) +

    # Layer 3: HATT Lines mapped to 'color' for legend
    geom_line(data = hatt_data,
              aes(x = relative_time, y = hatt_value, group = cluster_id, color = "HATT"),
              size = 1.1) +

    # --- 3. REFERENCE LINES ---
    geom_hline(yintercept = 0, linetype = "solid", color = color_ref) +
    geom_vline(xintercept = 0, linetype = "dashed", color = color_treat) +

    # --- 4. MANUAL SCALES FOR LEGEND ---
    # Setting the labels as used in your original script
    scale_fill_manual(name = NULL,
                      values = c("HDR" = color_hdr),
                      labels = c("HDR" = hdr_label)) +
    scale_color_manual(name = NULL,
                       values = c("ITE" = color_ite, "HATT" = color_hatt),
                       labels = c("ITE" = hatt_label,
                                  "HATT" = hatt_label)) +

    # --- 5. THEME AND DESIGN ---
    theme_minimal() +
    labs(
      title = title,
      x = "Time relative to treatment time",
      y = "Effect"
    ) +
    theme(
      legend.position = "top", # Legend on top as per academic standard
      legend.direction = "horizontal",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(size = 11)
    ) +
    # Guiding the legend to display lines and boxes correctly
    guides(
      color = guide_legend(override.aes = list(size = c(1.1, 0.4))),
      fill = guide_legend(override.aes = list(alpha = 0.3))
    )

  return(p)
}
