#' Plot Donor Pool Optimization Fit
#'
#' Plots the summed normalized mean squared error used to select the optimal
#' donor pool size.
#'
#' @param scm_output Output from \code{prep_scm_panel()}.
#' @param line_color Line color.
#' @param dot_color Point color.
#'
#' @return A ggplot object.
#'
#' @export


plot_scm_fit <- function(scm_output, line_color = "black", dot_color = "black") {

  ggplot(scm_output$fit_table, aes(x = n_donors, y = summed_nmse)) +
    geom_line(color = line_color, size = 0.8) +
    geom_point(color = dot_color, size = 3) +
    scale_x_continuous(
      breaks = scm_output$fit_table$n_donors
    ) +
    theme_classic() +
    labs(
      title = "Model Fit Optimization",
      subtitle = paste("Optimal Size:", scm_output$optimal_n),
      x = "Number of Units in Donor Pool",
      y = "Summed nMSE (Pre-treatment)"
    )
}
