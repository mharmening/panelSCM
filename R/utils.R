# Global variables to avoid R CMD check NOTEs
# This file should be saved as R/utils.R

#' @importFrom stats cor density var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom dplyr left_join group_by ungroup %>%
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon geom_hline 
#'   geom_vline scale_x_continuous scale_fill_manual scale_color_manual 
#'   theme_classic theme_minimal labs theme element_blank element_text 
#'   guides guide_legend
#' @importFrom tidyr fill
NULL

# Declare global variables used in NSE (non-standard evaluation) from dplyr/ggplot2
utils::globalVariables(c(
  "n_donors",
  "summed_nmse", 
  "relative_time",
  "hatt_value",
  "lower_bound",
  "upper_bound",
  "cluster_id",
  "ite_gap",
  "unit"
))
