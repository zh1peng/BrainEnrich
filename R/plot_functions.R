#' Plot Brain Data
#'
#' This function creates a brain plot using the ggplot2 and ggseg packages.
#'
#' @param df2plot A data frame containing the data to plot.
#' @param ats A character string indicating the atlas to use ('dx', 'dk', 'aseg').
#' @param what2plot A character string indicating the variable to plot ('statistic').
#' @param filterby A character string indicating the filter to apply ('p.value', 'p.adj', 'none').
#' @param title2show A character string indicating the title of the plot.
#' @param limit2show A numeric vector of length 2 indicating the limits for the color scale.
#' @param legend2show A character string indicating the legend title.
#' @param hide_legend A logical value indicating whether to hide the legend.
#' @param hem A character string indicating which hemisphere to plot ('both', 'left', 'right').
#' @param low A character string indicating the color for the low end of the scale.
#' @param mid A character string indicating the color for the midpoint of the scale.
#' @param high A character string indicating the color for the high end of the scale.
#' @param sufix2remove A character string indicating the suffix to remove from labels.
#' @return A ggplot2 object.
#' @import ggplot2
#' @import ggseg
#' @importFrom dplyr filter mutate %>% recode
#' @importFrom stats as.formula
#' @export
plot_brain <- function(df2plot,
                       ats = c("dx", "dk", "aseg"),
                       what2plot = "statistic",
                       filterby = c("p.value", "p.adj", "none"),
                       title2show = "",
                       limit2show = c(-15, 15),
                       legend2show = "Stat",
                       hide_legend = FALSE,
                       hem = "both",
                       low = "steelblue1",
                       mid = "white",
                       high = "firebrick1",
                       sufix2remove = "_thickavg") {
  # Match arguments with allowed values
  ats <- match.arg(ats)
  filterby <- match.arg(filterby)

  # Apply filters based on p-values
  if (filterby == "p.value") {
    df2plot <- df2plot %>% filter(.data$p.value < 0.05)
  } else if (filterby == "p.adj") {
    df2plot <- df2plot %>% filter(.data$p.adj < 0.05)
  } else if (filterby == "none") {
    df2plot <- df2plot
  }

  # Modify labels and join data based on anatomical terms
  if (ats == "dx" || ats == "dk") {
    atlas <- ifelse(ats == "dx", "desterieux", "dk")
    atlas_data <- getExportedValue("ggseg", atlas)
    df2plot <- df2plot %>%
      mutate(
        label = sub(".*L_", "lh_", .data$label),
        label = sub(".*R_", "rh_", .data$label),
        label = sub(sufix2remove, "", .data$label)
      ) %>%
      brain_join(atlas_data) %>%
      reposition_brain(as.formula(". ~ hemi + side"))
  } else if (ats == "aseg") {
    atlas_data <- getExportedValue("ggseg", "aseg") # ggseg::aseg
    df2plot <- df2plot %>%
      mutate(label = recode(.data$key,
        "SV_L_" = "Left-",
        "SV_R_" = "Right-",
        "thal" = "Thalamus-Proper",
        "caud" = "Caudate",
        "put" = "Putamen",
        "pal" = "Pallidum",
        "hippo" = "Hippocampus",
        "amyg" = "Amygdala",
        "accumb" = "Accumben-area",
        "LatVent" = "Lateral-Ventricle"
      )) %>%
      filter(!.data$key %in% c("SV_L_accumb", "SV_R_accumb")) %>%
      brain_join(atlas_data) %>%
      filter(.data$side == "coronal")
  }

  # Apply conditional filtering based on hemisphere
  if (hem %in% c("left", "right")) {
    df2plot <- df2plot %>% filter(.data$hemi == hem)
  }

  # Construct the plot
  p <- df2plot %>%
    ggplot() +
    geom_sf(aes(fill = .data[[what2plot]])) +
    scale_fill_gradient2(midpoint = 0, low = low, mid = mid, high = high, space = "Lab", limits = limit2show) +
    ggtitle(title2show) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill = legend2show) +
    {
      if (hide_legend) theme(legend.position = "none")
    }

  return(p)
}
