#' Quick (gg)plot for DCA
#'
#' @param dca A dca object to plot, typically created by vegan::decorana.
#' @param factor Either a factor or a character string of the same length than 
#'     the number of 'sites' in the dca object. Will be used as a grouping 
#'     factor for colors, hulls and ellipses.
#' @param factor.name Name of the grouping factor (in plot legend). Defaults 
#'     to 'Group' if not provided.
#' @param display What to display on the plot. Either 'sites', 'species' 
#'     or both.
#' @param hull Logical. Whether or not to draw hulls around the groups.
#' @param ellipses Logical. Whether or not to draw ellipses around the groups.
#' @param plot Logical. Whether or not to draw the plot. If FALSE, the plot 
#'     is stored, but not displayed. It is recommended for large data.
#' @param colors Set of colors to be used for the different groups. Must
#'     be of the same length than the number of groups.
#' @param alpha.text Transparency of the text (if display == 'species').
#' @param lty.ellipse Type of lines for the ellipses (see graphics::par).
#' @param alpha.ellipse Transparency of the ellipses.
#' @param size.text Size of the text (if display == 'species').
#' @param alpha.point Transparency of the points (if display == 'sites').
#' @param size.point Size of the points (if display == 'sites').
#' @param alpha.hull Transparency of the hulls.
#' @param vjust.label Adjust the vertical position of the label (displaying the
#'     total inertia and the length of the first DCA axis). 
#' @param hjust.label Adjust the horizontal position of the label (displaying the
#'     total inertia and the length of the first DCA axis). 
#' @param size.label Size of the label.
#' @param legend.position Where to place the legend (see ggplot2::theme).
#' @param size.legend Size of the text in legend.
#' @param ... Additional arguments (not currently used).
#'
#' @return A ggplot object containing the requested layers.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' Community <- EcotoneFinder::SyntheticData(SpeciesNum = 100, CommunityNum = 5, Length = 100,
#'                                           Parameters = list(a = c(60, 20, 50, 40, 60),
#'                                                             b = c(0, 0, 50, 100, 100),
#'                                                             c = c(0.05, 0.12, 0.14, 0.3, 0.05)),
#'                                           dev.a = 30, dev.b = 5, dev.c = 0.01,
#'                                           pal = colorspace::sequential_hcl(5, palette = "Turku"))
#' 
#' Comm_dat <- vegan::decostand(Community[,-1], method = "hellinger")
#' Community_DCA <- vegan::decorana(Community[,-1])
#' 
#' Community_group <- c(rep("1", times = 33),
#'                      rep("2", times = 34),
#'                      rep("3", times = 33))
#' 
#' GGDCA_Community <- Visualise_DCA(Community_DCA, factor = Community_group, 
#'                                  factor.name = "Community type",
#'                                  display = c("sites"),
#'                                  facet.factor = NULL, facet.name = NULL, 
#'                                  hull = TRUE, ellipses = TRUE,
#'                                  colors = colorspace::sequential_hcl(3, palette = "Turku"),
#'                                  plot = FALSE, alpha.text = 0.5,
#'                                  size.text = 2, alpha.point = 0.80, size.point = 1,
#'                                  alpha.hull = 0.30, label.x.offset = .2, label.y.offset = .5,
#'                                  size.label = 2.3, stress.digits = 3, legend.position = "right",
#'                                  size.legend = 12)
#' GGDCA_Community
#' 
#' 
#' }
Visualise_DCA <- function(dca, factor = NULL, factor.name = NULL, 
                          display = c("sites", "species"), hull = FALSE, ellipses = FALSE,
                          plot = FALSE, colors = NULL,
                          alpha.text = 0.5, lty.ellipse = 2, alpha.ellipse = 0.5,
                          size.text = 2, alpha.point = 0.80, size.point = 1,
                          alpha.hull = 0.30, vjust.label = 1.1, hjust.label = 1.1,
                          size.label = 2.3, legend.position = "right",
                          size.legend = 12, ...) {
  # Extract site scores:
  if (any(display == "sites")) {
    DCA_scores <- data.frame("DCA1" = dca$rproj[,1], 
                             "DCA2" = dca$rproj[,2])
    if(!is.null(factor)) {
      factor <- as.factor(factor) 
      if (length(factor) != nrow(DCA_scores)) {
        stop("factor must be of the same lenght than dca. 
             Is there any empty sites in the original dataframe?")
      }
      if (is.null(factor.name)) {
        factor.name = "Group"
      }
      # Add the factor column:
      DCA_scores[, factor.name] <- factor
    }
  }
  # Add the scores for species data:
  if (any(display == "species")) {
    DCA_species_scores <- as.data.frame(dca$cproj)
    # Add a column equivalent to the row name to create species labels
    DCA_species_scores$species <- rownames(DCA_species_scores)
  }
  # Hulls:
  if (hull) {
    if (all(display != "sites")) {
      stop("hulls can only be drawn when 'display' contains 'sites'")
    }
    if (is.null(factor)) {
      stop("'factor' must be specified when hull = TRUE")
    }
    # Building hulls: (according to factor)
    Hull <- list()
    for (i in levels(factor)) {
      Hull[[i]] <- DCA_scores[DCA_scores[,factor.name] == i, ][grDevices::chull(DCA_scores[DCA_scores[,factor.name] == i, 
                                                                                                c("DCA1", "DCA2")]), ]
    }
    hull.data <- purrr::reduce(Hull, rbind)
  }
  ## Plot:
  require(ggplot2)
  require(ggrepel)
  if (is.null(factor) & is.null(colors)) {
    colors = "black"
  }
  if (is.null(factor) & length(colors) > 1) {
    colors = colors[1]
  }
  if (length(colors) != length(levels(factor))) {
    stop("colors should be of the same length than the grouping factor")
  }
  # Plot:
  DCA_plot <- ggplot()
  # Species:
  if (any(display == "species")) {
    DCA_plot <- DCA_plot +
      geom_text_repel(data = DCA_species_scores, aes(x = DCA1, y = DCA2, label = species),
                      alpha = alpha.text, size = size.text) 
  }
  # Sites:
  if (any(display == "sites")) {
    DCA_plot <- DCA_plot +
      geom_point(data = DCA_scores, aes(x = DCA1, y = DCA2, col = DCA_scores[,factor.name]), 
                 size = size.point, alpha = alpha.point) +
      scale_color_manual(values = colors,
                         name = factor.name) 
  }
  # Add hulls & ellipses:
  if (hull == TRUE) { 
    DCA_plot <-   DCA_plot + 
      geom_polygon(data = hull.data, aes(x = DCA1, y = DCA2, 
                                         fill = hull.data[,factor.name], 
                                         group = hull.data[,factor.name]), 
                   alpha = alpha.hull) +
      scale_fill_manual(values = colors,
                        name = factor.name) 
  } 
  if (ellipses == TRUE) {
    DCA_plot <-   DCA_plot + 
      stat_ellipse(data = DCA_scores, aes(x = DCA1, y = DCA2,
                                           col = DCA_scores[,factor.name]), 
                   lty = lty.ellipse, alpha = alpha.ellipse, level = 0.95) 
  }
  # Graphical:
  DCA_plot <-   DCA_plot + 
    annotate(geom = "label", x = Inf, y = Inf, size = size.label,
             vjust = vjust.label, hjust = hjust.label,
             label = paste("Total inertia:", round(dca$totchi, digits = 3),"\n",
               "First axis length:", round(max(dca$rproj[,1]) - min(dca$rproj[,1]), digits = 3))) +
    theme_bw() +
    theme(legend.position = legend.position,
          text = element_text(size = size.legend))
  # Return lists:
  if (plot) {
    print(DCA_plot)
  }
  return(DCA_plot)
}
