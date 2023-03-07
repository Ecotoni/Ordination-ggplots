#' Quick (gg)plot for nMDS objects:
#'
#' @param nMDS A nMDS object to plot, typically created by vegan::metaMDS.
#' @param factor Either a factor or a character string of the same length than 
#'     the number of 'sites' in the nMDS object. Will be used as a grouping 
#'     factor for colors, hulls and ellipses.
#' @param factor.name Name of the grouping factor (in plot legend). Defaults 
#'     to 'Group' if not provided.
#' @param display What to display on the plot. Either 'sites', 'species' 
#'     or both.
#' @param facet.factor An additional factor (or character string) to 
#'     create facets. Must be of the same length than the number of 'sites' 
#'     in the nMDS object.
#' @param facet.name Name of the faceting factor. Defaults 
#'     to 'Facet' if not provided.
#' @param scales Scale dependancy between the facets. Can be "fixed", "free_x",
#'     "free_y" or "free".
#' @param hull Logical. Whether or not to draw hulls around the groups.
#' @param ellipses Logical. Whether or not to draw ellipses around the groups.
#' @param colors Set of colors to be used for the different groups. Must
#'     be of the same length than the number of groups.
#' @param plot Logical. Whether or not to draw the plot. If FALSE, the plot 
#'     is stored, but not displayed. It is recomended for large data.
#' @param alpha.text Transparency of the text (if display == 'species').
#' @param lty.ellipse Type of lines for the ellipses (see graphics::par).
#' @param alpha.ellipse Transparency of the ellipses.
#' @param size.text Size of the text (if display == 'species').
#' @param alpha.point Transparency of the points (if display == 'sites').
#' @param size.point Size of the points (if display == 'sites').
#' @param alpha.hull Transparency of the hulls.
#' @param vjust.label Adjust the vertical position of the label (displaying the
#'     stress value). 
#' @param hjust.label Adjust the horizontal position of the label (displaying the
#'     stress value). 
#' @param size.label Size of the label (displaying the stress value). 
#' @param stress.digits Number of digits to display for the stress value.
#' @param legend.position Where to place the legend (see ggplot2::theme)
#' @param size.legend Size of the text in legend.
#' @param ... Additional arguments (not currently used).
#'
#' @return A ggplot object containing the requested layers.
#' @export
#'
#' @examples
#' \dontrun{
#' require(EcotoneFinder)
#' require(ggplot2)
#' require(ggrepel)
#' 
#' Community <- EcotoneFinder::SyntheticData(SpeciesNum = 100, CommunityNum = 5, Length = 100,
#'                                           Parameters = list(a = c(60, 20, 50, 40, 60),
#'                                                             b = c(0, 0, 50, 100, 100),
#'                                                             c = c(0.05, 0.12, 0.14, 0.3, 0.05)),
#'                                           dev.a = 30, dev.b = 5, dev.c = 0.01,
#'                                           pal = colorspace::sequential_hcl(5, palette = "Turku"))
#'                                           
#' Community_nMDS <- vegan::metaMDS(Community[,-1], 
#'                                  distance = "bray", k = 2, plot = F,
#'                                  autotransform = TRUE,  wascores = TRUE,
#'                                  trymax = 100)
#'                                  
#' Community_group <- c(rep("1", times = 33),
#'                      rep("2", times = 34),
#'                      rep("3", times = 33))
#'
#' GGnMDS_Community <- Visualise_nMDS(Community_nMDS, factor = Community_group, factor.name = "Community type",
#'                                   display = c("sites", "species"),
#'                                   facet.factor = NULL, facet.name = NULL, 
#'                                   hull = TRUE, ellipses = TRUE,
#'                                   colors = colorspace::sequential_hcl(3, palette = "Turku"),
#'                                   plot = FALSE, alpha.text = 0.5,
#'                                   size.text = 2, alpha.point = 0.80, size.point = 1,
#'                                   alpha.hull = 0.30, label.x.offset = .2, label.y.offset = .5,
#'                                   size.label = 2.3, stress.digits = 3, legend.position = "right",
#'                                   size.legend = 12)
#' GGnMDS_Community
#' }
Visualise_nMDS <- function(nMDS, factor = NULL, factor.name = NULL, display = c("sites", "species"),
                           facet.factor = NULL, facet.name = NULL, scales = "free",
                           hull = FALSE, ellipses = FALSE,
                           colors = NULL,
                           plot = TRUE, alpha.text = 0.5, lty.ellipse = 2, alpha.ellipse = 0.5,
                           size.text = 2, alpha.point = 0.80, size.point = 1,
                           alpha.hull = 0.30, vjust.label = 1.1, hjust.label = 1.1,
                           size.label = 2.3, stress.digits = 3, legend.position = "right",
                           size.legend = 12, ...) {
  # Extract site scores:
  if (any(display == "sites")) {
    nMDS_scores <- as.data.frame(vegan::scores(nMDS, "site"))
    if(!is.null(factor)) {
      factor <- as.factor(factor) 
      if (length(factor) != nrow(nMDS_scores)) {
        stop("factor must be of the same lenght than nMDS. 
             Is there any empty sites in the original dataframe?")
      }
      if (is.null(factor.name)) {
        factor.name = "Group"
      }
      # Add the factor column:
      nMDS_scores[, factor.name] <- factor
    }
    if (!is.null(facet.factor)) {
      if (length(facet.factor) != nrow(nMDS_scores)) {
        stop("facet.factor must be of the same lenght than nMDS.")
      }
      if (is.null(facet.name)) {
        facet.name = "Facet"
      }
      # Add the facet column:
      nMDS_scores[, facet.name] <- facet.factor
    }
  }
  # Add the scores for species data:
  if (any(display == "species")) {
    nMDS_species_scores <- as.data.frame(vegan::scores(nMDS, "species"))
    nMDS_species_scores <- na.omit(nMDS_species_scores)
    # Add a column equivalent to the row name to create species labels
    nMDS_species_scores$species <- rownames(nMDS_species_scores)
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
    if (is.null(facet.factor)) {
      Hull <- list()
      for (i in levels(factor)) {
        Hull[[i]] <- nMDS_scores[nMDS_scores[,factor.name] == i, ][grDevices::chull(nMDS_scores[nMDS_scores[,factor.name] == i, 
                                       c("NMDS1", "NMDS2")]), ]
      }
      hull.data <- purrr::reduce(Hull, rbind)
    } else {
      facet.factor <- as.factor(facet.factor)
      Hull <- list()
      hull.data <- list()
      for (j in levels(facet.factor)) {
        for (i in levels(factor)) {
          Hull[[j]][[i]] <- nMDS_scores[nMDS_scores[,factor.name] == i & 
                                          nMDS_scores[, facet.name] == j, ][grDevices::chull(nMDS_scores[nMDS_scores[,factor.name] == i &
                                                                                                           nMDS_scores[, facet.name] == j, 
                                                                                                  c("NMDS1", "NMDS2")]), ]
        }
        hull.data[[j]] <- purrr::reduce(Hull[[j]], rbind)
        hull.data[[j]][, facet.name] <- j
      }
      hull.data <- purrr::reduce(hull.data, rbind)
    }
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
  nMDS_plot <- ggplot()
  # Species:
  if (any(display == "species")) {
    nMDS_plot <- nMDS_plot +
      geom_text_repel(data = nMDS_species_scores, aes(x = NMDS1, y = NMDS2, label = species),
                      alpha = alpha.text, size = size.text) 
  }
  # Sites:
  if (any(display == "sites")) {
    nMDS_plot <- nMDS_plot +
      geom_point(data = nMDS_scores, aes(x = NMDS1, y = NMDS2, col = nMDS_scores[,factor.name]), 
                 size = size.point, alpha = alpha.point) +
      scale_color_manual(values = colors,
                         name = factor.name) 
  }
   # Faceting:
  if (!is.null(facet.factor)) {
    nMDS_plot <-   nMDS_plot + 
      facet_wrap(facet.name, scales = scales) +
      labs(subtitle = paste("Stress: ", round(nMDS$stress, digits = stress.digits))) 
  }
  if (is.null(facet.factor)) {
    nMDS_plot <-   nMDS_plot + 
      annotate(geom = "label", x = Inf, y = Inf, size = size.label,
               vjust = vjust.label, hjust = hjust.label,
               label = paste("Stress: ", round(nMDS$stress, digits = stress.digits))) 
  }
    # Add Hulls
    if (hull == TRUE) { 
      nMDS_plot <-   nMDS_plot + 
      geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, 
                                         fill = hull.data[,factor.name], 
                                         group = hull.data[,factor.name]), 
                   alpha = alpha.hull) +
        scale_fill_manual(values = colors,
                          name = factor.name) 
    } 
  if (ellipses == TRUE) {
    nMDS_plot <-   nMDS_plot + 
    stat_ellipse(data = nMDS_scores, aes(x = NMDS1, y = NMDS2,
                                            col = nMDS_scores[,factor.name]), 
                 lty = lty.ellipse, alpha = alpha.ellipse, level = 0.95) 
  }
  # Graphical:
  nMDS_plot <-   nMDS_plot + 
    theme_bw() +
    theme(legend.position = legend.position,
          text = element_text(size = size.legend),
          plot.subtitle = element_text(size = size.label))
  # Return lists:
    if (plot) {
      print(nMDS_plot)
    }
    return(nMDS_plot)
}
