#' Quick (gg)plot for RDA
#'
#' @param rda A rda object to plot, typically created by vegan::rda.
#' @param factor Either a factor or a character string of the same length than 
#'     the number of 'sites' in the rda object. Will be used as a grouping 
#'     factor for colors, hulls and ellipses.
#' @param factor.name Name of the grouping factor (in plot legend). Defaults 
#'     to 'Group' if not provided.
#' @param axis Which combination of rda axes to be chosen for the plot. May 
#'     typically be c("RAD1", "RDA2"), c("RAD1", "PC1"), or c("PC1", "PC2").
#' @param display What to display on the plot. Can include 'sites', 'species', 
#'     and 'arrows'.
#' @param hull Logical. Whether or not to draw hulls around the groups.
#' @param ellipses Logical. Whether or not to draw ellipses around the groups.
#' @param colors Set of colors to be used for the different groups. Must
#'     be of the same length than the number of groups.
#' @param arrows.col Color to use to draw the arrows (if display == 'arrows').
#' @param plot Logical. Whether or not to draw the plot. If FALSE, the plot 
#'     is stored, but not displayed. It is recommended for large data.
#' @param components.plot Logical. If TRUE, the function additionally returns
#'     a barchart of the relative importance of each axis, in percent.
#' @param component.cut Percentage under which axes will be discarded from the 
#'     component plot. Usefull when the rda returns many axes. Default to 1. 
#'     Set to 0 to display all the axes.
#' @param col.component Colors to be used in the component plot. The first 
#'     color will be used for RDA axes, and the second for PC axes.
#' @param lab.component Size of the component plot label.
#' @param alpha.text Transparency of the text (if display == 'species').
#' @param lty.ellipse Type of lines for the ellipses (see graphics::par).
#' @param alpha.ellipse Transparency of the ellipses.
#' @param size.text Size of the text (if display == 'species').
#' @param alpha.point Transparency of the points (if display == 'sites').
#' @param size.point Size of the points (if display == 'sites').
#' @param alpha.hull Transparency of the hulls.
#' @param vjust.label Adjust the vertical position of the label (displaying the
#'     importance of the selected axes in percent). 
#' @param hjust.label Adjust the horizontal position of the label (displaying the
#'     importance of the selected axes in percent). 
#' @param size.label Size of the label.
#' @param legend.position Where to place the legend (see ggplot2::theme).
#' @param size.legend Size of the text in legend.
#' @param size.axis.comp Size of the x-axis text in the component plot.
#' @param ... Additional arguments (not currently used).
#'
#' @return A ggplot object containing the requested layers to plot the
#'     results of the rda ("RDA"), and (if components.plot == TRUE) an
#'     additional ggplot object containing the component plot 
#'     ("Components")
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
#' Var1_random <- rnorm(100, mean = 132, sd = 23.7)
#' Var2_random <- rnorm(100, mean = 0.03, sd = 0.001)                  
#' Var3_random <- rnorm(100, mean = 3.78, sd = 1.556)
#'
#' Env <- data.frame("Var1" = Var1_random,
#'                  "Var2" = Var2_random,
#'                  "Var3" = Var3_random)      
#' 
#' 
#' C <- as.matrix(vegan::decostand(Community[,-1], method = "hellinger"))
#' E <- as.matrix(vegan::decostand(Env, method = "standardize"))
#' Community_RDA <- vegan::rda(C ~ E)                                   
#' 
#' Community_group <- c(rep("1", times = 33),
#'                      rep("2", times = 34),
#'                      rep("3", times = 33))
#'                      
#'  GGRDA_Community <- Visualise_RDA(Community_RDA, factor = Community_group, 
#'                                   factor.name = "Community type",
#'                                   axis = c("RDA1", "PC1"),
#'                                   display = c("sites", "species", "arrows"),
#'                                   components.plot = FALSE,
#'                                   hull = TRUE, ellipses = TRUE,
#'                                   colors = colorspace::sequential_hcl(3, palette = "Turku"),
#'                                   plot = FALSE, alpha.text = 0.5,
#'                                   size.text = 2, alpha.point = 0.80, size.point = 1,
#'                                   alpha.hull = 0.30, label.x.offset = .2, label.y.offset = .5,
#'                                   size.label = 2.3, stress.digits = 3, legend.position = "right",
#'                                   size.legend = 12)
#' GGRDA_Community
#' }
Visualise_RDA <- function(rda, factor = NULL, factor.name = NULL, axis = c("RDA1", "RDA2"),
                          display = c("sites", "species", "arrows"), hull = FALSE, ellipses = FALSE,
                          colors = NULL,  arrows.col = "red", plot = FALSE, 
                          components.plot = FALSE, component.cut = 1, 
                          col.component = c("darkred", "grey20"), lab.component = 3,
                          alpha.text = 0.5, lty.ellipse = 2, alpha.ellipse = 0.5,
                          size.text = 2, alpha.point = 0.80, size.point = 1,
                          alpha.hull = 0.30, vjust.label = 1.1, hjust.label = 1.1,
                          size.label = 2.3, legend.position = "right",
                          size.legend = 12, size.axis.comp = 6, ...) {
  ## Sites:
  if (any(display == "sites")) {
    if (!all(axis %in% colnames(as.data.frame(summary(rda)$sites)))) {
      stop("All specified axis must appear in the RDA object")
    }
    RDA_sites <- as.data.frame(cbind(as.data.frame(summary(rda)$sites)[,axis[1]], 
                                     as.data.frame(summary(rda)$sites)[,axis[2]]))
    colnames(RDA_sites) <- axis
    if(!is.null(factor)) {
      factor <- as.factor(factor) 
      if (length(factor) != nrow(RDA_sites)) {
        stop("factor must be of the same lenght than 'rda'. 
             Is there any empty sites in the original dataframe?")
      }
      if (is.null(factor.name)) {
        factor.name = "Group"
      }
      # Add the factor column:
      RDA_sites[, factor.name] <- factor
    }
  }
  # Add the scores for species data:
  if (any(display == "species")) {
    RDA_species <- data.frame(summary(rda)$species)
    # Add a column equivalent to the row name to create species labels
    RDA_species$species <- rownames(RDA_species)
  }
  # Add environment for arrows:
  if (any(display == "arrows")) {
    RDA_env <- data.frame(summary(rda)$biplot)
    RDA_env$variables <- rownames(RDA_env)
  }
  # Build hulls:
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
      Hull[[i]] <- RDA_sites[RDA_sites[,factor.name] == i, ][grDevices::chull(RDA_sites[RDA_sites[,factor.name] == i, axis]), ]
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
  RDA_plot <- ggplot()
  # Species:
  if (any(display == "species")) {
    RDA_plot <- RDA_plot +
      geom_text_repel(data = RDA_species, aes(x = RDA_species[,1], y = RDA_species[,2], 
                                              label = species),
                      alpha = alpha.text, size = size.text) 
  }
  # Sites:
  if (any(display == "sites")) {
    RDA_plot <- RDA_plot +
      geom_point(data = RDA_sites, aes(x = RDA_sites[,1], y = RDA_sites[,2], col = RDA_sites[,factor.name]), 
                 size = size.point, alpha = alpha.point) +
      scale_color_manual(values = colors,
                         name = factor.name) 
  }
  # Arrows:
  if (any(display == "arrows")) {
    RDA_plot <- RDA_plot +
      geom_segment(data = RDA_env, aes(x = 0, y = 0, xend = RDA_env[,1], yend = RDA_env[,2]), 
                   col = arrows.col) +
      geom_text(data = RDA_env, aes(x = RDA_env[,1] + .05, y = RDA_env[,2] + .05, label = variables),
                      size = 2, col = "red")
  }
  # Hulls:
  if (hull == TRUE) { 
    RDA_plot <-   RDA_plot + 
      geom_polygon(data = hull.data, aes(x = hull.data[,1], y = hull.data[,2], 
                                         fill = hull.data[,factor.name], 
                                         group = hull.data[,factor.name]), 
                   alpha = alpha.hull) +
      scale_fill_manual(values = colors,
                        name = factor.name) 
  } 
  if (ellipses == TRUE) {
    RDA_plot <-   RDA_plot + 
      stat_ellipse(data = RDA_sites, aes(x = RDA_sites[,1], y = RDA_sites[,2],
                                           col = RDA_sites[,factor.name]), 
                   lty = lty.ellipse, alpha = alpha.ellipse, level = 0.95) 
  }
  # Graphical:
  RDA_plot <-   RDA_plot + 
    annotate(geom = "label", x = Inf, y = Inf, 
             vjust = vjust.label, hjust = hjust.label,
             size = size.label,
             label = paste("Axes sores: \n", axis[1], 
                           round(summary(rda)$cont$importance[2, axis[1]]*100, digits = 2), "% \n",
                           axis[2], round(summary(rda)$cont$importance[2,axis[2]]*100, digits = 2), "%")) +
    xlab(axis[1]) +
    ylab(axis[2]) +
    theme_bw() +
    theme(legend.position = legend.position,
          text = element_text(size = size.legend))
  
  # Component barchart:
  if (components.plot) {
    GGdata <- data.frame(summary(rda)$cont$importance[2,]*100)
    GGdata$axis <- ordered(rownames(GGdata), levels = rownames(GGdata))
    GGdata <- GGdata[GGdata[,1] >= component.cut,]
    
    Component_Bchart <- ggplot() +
      geom_col(data = GGdata, aes(x = axis, y = GGdata[,1], fill = axis),
               position = "dodge", colour = NA) +
      scale_fill_manual(values = ifelse(substr(GGdata$axis, 1, 1) == "R",
                                        col.component[1], col.component[2])) +
      ggtitle("Importance of components") +
      ylab("% variation") +
      annotate(geom = "label", x = Inf, y = Inf, size = lab.component, vjust = 1.2, hjust = 1.1,
               label = paste("RDA total =", round(sum(GGdata[substr(GGdata$axis, 1, 1) == "R", 1]), digits = 2),
                             "% \n PC total =", round(sum(GGdata[substr(GGdata$axis, 1, 1) == "P", 1]), digits = 2), "%")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, size = size.axis.comp)) +
      theme(legend.position = "none")
  }
  
  #Return lists:
  if (plot) {
    print(RDA_plot)
    if (components.plot) {
      print(Component_Bchart)
    }
  }
  if (components.plot == FALSE) {
    return(RDA_plot)
  }
  if (components.plot) {
    list_output <- list()
    list_output[["RDA"]] <- RDA_plot
    list_output[["Components"]] <- Component_Bchart
    return(list_output)
  }
}
