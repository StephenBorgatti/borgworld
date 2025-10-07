#' Create MDS-style plot with optional clustering colors
#'
#' @description
#' Creates a scatter plot of MDS coordinates with labeled points.
#' Optionally colors points and labels by a grouping variable (e.g., clusters).
#' Automatically adjusts label positions to avoid overlap with plot edges.
#'
#' @param coords A matrix or data frame of coordinates (n x 2), or an MDS object
#'   with a \code{points} component (from cmdscale, isoMDS, etc.)
#' @param labels Character vector of point labels. If NULL, uses row names of
#'   coords or numeric indices.
#' @param color Optional factor or character vector for coloring points by group.
#'   Typically cluster assignments or categorical variables.
#' @param main Character string for plot title. Default is "MDS Plot".
#' @param cex Numeric expansion factor for point size. Default is 1.
#' @param cex.text Numeric expansion factor for label text size. Default is 0.8.
#' @param pch Plotting character (symbol) for points. Default is 19 (solid circle).
#' @param pos Position of labels relative to points. Can be:
#'   \itemize{
#'     \item "auto" (default): Smart positioning based on location in plot
#'     \item 1: Below point
#'     \item 2: Left of point
#'     \item 3: Above point
#'     \item 4: Right of point
#'   }
#' @param ... Additional graphical parameters passed to \code{plot()}.
#'
#' @return Invisibly returns a list containing:
#'   \item{coords}{The coordinate matrix used for plotting}
#'   \item{colors}{Vector of colors used for each point}
#'   \item{labels}{Vector of labels used}
#'
#' @examples
#' \dontrun{
#' # Basic MDS plot
#' mds_coords <- cmdscale(cities)
#' bplotcoord(mds_coords)
#'
#' # With clustering colors
#' bplotcoord(mds_coords, color = citiesattr$cluster)
#'
#' # Custom labels and parameters
#' bplotcoord(mds_coords,
#'           labels = citiesattr$name,
#'           color = citiesattr$region,
#'           main = "US Cities by Region",
#'           cex = 1.2,
#'           cex.text = 0.7)
#'
#' # Manual label positioning
#' bplotcoord(mds_coords,
#'           color = citiesattr$cluster,
#'           pos = 1)  # All labels below points
#' }
#'
#' @export
bplotcoord <- function(coords, labels = NULL, color = NULL,
                       main = "MDS Plot", cex = 1, cex.text = 0.8,
                       pch = 19, pos = "auto", ...) {

  # Extract coordinates if coords is an MDS object
  if(inherits(coords, c("cmdscale", "isoMDS", "metaMDS"))) {
    coord_matrix <- coords$points
  } else {
    coord_matrix <- as.matrix(coords)
  }

  # Validate dimensions
  if(ncol(coord_matrix) < 2) {
    stop("coords must have at least 2 columns")
  }

  # Set up labels
  if(is.null(labels)) {
    labels <- rownames(coord_matrix)
    if(is.null(labels)) {
      labels <- 1:nrow(coord_matrix)
    }
  }

  # Validate label length
  if(length(labels) != nrow(coord_matrix)) {
    stop("Length of labels must match number of rows in coords")
  }

  # Handle colors
  if(is.null(color)) {
    # No clustering - use black
    col_vector <- rep("black", nrow(coord_matrix))
    col_palette <- "black"
  } else {
    # Validate color length
    if(length(color) != nrow(coord_matrix)) {
      stop("Length of color must match number of rows in coords")
    }

    # Convert to factor if not already
    color_factor <- as.factor(color)
    n_colors <- length(unique(color_factor))

    # Create color palette
    if(n_colors <= 8) {
      # ColorBrewer Set1 colors for small number of groups
      col_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                       "#FF7F00", "#FFFF33", "#A65628", "#F781BF")[1:n_colors]
    } else {
      # Rainbow palette for many groups
      col_palette <- rainbow(n_colors)
    }

    col_vector <- col_palette[as.numeric(color_factor)]
  }

  # Smart positioning: adjust based on location in plot
  if(identical(pos, "auto")) {
    # Get position for each point based on its location
    y_coords <- coord_matrix[,2]
    y_range <- range(y_coords)
    y_top_third <- y_range[1] + 2 * diff(y_range)/3

    # Points in top third get labels below, others above
    pos_vector <- ifelse(y_coords > y_top_third,
                         1,  # below for points in top third
                         3)  # above for everything else
  } else {
    # Use specified position for all points
    pos_vector <- rep(pos, nrow(coord_matrix))
  }

  # Calculate plot limits with padding for labels
  x_range <- range(coord_matrix[,1])
  y_range <- range(coord_matrix[,2])
  x_expand <- diff(x_range) * 0.05
  y_expand <- diff(y_range) * 0.08

  # Create the plot
  plot(coord_matrix[,1], coord_matrix[,2],
       col = col_vector, pch = pch, cex = cex,
       xlim = c(x_range[1] - x_expand, x_range[2] + x_expand),
       ylim = c(y_range[1] - y_expand, y_range[2] + y_expand),
       xlab = "Dimension 1", ylab = "Dimension 2",
       main = main, ...)

  # Add labels with smart positioning
  text(coord_matrix[,1], coord_matrix[,2],
       labels = labels,
       col = col_vector,
       cex = cex.text,
       pos = pos_vector)

  # Add legend if clustering provided
  if(!is.null(color)) {
    legend("topright",
           legend = levels(color_factor),
           col = col_palette,
           pch = pch,
           title = "Cluster",
           cex = 0.8)
  }

  # Add grid for reference
  grid(col = "lightgray", lty = "dotted")

  # Return invisible list
  invisible(list(coords = coord_matrix,
                 colors = col_vector,
                 labels = labels))
}
