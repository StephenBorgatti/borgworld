#' Create an Enhanced Biplot from a BPCA Object
#'
#' @description
#' Produces a biplot visualization from a \code{bpca} object with improved label
#' placement and customization options. The biplot displays both observations
#' (as points) and variables (as arrows) in the same plot using the first two
#' principal components.
#'
#' @param bpca_obj An object of class "bpca" produced by either the borgworld
#'   \code{bpca()} function or the original bpca package
#' @param scale_factor Numeric. Scaling factor for variable arrows to improve
#'   visibility. Default is 1.5
#' @param label_offset Numeric. Multiplicative offset for label positions from
#'   arrow tips to reduce overlap. Default is 1.15
#' @param point_col Color for observation points. Default is "steelblue"
#' @param arrow_col Color for variable arrows. Default is "red"
#' @param label_col Color for variable labels. Default is "darkred"
#' @param point_cex Numeric. Size of observation points. Default is 1
#' @param label_cex Numeric. Size of variable labels. Default is 0.8
#' @param main Character string. Main title for the plot. Default is "PCA Biplot"
#' @param show_grid Logical. Whether to display grid lines. Default is TRUE
#' @param ... Additional graphical parameters passed to \code{plot()}
#'
#' @return Invisibly returns a list containing:
#'   \item{scores}{Matrix of PC scores for observations (first 2 PCs)}
#'   \item{loadings}{Matrix of scaled loadings for variables}
#'   \item{var_names}{Character vector of variable names}
#'
#' @details
#' The function extracts PC scores and loadings from the bpca object and creates
#' a biplot with the following features:
#' \itemize{
#'   \item Smart label positioning based on quadrant to minimize overlap
#'   \item Adjustable arrow scaling for better visibility
#'   \item Percentage of variance explained shown on axis labels
#'   \item Optional grid lines and legend
#'   \item Customizable colors and sizes for all plot elements
#' }
#'
#' @examples
#' \dontrun{
#' # For borgworld bpca:
#' data(iris)
#' p <- bpca(iris[,1:4])
#' bpcabiplot(p)
#'
#' # Customized biplot with longer arrows and smaller labels
#' bpcabiplot(p,
#'            scale_factor = 2,
#'            label_cex = 0.7,
#'            point_col = "blue",
#'            arrow_col = "darkgreen")
#'
#' # Without grid and custom title
#' bpcabiplot(p,
#'            show_grid = FALSE,
#'            main = "Iris Dataset PCA Biplot")
#' }
#'
#' @seealso \code{\link{bpca}} for performing the PCA analysis
#'
#' @export
#' @importFrom graphics plot arrows text abline grid legend
bpcabiplot <- function(bpca_obj,
                       scale_factor = 1.5,
                       label_offset = 1.15,
                       point_col = "steelblue",
                       arrow_col = "red",
                       label_col = "darkred",
                       point_cex = 1,
                       label_cex = 0.8,
                       main = "PCA Biplot",
                       show_grid = TRUE,
                       ...) {

  # Extract components from bpca object
  # Support both borgworld bpca and original bpca package structures
  if (!is.null(bpca_obj$scores)) {
    scores <- bpca_obj$scores[, 1:2]  # borgworld bpca
    loadings <- bpca_obj$loadings[, 1:2]
    var_names <- rownames(bpca_obj$loadings)
    var_explained <- bpca_obj$var_explained[1:2]
  } else if (!is.null(bpca_obj$scor)) {
    scores <- bpca_obj$scor[, 1:2]  # original bpca package
    loadings <- bpca_obj$load[, 1:2]
    var_names <- rownames(bpca_obj$load)
    importance <- bpca_obj$import
    var_explained <- importance[2, 1:2] * 100
  } else {
    stop("Unrecognized bpca object structure")
  }

  # Scale loadings for visibility
  loadings_scaled <- loadings * scale_factor

  # Calculate plot limits to accommodate both scores and scaled loadings
  xlim <- range(c(scores[,1], loadings_scaled[,1] * 1.2))
  ylim <- range(c(scores[,2], loadings_scaled[,2] * 1.2))

  # Create the plot
  plot(scores[,1], scores[,2],
       xlim = xlim,
       ylim = ylim,
       xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
       main = main,
       pch = 19,
       col = point_col,
       cex = point_cex,
       ...)

  # Add grid if requested
  if(show_grid) {
    abline(h = 0, v = 0, col = "gray", lty = 2)
    grid(col = "lightgray", lty = "dotted")
  }

  # Add arrows for variables
  arrows(0, 0,
         loadings_scaled[,1],
         loadings_scaled[,2],
         length = 0.1,
         col = arrow_col,
         lwd = 2)

  # Calculate label positions with offset to avoid overlap
  label_x <- loadings_scaled[,1] * label_offset
  label_y <- loadings_scaled[,2] * label_offset

  # Adjust text positioning based on quadrant
  # This helps prevent labels from overlapping with arrows
  pos_adjust <- numeric(length(var_names))
  for(i in 1:length(var_names)) {
    if(label_x[i] >= 0 & label_y[i] >= 0) pos_adjust[i] <- 3  # top-right
    else if(label_x[i] < 0 & label_y[i] >= 0) pos_adjust[i] <- 2  # top-left
    else if(label_x[i] < 0 & label_y[i] < 0) pos_adjust[i] <- 1  # bottom-left
    else pos_adjust[i] <- 4  # bottom-right
  }

  # Add variable labels with improved positioning
  text(label_x, label_y,
       labels = var_names,
       col = label_col,
       cex = label_cex,
       pos = pos_adjust,
       font = 2)  # Bold font

  # Add legend
  legend("topright",
         legend = c("Observations", "Variables"),
         col = c(point_col, arrow_col),
         pch = c(19, NA),
         lty = c(NA, 1),
         lwd = c(NA, 2),
         cex = 0.8,
         bg = "white")

  # Return the scaled loadings invisibly for potential further use
  invisible(list(scores = scores,
                 loadings = loadings_scaled,
                 var_names = var_names))
}

#' Create a Biplot with Automatic Label Spreading
#'
#' @description
#' An alternative version of \code{bpcabiplot} that uses automatic label spreading
#' algorithms to prevent label overlap when many variables are present. If the
#' \code{plotrix} package is available, it will use its \code{spread.labels}
#' function for optimal label placement.
#'
#' @inheritParams bpcabiplot
#' @param spread_labels Logical. Whether to use automatic label spreading
#'   algorithms. Default is TRUE
#'
#' @return Invisibly returns a list containing:
#'   \item{scores}{Matrix of PC scores for observations (first 2 PCs)}
#'   \item{loadings}{Matrix of scaled loadings for variables}
#'   \item{var_names}{Character vector of variable names}
#'
#' @details
#' This function provides better label placement for crowded biplots by:
#' \itemize{
#'   \item Using \code{plotrix::spread.labels} if available for optimal spreading
#'   \item Falling back to a custom smart positioning algorithm if plotrix is not installed
#'   \item Detecting label clusters and automatically adjusting positions
#' }
#'
#' @examples
#' \dontrun{
#' # For datasets with many variables
#' data(USArrests)
#' p <- bpca(USArrests)
#' bpcabiplot_repel(p)
#'
#' # With custom settings
#' bpcabiplot_repel(p,
#'                  scale_factor = 1.8,
#'                  label_cex = 0.6)
#' }
#'
#' @seealso
#' \code{\link{bpcabiplot}} for the standard version,
#' \code{\link{bpca}} for performing the PCA analysis
#'
#' @export
#' @importFrom graphics plot arrows text abline grid legend
bpcabiplot_repel <- function(bpca_obj,
                             scale_factor = 1.5,
                             point_col = "steelblue",
                             arrow_col = "red",
                             label_col = "darkred",
                             point_cex = 1,
                             label_cex = 0.8,
                             main = "PCA Biplot",
                             show_grid = TRUE,
                             spread_labels = TRUE,
                             ...) {

  # Extract components
  # Support both borgworld bpca and original bpca package structures
  if (!is.null(bpca_obj$scores)) {
    scores <- bpca_obj$scores[, 1:2]  # borgworld bpca
    loadings <- bpca_obj$loadings[, 1:2]
    var_names <- rownames(bpca_obj$loadings)
    var_explained <- bpca_obj$var_explained[1:2]
  } else if (!is.null(bpca_obj$scor)) {
    scores <- bpca_obj$scor[, 1:2]  # original bpca package
    loadings <- bpca_obj$load[, 1:2]
    var_names <- rownames(bpca_obj$load)
    importance <- bpca_obj$import
    var_explained <- importance[2, 1:2] * 100
  } else {
    stop("Unrecognized bpca object structure")
  }

  # Scale loadings
  loadings_scaled <- loadings * scale_factor

  # Calculate plot limits
  xlim <- range(c(scores[,1], loadings_scaled[,1] * 1.3))
  ylim <- range(c(scores[,2], loadings_scaled[,2] * 1.3))

  # Create the plot
  plot(scores[,1], scores[,2],
       xlim = xlim,
       ylim = ylim,
       xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
       main = main,
       pch = 19,
       col = point_col,
       cex = point_cex,
       ...)

  # Add grid
  if(show_grid) {
    abline(h = 0, v = 0, col = "gray", lty = 2)
    grid(col = "lightgray", lty = "dotted")
  }

  # Add arrows
  arrows(0, 0,
         loadings_scaled[,1],
         loadings_scaled[,2],
         length = 0.1,
         col = arrow_col,
         lwd = 2)

  # Smart label placement
  if(spread_labels && requireNamespace("plotrix", quietly = TRUE)) {
    # Use plotrix::spread.labels for automatic label spreading
    plotrix::spread.labels(loadings_scaled[,1],
                           loadings_scaled[,2],
                           labels = var_names,
                           ony = FALSE,
                           offsets = 0.05,
                           col = label_col,
                           cex = label_cex,
                           font = 2)
  } else {
    # Manual smart positioning
    label_positions <- smart_label_position(loadings_scaled, var_names)
    text(label_positions$x, label_positions$y,
         labels = var_names,
         col = label_col,
         cex = label_cex,
         font = 2)
  }

  # Add legend
  legend("topright",
         legend = c("Observations", "Variables"),
         col = c(point_col, arrow_col),
         pch = c(19, NA),
         lty = c(NA, 1),
         lwd = c(NA, 2),
         cex = 0.8,
         bg = "white")

  invisible(list(scores = scores,
                 loadings = loadings_scaled,
                 var_names = var_names))
}

#' Smart Label Positioning for Biplots
#'
#' @description
#' Internal helper function that calculates optimal label positions to minimize
#' overlap in biplots. It detects clusters of nearby labels and adjusts their
#' positions accordingly.
#'
#' @param coords Matrix with 2 columns containing x and y coordinates of points
#' @param labels Character vector of labels (used only for length)
#' @param offset Numeric. Base offset multiplier for label distance from points.
#'   Default is 1.15
#'
#' @return A list with components:
#'   \item{x}{Numeric vector of adjusted x coordinates for labels}
#'   \item{y}{Numeric vector of adjusted y coordinates for labels}
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Calculates angles for each point from the origin
#'   \item Applies base offset to all points
#'   \item Detects pairs of labels that would overlap (distance < 0.15)
#'   \item Adjusts overlapping labels by increasing their offset
#' }
#'
#' @keywords internal
#' @noRd
smart_label_position <- function(coords, labels, offset = 1.15) {
  n <- nrow(coords)
  new_x <- numeric(n)
  new_y <- numeric(n)

  # Calculate angles for each point
  angles <- atan2(coords[,2], coords[,1])

  # Sort by angle to detect clusters
  angle_order <- order(angles)

  for(i in 1:n) {
    # Base position with offset
    new_x[i] <- coords[i,1] * offset
    new_y[i] <- coords[i,2] * offset

    # Check for nearby points and adjust
    if(i > 1) {
      for(j in 1:(i-1)) {
        dist <- sqrt((new_x[i] - new_x[j])^2 + (new_y[i] - new_y[j])^2)
        if(dist < 0.15) {  # Too close, need adjustment
          # Push away from neighbor
          angle <- atan2(coords[i,2], coords[i,1])
          new_x[i] <- coords[i,1] * (offset + 0.1)
          new_y[i] <- coords[i,2] * (offset + 0.1)
        }
      }
    }
  }

  return(list(x = new_x, y = new_y))
}
