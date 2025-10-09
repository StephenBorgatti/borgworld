#' Non-metric Multidimensional Scaling with Similarity/Dissimilarity Input
#'
#' Performs Kruskal's non-metric multidimensional scaling on a similarity or
#' dissimilarity matrix and optionally produces a plot of the results.
#'
#' @param x A square matrix of similarities or dissimilarities
#' @param type Character string specifying the input type. Must be one of
#'   "similarities", "dissimilarities", or any unambiguous abbreviation
#'   (e.g., "s", "sim", "d", "dis", "diss"). If missing, the user will be
#'   prompted to specify.
#' @param dim Number of dimensions for the MDS solution (default = 2)
#' @param max_iter Maximum number of iterations (default = 100)
#' @param plot Logical, whether to produce a plot (default = TRUE)
#' @param labels Optional character vector of labels for points
#' @param trace Logical/numeric for iteration trace info (0 = none, 1 = minimal,
#'   2 = full; default = FALSE which equals 0)
#' @param ... Additional arguments passed to MASS::isoMDS
#'
#' @return An object of class "bnmds" containing:
#'   \item{points}{Matrix of MDS coordinates}
#'   \item{stress}{Final stress value (Kruskal's stress formula 1)}
#'   \item{stress.type}{Type of stress measure used}
#'   \item{n_iter}{Number of iterations performed}
#'   \item{dim}{Number of dimensions}
#'   \item{n}{Number of objects}
#'   \item{type}{Input type used}
#'   \item{call}{The function call}
#'
#' @details
#' This function uses MASS::isoMDS to perform Kruskal's non-metric MDS.
#' Non-metric MDS attempts to preserve the rank order of dissimilarities
#' rather than the exact values, making it more flexible than classical MDS.
#'
#' When a similarity matrix is provided (type = "similarities"), it is converted
#' to dissimilarity using: d = max(s) - s
#'
#' The stress value reported is Kruskal's stress formula 1, which ranges from 0
#' (perfect fit) to 1 (worst fit). General guidelines for stress interpretation:
#' \itemize{
#'   \item < 0.025: Excellent
#'   \item 0.025 - 0.05: Good
#'   \item 0.05 - 0.10: Fair
#'   \item 0.10 - 0.20: Poor
#'   \item > 0.20: Unacceptable
#' }
#'
#' @importFrom MASS isoMDS
#' @importFrom stats cmdscale dist
#' @importFrom graphics plot text grid abline
#' @export
#'
#' @examples
#' # Using dissimilarity matrix
#' dist_mat <- as.matrix(dist(USArrests))
#' bnonmetricmds(dist_mat, "d", dim = 2)
#' bnonmetricmds(dist_mat, "dissimilarities", dim = 2)  # Same as above
#'
#' # Using similarity matrix (correlation)
#' cor_mat <- cor(mtcars)
#' bnonmetricmds(cor_mat, "s", dim = 2)
#' bnonmetricmds(cor_mat, "similarities", dim = 3, plot = FALSE)
#'
#' # With custom labels
#' bnonmetricmds(cor_mat, "sim", labels = colnames(cor_mat))
#'
bnonmetricmds <- function(x, type, dim = 2, max_iter = 100,
                          plot = TRUE, labels = NULL, trace = FALSE, ...) {

  # Check if type is provided, if not, prompt user
  if (missing(type)) {
    if (interactive()) {
      cat("Please specify the type of your input matrix:\n")
      cat("  1. similarities\n")
      cat("  2. dissimilarities\n")
      response <- readline("Enter your choice (1 or 2, or 's'/'d'): ")

      response <- tolower(trimws(response))
      if (response %in% c("1", "s", "sim", "similarities")) {
        type <- "similarities"
      } else if (response %in% c("2", "d", "dis", "diss", "dissimilarities")) {
        type <- "dissimilarities"
      } else {
        stop("Invalid choice. Please specify 'similarities' or 'dissimilarities'")
      }
      cat("Using type:", type, "\n")
    } else {
      stop("Argument 'type' is missing with no default. Please specify 'similarities' or 'dissimilarities'")
    }
  }

  # Match type argument with abbreviations
  type <- match.arg(type, c("similarities", "dissimilarities"))
  is_similarity <- (type == "similarities")

  # Check inputs
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  if (nrow(x) != ncol(x)) {
    stop("Input must be a square matrix")
  }

  n <- nrow(x)

  # Check for dim validity
  if (dim >= n) {
    stop("Number of dimensions 'dim' must be less than the number of objects")
  }
  if (dim < 1) {
    stop("Number of dimensions 'dim' must be at least 1")
  }

  # Store original row names for labels
  if (is.null(labels)) {
    if (!is.null(rownames(x))) {
      labels <- rownames(x)
    } else {
      labels <- as.character(1:n)
    }
  } else if (length(labels) != n) {
    stop("Length of labels must match the number of rows/columns in x")
  }

  # Convert similarity to dissimilarity if needed
  if (is_similarity) {
    # Convert similarity to dissimilarity using max(s) - s
    max_sim <- max(x, na.rm = TRUE)
    d_matrix <- max_sim - x
    diag(d_matrix) <- 0
  } else {
    d_matrix <- x
    diag(d_matrix) <- 0
  }

  # Check for non-negative values
  if (any(d_matrix < 0, na.rm = TRUE)) {
    warning("Negative values detected in dissimilarity matrix. Setting to 0.")
    d_matrix[d_matrix < 0] <- 0
  }

  # Convert matrix to dist object for isoMDS
  # We need to extract the lower triangle
  d_dist <- as.dist(d_matrix)

  # Get initial configuration using classical MDS
  # This provides a good starting point for isoMDS
  init_config <- cmdscale(d_dist, k = dim)

  # Convert trace parameter to numeric if needed
  if (is.logical(trace)) {
    trace_val <- ifelse(trace, 1, 0)
  } else {
    trace_val <- as.numeric(trace)
  }

  # Perform non-metric MDS using MASS::isoMDS
  nmds_result <- MASS::isoMDS(d_dist,
                              y = init_config,
                              k = dim,
                              maxit = max_iter,
                              trace = trace_val > 0,
                              tol = 1e-3,
                              ...)

  # Extract coordinates
  coords <- nmds_result$points

  # Create dimension names
  colnames(coords) <- paste0("Dim", 1:dim)
  rownames(coords) <- labels

  # Create plot if requested
  if (plot && dim >= 2) {
    # Extract first two dimensions for plotting
    x_coord <- coords[, 1]
    y_coord <- coords[, 2]

    # Create stress label for plot
    stress_label <- sprintf("Stress = %.3f", nmds_result$stress / 100)

    # Determine stress quality for subtitle
    stress_pct <- nmds_result$stress / 100
    if (stress_pct < 0.025) {
      quality <- "Excellent fit"
    } else if (stress_pct < 0.05) {
      quality <- "Good fit"
    } else if (stress_pct < 0.10) {
      quality <- "Fair fit"
    } else if (stress_pct < 0.20) {
      quality <- "Poor fit"
    } else {
      quality <- "Questionable fit"
    }

    # Create the plot
    plot(x_coord, y_coord,
         xlab = "Dimension 1",
         ylab = "Dimension 2",
         main = "Non-metric MDS",
         sub = paste(stress_label, "-", quality),
         pch = 16, col = "steelblue",
         xlim = range(x_coord) * 1.1,
         ylim = range(y_coord) * 1.1)

    # Add labels
    text(x_coord, y_coord,
         labels = labels,
         pos = 3, cex = 0.8)

    # Add grid and reference lines
    grid(col = "lightgray", lty = "dotted")
    abline(h = 0, v = 0, col = "gray", lty = 2)
  }

  # Prepare result object
  result <- list(
    points = coords,
    stress = nmds_result$stress / 100,  # Convert to 0-1 scale
    stress.type = "Kruskal's Stress Formula 1",
    n_iter = nmds_result$n_iter,
    dim = dim,
    n = n,
    type = type,
    call = match.call()
  )

  class(result) <- "bnmds"

  return(result)
}


#' Print method for bnmds objects
#'
#' @param x A bnmds object
#' @param ... Additional arguments (ignored)
#'
#' @importFrom utils head
#' @export
#' @method print bnmds
print.bnmds <- function(x, ...) {
  cat("\nNon-metric Multidimensional Scaling\n")
  cat("====================================\n")
  cat("\nCall:\n")
  print(x$call)
  cat("\nNumber of objects:", x$n, "\n")
  cat("Number of dimensions:", x$dim, "\n")
  cat("Input type:", x$type, "\n")

  cat("\nStress:", sprintf("%.4f", x$stress), "\n")
  cat("Stress type:", x$stress.type, "\n")

  # Interpret stress
  if (x$stress < 0.025) {
    interpretation <- "Excellent"
  } else if (x$stress < 0.05) {
    interpretation <- "Good"
  } else if (x$stress < 0.10) {
    interpretation <- "Fair"
  } else if (x$stress < 0.20) {
    interpretation <- "Poor"
  } else {
    interpretation <- "Questionable"
  }
  cat("Stress interpretation:", interpretation, "\n")

  if (!is.null(x$n_iter)) {
    cat("Iterations:", x$n_iter, "\n")
  }

  cat("\nCoordinates:\n")
  print(head(x$points))

  if (nrow(x$points) > 6) {
    cat("...\n")
    cat("(Showing first 6 of", nrow(x$points), "points)\n")
  }

  invisible(x)
}


#' Shepard plot for bnmds objects
#'
#' Creates a Shepard diagram to assess the fit of non-metric MDS
#'
#' @param x A bnmds object
#' @param original_dist Optional original dissimilarity matrix. If not provided,
#'   the function will not be able to create the plot.
#' @param ... Additional arguments passed to plot()
#'
#' @importFrom stats dist
#' @importFrom graphics plot points legend
#' @export
shepard.bnmds <- function(x, original_dist = NULL, ...) {
  if (is.null(original_dist)) {
    stop("Please provide the original dissimilarity matrix to create a Shepard plot")
  }

  # Convert to matrix if needed
  if (!is.matrix(original_dist)) {
    original_dist <- as.matrix(original_dist)
  }

  # Get distances in reduced space
  config_dist <- as.matrix(dist(x$points))

  # Extract lower triangle (excluding diagonal)
  orig_vec <- original_dist[lower.tri(original_dist)]
  config_vec <- config_dist[lower.tri(config_dist)]

  # Create Shepard plot
  plot(orig_vec, config_vec,
       xlab = "Original Dissimilarities",
       ylab = "Configuration Distances",
       main = "Shepard Diagram",
       pch = 16, col = rgb(0, 0, 1, 0.5),
       ...)

  # Add 1:1 line for reference
  abline(0, 1, col = "red", lty = 2)

  # Add smooth line if possible
  if (requireNamespace("stats", quietly = TRUE)) {
    lo <- stats::lowess(orig_vec, config_vec)
    points(lo$x, lo$y, type = "l", col = "darkgreen", lwd = 2)

    legend("bottomright",
           legend = c("Data points", "Perfect fit", "Smooth fit"),
           col = c(rgb(0, 0, 1, 0.5), "red", "darkgreen"),
           pch = c(16, NA, NA),
           lty = c(NA, 2, 1),
           lwd = c(NA, 1, 2))
  } else {
    legend("bottomright",
           legend = c("Data points", "Perfect fit"),
           col = c(rgb(0, 0, 1, 0.5), "red"),
           pch = c(16, NA),
           lty = c(NA, 2))
  }

  invisible(list(original = orig_vec, configuration = config_vec))
}
