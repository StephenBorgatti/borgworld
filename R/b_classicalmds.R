#' Classical Multidimensional Scaling with Similarity/Dissimilarity Input
#'
#' Performs classical multidimensional scaling on a similarity or dissimilarity matrix
#' and optionally produces a plot of the results.
#'
#' @param x A square matrix of similarities or dissimilarities
#' @param type Character string specifying the input type. Must be one of
#'   "similarities", "dissimilarities", or any unambiguous abbreviation
#'   (e.g., "s", "sim", "d", "dis", "diss"). If missing, the user will be
#'   prompted to specify.
#' @param dim Number of dimensions for the MDS solution (default = 2)
#' @param add Logical, whether to add a constant to make eigenvalues non-negative (default = FALSE)
#' @param plot Logical, whether to produce a plot (default = TRUE)
#' @param labels Optional character vector of labels for points
#'
#' @return An object of class "bmds" containing:
#'   \item{points}{Matrix of MDS coordinates}
#'   \item{eig}{Eigenvalues}
#'   \item{GOF}{Goodness of fit statistics}
#'   \item{stress}{Stress measure}
#'
#' @importFrom stats cmdscale dist
#' @importFrom graphics plot text grid abline
#' @export
#'
#' @examples
#' # Using dissimilarity matrix
#' bclassicalmds(as.matrix(dist(USArrests)), "d", dim=2)
#' bclassicalmds(as.matrix(dist(USArrests)), "dissimilarities", dim=2)
#'
#' # Using similarity matrix (correlation)
#' cor_mat <- cor(mtcars)
#' bclassicalmds(cor_mat, "s", dim=2)
#' bclassicalmds(cor_mat, "similarities", dim=3, plot=FALSE)
#'
bclassicalmds <- function(x, type, dim = 2, add = FALSE,
                          plot = TRUE, labels = NULL) {

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

  # Perform classical MDS using cmdscale
  mds_result <- cmdscale(d_matrix, k = dim, eig = TRUE, add = add)

  # Extract coordinates
  coords <- mds_result$points

  # Create dimension names
  colnames(coords) <- paste0("Dim", 1:dim)
  rownames(coords) <- labels

  # Calculate goodness of fit
  eig <- mds_result$eig

  # GOF measures
  if (length(eig) > 0) {
    # Keep only positive eigenvalues for GOF calculation
    pos_eig <- eig[eig > 0]

    if (length(pos_eig) >= dim) {
      gof1 <- sum(abs(eig[1:dim])) / sum(abs(eig))
      gof2 <- sum(pos_eig[1:dim]^2) / sum(pos_eig^2)
    } else {
      gof1 <- NA
      gof2 <- NA
    }
  } else {
    gof1 <- NA
    gof2 <- NA
  }

  # Calculate stress (Kruskal's stress formula)
  # Reconstruct distances from dim-dimensional solution
  fitted_dist <- as.matrix(dist(coords))

  # Calculate stress
  stress_num <- sum((as.vector(d_matrix) - as.vector(fitted_dist))^2)
  stress_denom <- sum(as.vector(d_matrix)^2)

  if (stress_denom > 0) {
    stress <- sqrt(stress_num / stress_denom)
  } else {
    stress <- NA
  }

  # Create plot if requested
  if (plot && dim >= 2) {
    # Extract first two dimensions for plotting
    x_coord <- coords[, 1]
    y_coord <- coords[, 2]

    # Calculate axis labels with variance explained
    if (!is.na(gof2) && length(pos_eig) >= 2) {
      var_explained <- pos_eig[1:2]^2 / sum(pos_eig^2) * 100
      xlab <- sprintf("Dimension 1 (%.1f%%)", var_explained[1])
      ylab <- sprintf("Dimension 2 (%.1f%%)", var_explained[2])
    } else {
      xlab <- "Dimension 1"
      ylab <- "Dimension 2"
    }

    # Create the plot
    plot(x_coord, y_coord,
         xlab = xlab, ylab = ylab,
         main = "Classical MDS",
         pch = 16, col = "steelblue",
         xlim = range(x_coord) * 1.1,
         ylim = range(y_coord) * 1.1)

    # Add labels
    text(x_coord, y_coord,
         labels = labels,
         pos = 3, cex = 0.8)

    # Add grid
    grid(col = "lightgray", lty = "dotted")
    abline(h = 0, v = 0, col = "gray", lty = 2)
  }

  # Prepare result object
  result <- list(
    points = coords,
    eig = eig,
    GOF = c(GOF1 = gof1, GOF2 = gof2),
    stress = stress,
    dim = dim,
    n = n,
    type = type,
    call = match.call()
  )

  class(result) <- "bmds"

  return(result)
}


#' Print method for bmds objects
#'
#' @param x A bmds object
#' @param ... Additional arguments (ignored)
#'
#' @importFrom utils head
#' @export
#' @method print bmds
print.bmds <- function(x, ...) {
  cat("\nClassical Multidimensional Scaling\n")
  cat("==================================\n")
  cat("\nCall:\n")
  print(x$call)
  cat("\nNumber of objects:", x$n, "\n")
  cat("Number of dimensions:", x$dim, "\n")
  cat("Input type:", x$type, "\n")

  if (!is.na(x$stress)) {
    cat("\nStress:", sprintf("%.4f", x$stress), "\n")
  }

  if (!any(is.na(x$GOF))) {
    cat("\nGoodness of fit:\n")
    cat("  GOF1:", sprintf("%.4f", x$GOF[1]), "\n")
    cat("  GOF2:", sprintf("%.4f", x$GOF[2]), "\n")
  }

  cat("\nCoordinates:\n")
  print(head(x$points))

  if (nrow(x$points) > 6) {
    cat("...\n")
    cat("(Showing first 6 of", nrow(x$points), "points)\n")
  }

  invisible(x)
}
