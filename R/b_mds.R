#' Classical Multidimensional Scaling with Similarity/Dissimilarity Input
#'
#' Performs classical multidimensional scaling on a similarity or dissimilarity matrix
#' and optionally produces a plot of the results.
#'
#' @param x A square matrix of similarities or dissimilarities
#' @param type Character string, either "s" for similarity or "d" for dissimilarity
#' @param k Number of dimensions for the MDS solution (default = 2)
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
#' @importFrom graphics plot text
#' @export
#'
#' @examples
#' # Using dissimilarity matrix
#' bclassicalmds(as.matrix(dist(USArrests)), type="d", k=2)
#'
#' # Using similarity matrix (convert correlation to similarity)
#' cor_mat <- cor(mtcars)
#' bclassicalmds(cor_mat, type="s", k=2)
#'
bclassicalmds <- function(x, type = c("s", "d"), k = 2, add = FALSE,
                          plot = TRUE, labels = NULL) {

  # Check inputs
  type <- match.arg(type)

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  if (nrow(x) != ncol(x)) {
    stop("Input must be a square matrix")
  }

  n <- nrow(x)

  # Convert similarity to dissimilarity if needed
  if (type == "s") {
    # Convert similarity to dissimilarity
    # Common approach: d = sqrt(2(1-s)) for correlations
    # or d = max(s) - s for general similarities
    max_sim <- max(x, na.rm = TRUE)
    d_matrix <- max_sim - x
    diag(d_matrix) <- 0
  } else {
    d_matrix <- x
    diag(d_matrix) <- 0
  }

  # Perform classical MDS using cmdscale
  mds_result <- cmdscale(d_matrix, k = k, eig = TRUE, add = add)

  # Extract coordinates
  coords <- mds_result$points

  # Create dimension names
  colnames(coords) <- paste0("Dim", 1:k)

  # Add row labels
  if (!is.null(labels)) {
    rownames(coords) <- labels
  } else if (!is.null(rownames(x))) {
    rownames(coords) <- rownames(x)
  } else {
    rownames(coords) <- 1:n
  }

  # Calculate goodness of fit
  eig <- mds_result$eig

  # GOF measures
  if (length(eig) > 0) {
    # Keep only positive eigenvalues for GOF calculation
    pos_eig <- eig[eig > 0]

    if (length(pos_eig) >= k) {
      gof1 <- sum(abs(eig[1:k])) / sum(abs(eig))
      gof2 <- sum(pos_eig[1:k]^2) / sum(pos_eig^2)
    } else {
      gof1 <- NA
      gof2 <- NA
    }
  } else {
    gof1 <- NA
    gof2 <- NA
  }

  # Calculate stress (Kruskal's stress formula)
  # Reconstruct distances from k-dimensional solution
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
  if (plot && k >= 2) {
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
    if (!is.null(rownames(coords))) {
      text(x_coord, y_coord,
           labels = rownames(coords),
           pos = 3, cex = 0.8)
    }

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
    k = k,
    n = n,
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
#' @export
#' @method print bmds
print.bmds <- function(x, ...) {
  cat("\nClassical Multidimensional Scaling\n")
  cat("==================================\n")
  cat("\nCall:\n")
  print(x$call)
  cat("\nNumber of objects:", x$n, "\n")
  cat("Number of dimensions:", x$k, "\n")

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
