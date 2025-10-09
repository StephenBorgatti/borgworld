#' Compare Multiple MDS Methods
#'
#' @description
#' Runs both classical and non-metric MDS on the same data for comparison.
#'
#' @param x Input matrix (similarity or dissimilarity)
#' @param type Character string specifying the input type. Must be one of
#'   "similarities", "dissimilarities", or any unambiguous abbreviation.
#'   If missing, the user will be prompted to specify.
#' @param dim Number of dimensions (default = 2)
#' @param plot Logical, whether to create comparison plots (default = TRUE)
#' @param ... Additional arguments passed to MDS functions
#'
#' @return A list of class "mds_comparison" containing both MDS results and comparison statistics
#'
#' @export
#' @examples
#' # Compare methods on dissimilarity matrix
#' dist_mat <- as.matrix(dist(USArrests))
#' comparison <- compare_mds(dist_mat, "d")
#'
#' # Compare methods on correlation matrix
#' cor_mat <- cor(mtcars)
#' comparison <- compare_mds(cor_mat, "similarities", dim = 3)
#'
compare_mds <- function(x, type, dim = 2, plot = TRUE, ...) {

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

  # Run both methods
  cat("Running Classical MDS...\n")
  classical_result <- bclassicalmds(x, type = type, dim = dim, plot = FALSE, ...)

  cat("Running Non-metric MDS...\n")
  nonmetric_result <- bnonmetricmds(x, type = type, dim = dim, plot = FALSE, ...)

  # Calculate Procrustes statistics for comparison
  proc_stats <- procrustes_compare(classical_result$points,
                                   nonmetric_result$points)

  # Create comparison plots if requested
  if (plot && dim >= 2) {
    oldpar <- par(mfrow = c(1, 2))
    on.exit(par(oldpar))

    # Plot classical MDS
    plot(classical_result$points[,1:2],
         main = paste("Classical MDS\nStress:",
                      ifelse(!is.na(classical_result$stress),
                             sprintf("%.3f", classical_result$stress), "N/A")),
         xlab = "Dimension 1", ylab = "Dimension 2",
         pch = 16, col = "steelblue")
    text(classical_result$points[,1:2],
         labels = rownames(classical_result$points),
         pos = 3, cex = 0.7)
    grid(col = "lightgray", lty = "dotted")
    abline(h = 0, v = 0, col = "gray", lty = 2)

    # Plot non-metric MDS
    plot(nonmetric_result$points[,1:2],
         main = paste("Non-metric MDS\nStress:",
                      sprintf("%.3f", nonmetric_result$stress)),
         xlab = "Dimension 1", ylab = "Dimension 2",
         pch = 16, col = "darkgreen")
    text(nonmetric_result$points[,1:2],
         labels = rownames(nonmetric_result$points),
         pos = 3, cex = 0.7)
    grid(col = "lightgray", lty = "dotted")
    abline(h = 0, v = 0, col = "gray", lty = 2)
  }

  # Create comparison object
  comparison <- list(
    classical = classical_result,
    nonmetric = nonmetric_result,
    procrustes = proc_stats,
    input_type = type
  )

  class(comparison) <- "mds_comparison"

  return(comparison)
}


#' Simple Procrustes Analysis
#'
#' @description
#' Compares two MDS configurations using Procrustes analysis.
#' This function aligns two configurations optimally and measures their similarity.
#'
#' @param X First configuration matrix
#' @param Y Second configuration matrix
#'
#' @return List containing Procrustes statistics:
#'   \item{distance}{Procrustes distance between configurations}
#'   \item{normalized_distance}{Normalized Procrustes distance (0-1 scale)}
#'   \item{correlation}{Correlation between aligned configurations}
#'   \item{rotation_matrix}{Optimal rotation matrix}
#'
#' @noRd
procrustes_compare <- function(X, Y) {
  # Center both configurations
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  Y_centered <- scale(Y, center = TRUE, scale = FALSE)

  # Calculate sum of squares
  ss_X <- sum(X_centered^2)
  ss_Y <- sum(Y_centered^2)

  # Procrustes rotation using SVD
  svd_result <- svd(t(Y_centered) %*% X_centered)
  R <- svd_result$v %*% t(svd_result$u)

  # Apply rotation
  Y_rotated <- Y_centered %*% R

  # Calculate Procrustes distance
  proc_dist <- sqrt(sum((X_centered - Y_rotated)^2))

  # Normalized Procrustes distance
  proc_dist_norm <- proc_dist / sqrt(ss_X)

  # Correlation between configurations
  config_cor <- cor(as.vector(X_centered), as.vector(Y_rotated))

  return(list(
    distance = proc_dist,
    normalized_distance = proc_dist_norm,
    correlation = config_cor,
    rotation_matrix = R
  ))
}


#' Print Method for MDS Comparison
#'
#' @param x An mds_comparison object
#' @param ... Additional arguments (ignored)
#'
#' @export
#' @method print mds_comparison
print.mds_comparison <- function(x, ...) {
  cat("\nMDS Methods Comparison\n")
  cat("======================\n")
  cat("Input type:", x$input_type, "\n\n")

  cat("Classical MDS:\n")
  cat("  Stress:", ifelse(!is.na(x$classical$stress),
                          sprintf("%.4f", x$classical$stress), "N/A"), "\n")
  if (!any(is.na(x$classical$GOF))) {
    cat("  GOF1:", sprintf("%.4f", x$classical$GOF[1]), "\n")
    cat("  GOF2:", sprintf("%.4f", x$classical$GOF[2]), "\n")
  }
  cat("  Dimensions:", x$classical$dim, "\n")
  cat("  Objects:", x$classical$n, "\n")

  cat("\nNon-metric MDS:\n")
  cat("  Stress:", sprintf("%.4f", x$nonmetric$stress), "\n")
  cat("  Iterations:", x$nonmetric$n_iter, "\n")
  cat("  Dimensions:", x$nonmetric$dim, "\n")
  cat("  Objects:", x$nonmetric$n, "\n")

  cat("\nProcrustes Comparison:\n")
  cat("  Distance:", sprintf("%.4f", x$procrustes$distance), "\n")
  cat("  Normalized distance:", sprintf("%.4f", x$procrustes$normalized_distance), "\n")
  cat("  Configuration correlation:", sprintf("%.4f", x$procrustes$correlation), "\n")

  # Interpretation of correlation
  if (x$procrustes$correlation > 0.95) {
    interpretation <- "Excellent agreement"
  } else if (x$procrustes$correlation > 0.90) {
    interpretation <- "Very good agreement"
  } else if (x$procrustes$correlation > 0.80) {
    interpretation <- "Good agreement"
  } else if (x$procrustes$correlation > 0.70) {
    interpretation <- "Moderate agreement"
  } else {
    interpretation <- "Poor agreement"
  }
  cat("  Interpretation:", interpretation, "\n")

  invisible(x)
}


#' Summary Method for MDS Comparison
#'
#' @param object An mds_comparison object
#' @param ... Additional arguments (ignored)
#'
#' @export
#' @method summary mds_comparison
summary.mds_comparison <- function(object, ...) {
  print.mds_comparison(object, ...)

  # Additional summary statistics
  cat("\nConfiguration Details:\n")
  cat("  Number of dimensions compared:", ncol(object$classical$points), "\n")
  cat("  Number of objects:", nrow(object$classical$points), "\n")

  # Stress comparison
  if (!is.na(object$classical$stress)) {
    stress_diff <- abs(object$classical$stress - object$nonmetric$stress)
    cat("\nStress difference:", sprintf("%.4f", stress_diff), "\n")

    if (stress_diff < 0.01) {
      cat("  Both methods achieve very similar stress\n")
    } else if (object$classical$stress < object$nonmetric$stress) {
      cat("  Classical MDS achieves lower stress\n")
    } else {
      cat("  Non-metric MDS achieves lower stress\n")
    }
  }

  invisible(object)
}


#' Plot Method for MDS Comparison
#'
#' @description
#' Creates side-by-side plots or overlay plot for MDS comparison
#'
#' @param x An mds_comparison object
#' @param type Type of plot: "separate" for side-by-side or "overlay" for superimposed
#' @param ... Additional arguments passed to plot functions
#'
#' @export
#' @method plot mds_comparison
plot.mds_comparison <- function(x, type = "separate", ...) {

  if (type == "separate") {
    # Side-by-side plots (default behavior)
    oldpar <- par(mfrow = c(1, 2))
    on.exit(par(oldpar))

    # Plot classical MDS
    plot(x$classical$points[,1:2],
         main = paste("Classical MDS\nStress:",
                      ifelse(!is.na(x$classical$stress),
                             sprintf("%.3f", x$classical$stress), "N/A")),
         xlab = "Dimension 1", ylab = "Dimension 2",
         pch = 16, col = "steelblue", ...)
    text(x$classical$points[,1:2],
         labels = rownames(x$classical$points),
         pos = 3, cex = 0.7)
    grid(col = "lightgray", lty = "dotted")
    abline(h = 0, v = 0, col = "gray", lty = 2)

    # Plot non-metric MDS
    plot(x$nonmetric$points[,1:2],
         main = paste("Non-metric MDS\nStress:",
                      sprintf("%.3f", x$nonmetric$stress)),
         xlab = "Dimension 1", ylab = "Dimension 2",
         pch = 16, col = "darkgreen", ...)
    text(x$nonmetric$points[,1:2],
         labels = rownames(x$nonmetric$points),
         pos = 3, cex = 0.7)
    grid(col = "lightgray", lty = "dotted")
    abline(h = 0, v = 0, col = "gray", lty = 2)

  } else if (type == "overlay") {
    # Overlay plot - align using Procrustes first
    X_centered <- scale(x$classical$points[,1:2], center = TRUE, scale = FALSE)
    Y_centered <- scale(x$nonmetric$points[,1:2], center = TRUE, scale = FALSE)

    # Apply Procrustes rotation to align
    Y_aligned <- Y_centered %*% x$procrustes$rotation_matrix

    # Determine plot limits
    xlim <- range(c(X_centered[,1], Y_aligned[,1])) * 1.1
    ylim <- range(c(X_centered[,2], Y_aligned[,2])) * 1.1

    # Create overlay plot
    plot(X_centered[,1:2],
         main = "MDS Comparison (Procrustes-aligned)",
         xlab = "Dimension 1", ylab = "Dimension 2",
         pch = 16, col = "steelblue",
         xlim = xlim, ylim = ylim, ...)
    points(Y_aligned, pch = 17, col = "darkgreen")

    # Add labels
    text(X_centered[,1:2],
         labels = rownames(x$classical$points),
         pos = 3, cex = 0.7, col = "steelblue")

    # Connect corresponding points
    for (i in 1:nrow(X_centered)) {
      lines(c(X_centered[i,1], Y_aligned[i,1]),
            c(X_centered[i,2], Y_aligned[i,2]),
            col = "gray", lty = 2)
    }

    # Add legend
    legend("topright",
           legend = c("Classical MDS", "Non-metric MDS"),
           pch = c(16, 17),
           col = c("steelblue", "darkgreen"))

    grid(col = "lightgray", lty = "dotted")
    abline(h = 0, v = 0, col = "gray", lty = 2)
  } else {
    stop("Plot type must be 'separate' or 'overlay'")
  }

  invisible(x)
}
