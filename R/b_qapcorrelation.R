#' QAP Correlation with Permutation Test
#'
#' Computes the correlation between two matrices using the Quadratic Assignment
#' Procedure (QAP) with permutation testing. If both matrices are symmetric,
#' only the lower triangle values are used for the correlation.
#'
#' @param y First matrix or data frame to correlate
#' @param x Second matrix or data frame to correlate
#' @param nperm Number of permutations for the permutation test (default: 10000)
#' @param method Correlation method: "pearson" (default), "spearman", or "kendall"
#'
#' @return A list of class "bqapcor" containing:
#' \describe{
#'   \item{observed.cor}{The observed correlation coefficient}
#'   \item{classical.p}{Classical parametric p-value}
#'   \item{permutation.p}{Permutation-based p-value (two-tailed)}
#'   \item{perm.mean}{Mean of the permutation distribution}
#'   \item{perm.sd}{Standard deviation of the permutation distribution}
#'   \item{perm.min}{Minimum value in the permutation distribution}
#'   \item{perm.max}{Maximum value in the permutation distribution}
#'   \item{n.obs}{Number of non-missing observations used}
#'   \item{nperm}{Number of permutations performed}
#'   \item{method}{Correlation method used}
#'   \item{symmetric}{Logical indicating if both matrices were symmetric}
#'   \item{var.names}{Names of the variables (y and x)}
#' }
#'
#' @details
#' The QAP correlation is used to test the correlation between two matrices
#' while accounting for the non-independence of matrix elements. The permutation
#' test randomly permutes the rows and columns of one matrix simultaneously,
#' preserving the structure of dependencies.
#'
#' For symmetric matrices (e.g., network adjacency matrices), only the lower
#' triangle values are used to avoid double-counting relationships.
#'
#' Missing values are handled automatically using pairwise deletion.
#'
#' @export
#' @examples
#' # Generate two random symmetric matrices
#' set.seed(123)
#' n <- 10
#' y <- matrix(rnorm(n*n), n, n)
#' y <- (y + t(y)) / 2  # Make symmetric
#' diag(y) <- 0
#'
#' x <- matrix(rnorm(n*n), n, n)
#' x <- (x + t(x)) / 2  # Make symmetric
#' diag(x) <- 0
#'
#' # Run QAP correlation with 1000 permutations
#' result <- bqapcorrelation(y, x, nperm = 1000)
#' print(result)
bqapcorrelation <- function(y, x, nperm = 10000, method = "pearson") {

  # Get variable names
  y_name <- deparse(substitute(y))
  x_name <- deparse(substitute(x))

  # Convert to matrix if needed
  if (is.data.frame(y)) {
    y <- as.matrix(y)
  }
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  # Validate inputs
  if (!is.matrix(y) || !is.matrix(x)) {
    stop("Both y and x must be matrices or data frames", call. = FALSE)
  }

  if (!all(dim(y) == dim(x))) {
    stop("Matrices y and x must have the same dimensions", call. = FALSE)
  }

  if (nperm < 1) {
    stop("nperm must be at least 1", call. = FALSE)
  }

  # Check if matrices are symmetric
  is_symmetric_y <- isSymmetric(y, tol = sqrt(.Machine$double.eps))
  is_symmetric_x <- isSymmetric(x, tol = sqrt(.Machine$double.eps))
  both_symmetric <- is_symmetric_y && is_symmetric_x

  # Extract values for correlation
  if (both_symmetric) {
    # Use lower triangle only (excluding diagonal)
    lower_idx <- lower.tri(y, diag = FALSE)
    y_vals <- y[lower_idx]
    x_vals <- x[lower_idx]
  } else {
    # Use all values
    y_vals <- as.vector(y)
    x_vals <- as.vector(x)
  }

  # Remove missing values pairwise
  complete_cases <- complete.cases(y_vals, x_vals)
  y_vals <- y_vals[complete_cases]
  x_vals <- x_vals[complete_cases]
  n_obs <- length(y_vals)

  if (n_obs < 3) {
    stop("Insufficient non-missing observations for correlation", call. = FALSE)
  }

  # Compute observed correlation
  obs_cor <- cor(y_vals, x_vals, method = method)

  # Compute classical p-value (for Pearson only)
  if (method == "pearson" && n_obs >= 3) {
    # Test statistic: t = r * sqrt(n-2) / sqrt(1-r^2)
    t_stat <- obs_cor * sqrt(n_obs - 2) / sqrt(1 - obs_cor^2)
    classical_p <- 2 * pt(abs(t_stat), df = n_obs - 2, lower.tail = FALSE)
  } else {
    classical_p <- NA_real_
  }

  # Permutation test
  perm_cors <- numeric(nperm)
  n <- nrow(y)

  for (i in 1:nperm) {
    # Generate random permutation
    perm <- sample(n)

    # Permute rows and columns of x simultaneously
    x_perm <- x[perm, perm, drop = FALSE]

    # Extract values
    if (both_symmetric) {
      x_perm_vals <- x_perm[lower_idx]
    } else {
      x_perm_vals <- as.vector(x_perm)
    }

    # Remove missing values
    x_perm_vals <- x_perm_vals[complete_cases]

    # Compute permuted correlation
    perm_cors[i] <- cor(y_vals, x_perm_vals, method = method)
  }

  # Compute permutation statistics
  perm_mean <- mean(perm_cors)
  perm_sd <- sd(perm_cors)
  perm_min <- min(perm_cors)
  perm_max <- max(perm_cors)

  # Compute two-tailed permutation p-value
  # Proportion of permuted correlations as extreme or more extreme than observed
  perm_p <- mean(abs(perm_cors) >= abs(obs_cor))

  # Create result object
  result <- structure(
    list(
      observed.cor = obs_cor,
      classical.p = classical_p,
      permutation.p = perm_p,
      perm.mean = perm_mean,
      perm.sd = perm_sd,
      perm.min = perm_min,
      perm.max = perm_max,
      n.obs = n_obs,
      nperm = nperm,
      method = method,
      symmetric = both_symmetric,
      var.names = c(y = y_name, x = x_name)
    ),
    class = "bqapcor"
  )

  return(result)
}


#' Print method for bqapcor objects
#'
#' @param x An object of class "bqapcor"
#' @param ... Additional arguments (ignored)
#' @export
print.bqapcor <- function(x, ...) {

  cat("\n")
  cat("QAP Correlation Analysis\n")
  cat(paste0(rep("=", 70), collapse = ""), "\n\n")

  # Variable information
  cat(sprintf("Variables:       %s vs %s\n", x$var.names["y"], x$var.names["x"]))
  cat(sprintf("Method:          %s correlation\n",
              tools::toTitleCase(x$method)))
  cat(sprintf("Observations:    %d (non-missing)\n", x$n.obs))
  cat(sprintf("Permutations:    %d\n", x$nperm))
  cat(sprintf("Matrix type:     %s\n\n",
              if(x$symmetric) "Symmetric (lower triangle used)" else "Non-symmetric"))

  # Results table
  cat("Results:\n")
  cat(paste0(rep("-", 70), collapse = ""), "\n")
  cat(sprintf("%-30s %12s\n", "Statistic", "Value"))
  cat(paste0(rep("-", 70), collapse = ""), "\n")

  cat(sprintf("%-30s %12.6f\n", "Observed correlation", x$observed.cor))

  if (!is.na(x$classical.p)) {
    if (x$classical.p < 0.0001) {
      cat(sprintf("%-30s %12s\n", "Classical p-value", "< 0.0001"))
    } else {
      cat(sprintf("%-30s %12.6f\n", "Classical p-value", x$classical.p))
    }
  }

  if (x$permutation.p < 0.0001) {
    cat(sprintf("%-30s %12s\n", "Permutation p-value", "< 0.0001"))
  } else {
    cat(sprintf("%-30s %12.6f\n", "Permutation p-value", x$permutation.p))
  }

  cat("\n")
  cat("Permutation Distribution:\n")
  cat(paste0(rep("-", 70), collapse = ""), "\n")
  cat(sprintf("%-30s %12.6f\n", "Mean", x$perm.mean))
  cat(sprintf("%-30s %12.6f\n", "Std. Deviation", x$perm.sd))
  cat(sprintf("%-30s %12.6f\n", "Minimum", x$perm.min))
  cat(sprintf("%-30s %12.6f\n", "Maximum", x$perm.max))
  cat(paste0(rep("=", 70), collapse = ""), "\n\n")

  invisible(x)
}
