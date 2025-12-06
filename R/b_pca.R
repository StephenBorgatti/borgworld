#' Perform PCA with Social Science Conventions
#'
#' @description
#' Performs Principal Component Analysis and returns results in the format
#' expected by social scientists familiar with Stata/SPSS conventions.
#' Unlike R's built-in functions, this returns loadings as correlations
#' between variables and components, and can optionally standardize scores.
#'
#' @param data A numeric matrix or data frame containing the variables to analyze.
#'   Non-numeric columns are automatically removed.
#' @param max Integer; maximum number of components to retain.
#'   Default is NULL (returns all components)
#' @param mineigen Numeric; minimum eigenvalue for component retention.
#'   Default is 1 (Kaiser criterion). Components with eigenvalues below this
#'   threshold are excluded from loadings and scores but shown in eigenvalue table.
#' @param standardize_scores Logical; if TRUE (default), standardizes component
#'   scores to have variance = 1 (matching Stata's default behavior).
#'   If FALSE, returns raw scores with variance equal to eigenvalues.
#' @param na_action How to handle missing values: "complete" (default) for listwise
#'   deletion, "pairwise" for pairwise complete observations in correlation.
#' @param verbose Logical; if TRUE, print messages about data cleaning.
#'
#' @return An object of class \code{bpca} containing:
#' \item{eigenvalues}{Numeric vector of eigenvalues of the correlation matrix}
#' \item{variance_explained}{Data frame with variance explained by each component}
#' \item{loadings}{Matrix of component loadings (correlations between variables and PCs)}
#' \item{eigenvectors}{Matrix of raw eigenvectors (component coefficients)}
#' \item{scores}{Matrix of component scores (standardized or raw based on parameter)}
#' \item{scores_raw}{Matrix of raw component scores (always unstandardized)}
#' \item{communalities}{Numeric vector of communalities for each variable}
#' \item{correlation_matrix}{The correlation matrix used for PCA}
#' \item{standardized_data}{The standardized data matrix}
#' \item{retained}{Number of components retained based on mineigen criterion}
#'
#' @details
#' This function bridges the gap between R's PCA implementations and the
#' conventions used in social science software like Stata and SPSS:
#'
#' \itemize{
#'   \item Loadings are computed as eigenvector * sqrt(eigenvalue), giving
#'     correlations between variables and components
#'   \item Scores are standardized by default to have unit variance
#'   \item Output includes communalities and Kaiser criterion
#'   \item Terminology matches social science conventions
#' }
#'
#' The function uses the correlation matrix (variables standardized with n-1)
#' which is standard in both R and Stata/SPSS for correlation-based PCA.
#'
#' @examples
#' # Create example data with correlation structure
#' set.seed(123)
#' n <- 100
#' data <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
#' colnames(data) <- paste0("Var", 1:5)
#' data[, 2] <- data[, 1] * 0.7 + rnorm(n) * 0.5
#' data[, 3] <- data[, 1] * 0.4 + data[, 2] * 0.4 + rnorm(n) * 0.5
#'
#' # Run PCA with Kaiser criterion (default)
#' result <- bpca(data)
#'
#' # Run PCA with custom eigenvalue threshold
#' result <- bpca(data, mineigen = 0.7)
#'
#' # Run PCA with maximum 3 components
#' result <- bpca(data, max = 3)
#'
#' # Access loadings (correlations)
#' result$loadings
#'
#' # Run PCA with raw scores
#' result_raw <- bpca(data, max = 3, standardize_scores = FALSE)
#'
#' # With data frame input (non-numeric columns automatically removed)
#' result <- bpca(iris)  # Species column removed automatically
#'
#' @seealso
#' \code{\link{prcomp}} for R's SVD-based PCA implementation,
#' \code{\link{princomp}} for R's eigen-based PCA implementation,
#' \code{\link{bcompare_pca}} for converting between R and social science conventions
#'
#' @references
#' Jolliffe, I. T. (2002). Principal Component Analysis (2nd ed.). Springer.
#'
#' Hair, J. F., Black, W. C., Babin, B. J., & Anderson, R. E. (2019).
#' Multivariate Data Analysis (8th ed.). Cengage Learning.
#'
#' @importFrom stats cor sd
#' @importFrom utils capture.output
#' @export
bpca <- function(data, max = NULL, mineigen = 1, standardize_scores = TRUE,
                 na_action = "complete", verbose = FALSE) {


  # Use standardized input handling
  data <- bprepare_data(data,
                        na_action = na_action,
                        numeric_only = TRUE,
                        output_format = "data.frame",
                        verbose = verbose)

  # Convert to matrix for computations
  X <- as.matrix(data)
  n <- nrow(X)
  p <- ncol(X)

  # Get variable names
  var_names <- bget_labels(X, dim = "cols", default_prefix = "Var")
  colnames(X) <- var_names

  # Default to all components for initial calculation
  if (is.null(max)) {
    max_components <- min(n - 1, p)
  } else {
    max_components <- min(max, n - 1, p)
  }

  # Standardize data using sample standard deviation (n-1)
  X_std <- scale(X, center = TRUE, scale = TRUE)

  # Using eigen decomposition of correlation matrix
  R <- cor(X)
  eigen_decomp <- eigen(R)

  # Get all eigenvalues for display
  all_eigenvalues <- eigen_decomp$values
  names(all_eigenvalues) <- paste0("PC", 1:length(all_eigenvalues))

  # Determine components to retain based on mineigen criterion
  retained_indices <- which(eigen_decomp$values >= mineigen)
  if (length(retained_indices) > 0) {
    n_retained <- min(max(retained_indices), max_components)
    retained_indices <- 1:n_retained
  } else {
    # If no components meet criterion, retain at least 1
    n_retained <- 1
    retained_indices <- 1
  }

  # Extract eigenvalues and eigenvectors for retained components
  eigenvalues <- eigen_decomp$values[retained_indices]
  names(eigenvalues) <- paste0("PC", retained_indices)

  eigenvectors <- eigen_decomp$vectors[, retained_indices, drop = FALSE]
  rownames(eigenvectors) <- var_names
  colnames(eigenvectors) <- paste0("PC", retained_indices)

  # Compute loadings as correlations (what social scientists expect)
  loadings <- matrix(NA, nrow = p, ncol = n_retained)
  for (i in 1:n_retained) {
    loadings[, i] <- eigenvectors[, i] * sqrt(eigenvalues[i])
  }
  rownames(loadings) <- var_names
  colnames(loadings) <- paste0("PC", retained_indices)

  # Compute scores
  scores_raw <- X_std %*% eigenvectors
  colnames(scores_raw) <- paste0("PC", retained_indices)
  if (!is.null(rownames(X))) {
    rownames(scores_raw) <- rownames(X)
  }

  # Standardize scores if requested
  if (standardize_scores) {
    scores <- scale(scores_raw, center = FALSE, scale = TRUE)
    colnames(scores) <- paste0("PC", retained_indices)
    if (!is.null(rownames(X))) {
      rownames(scores) <- rownames(X)
    }
  } else {
    scores <- scores_raw
  }

  # Compute variance explained (for all components)
  total_variance <- sum(all_eigenvalues)
  all_prop_variance <- all_eigenvalues / total_variance
  all_cum_variance <- cumsum(all_prop_variance)

  # Compute communalities (only for retained components)
  communalities <- rowSums(loadings^2)
  names(communalities) <- var_names

  # --- Print formatted output using shared utilities ---

  bprint_header("PCA RESULTS", subtitle = "Stata/SPSS Style Output")

  # Print eigenvalue table using shared utility
  bprint_eigenvalues(all_eigenvalues, threshold = mineigen)

  # Additional info
  n_kaiser <- sum(all_eigenvalues > 1)
  cat("\nMinimum eigenvalue threshold (mineigen):", mineigen, "\n")
  cat("Components retained:", n_retained, "\n")

  # Print loadings
  bprint_section("Component Loadings (Variable-Component Correlations)")
  cat("Showing loadings for", n_retained, "retained components\n\n")
  bprint_loadings(loadings, max_vars = 15)

  # Print communalities
  bprint_section("Communalities (Variance Explained per Variable)")
  cat("Based on", n_retained, "retained components\n\n")

  n_show <- min(15, p)
  comm_df <- data.frame(
    Variable = var_names[1:n_show],
    Communality = communalities[1:n_show]
  )
  bprint_table(comm_df)
  if (p > n_show) {
    cat("... [", p - n_show, " more variables]\n", sep = "")
  }

  cat("\n")

  # Return results with all names properly set
  result <- list(
    eigenvalues = eigenvalues,
    all_eigenvalues = all_eigenvalues,
    variance_explained = data.frame(
      component = 1:length(all_eigenvalues),
      eigenvalue = all_eigenvalues,
      prop_variance = all_prop_variance,
      cum_variance = all_cum_variance
    ),
    loadings = loadings,
    eigenvectors = eigenvectors,
    scores = scores,
    scores_raw = scores_raw,
    communalities = communalities,
    correlation_matrix = R,
    standardized_data = X_std,
    retained = n_retained,
    mineigen = mineigen,
    n = n,
    p = p
  )

  class(result) <- "bpca"
  invisible(result)
}


#' Print method for bpca objects
#'
#' @param x A bpca object
#' @param ... Additional arguments (unused)
#' @export
print.bpca <- function(x, ...) {
  cat("\nPCA Results (Stata/SPSS Style)\n")
  cat("==============================\n")
  cat("Observations:", x$n, "\n")
  cat("Variables:", x$p, "\n")
  cat("Components retained:", x$retained, "\n")
  cat("Mineigen threshold:", x$mineigen, "\n\n")

  cat("Eigenvalues (retained components):\n")
  print(round(x$eigenvalues, 4))
  cat("\n")

  cat("Component Loadings (first few variables):\n")
  print(round(head(x$loadings, 10), 3))

  invisible(x)
}
