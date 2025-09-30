#' Perform PCA with Social Science Conventions
#'
#' @description
#' Performs Principal Component Analysis and returns results in the format
#' expected by social scientists familiar with Stata/SPSS conventions.
#' Unlike R's built-in functions, this returns loadings as correlations
#' between variables and components, and can optionally standardize scores.
#'
#' @param data A numeric matrix or data frame containing the variables to analyze
#' @param n_components Integer; number of components to retain.
#'   Default is NULL (returns all components)
#' @param standardize_scores Logical; if TRUE (default), standardizes component
#'   scores to have variance = 1 (matching Stata's default behavior).
#'   If FALSE, returns raw scores with variance equal to eigenvalues.
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
#' # Run PCA with standardized scores (default)
#' result <- bpca(data, n_components = 3)
#'
#' # Access loadings (correlations)
#' result$loadings
#'
#' # Run PCA with raw scores
#' result_raw <- bpca(data, n_components = 3, standardize_scores = FALSE)
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
bpca <- function(data, n_components = NULL, standardize_scores = TRUE) {
  # Convert to clean matrix
  data <- bclean(data)
  X <- as.matrix(data)
  n <- nrow(X)
  p <- ncol(X)

  # Get variable names
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("Var", 1:p)
  }
  var_names <- colnames(X)

  # Default to all components
  if (is.null(n_components)) {
    n_components <- min(n - 1, p)
  }

  # Standardize data using sample standard deviation (n-1)
  X_std <- scale(X, center = TRUE, scale = TRUE)

  # Using eigen decomposition of correlation matrix
  R <- cor(X)
  eigen_decomp <- eigen(R)

  # Extract eigenvalues and eigenvectors
  eigenvalues <- eigen_decomp$values[1:n_components]
  names(eigenvalues) <- paste0("PC", 1:n_components)

  eigenvectors <- eigen_decomp$vectors[, 1:n_components, drop = FALSE]
  rownames(eigenvectors) <- var_names
  colnames(eigenvectors) <- paste0("PC", 1:n_components)

  # Compute loadings as correlations (what social scientists expect)
  loadings <- matrix(NA, nrow = p, ncol = n_components)
  for (i in 1:n_components) {
    loadings[, i] <- eigenvectors[, i] * sqrt(eigenvalues[i])
  }
  rownames(loadings) <- var_names
  colnames(loadings) <- paste0("PC", 1:n_components)

  # Compute scores
  scores_raw <- X_std %*% eigenvectors
  colnames(scores_raw) <- paste0("PC", 1:n_components)
  if (!is.null(rownames(X))) {
    rownames(scores_raw) <- rownames(X)
  }

  # Standardize scores if requested
  if (standardize_scores) {
    scores <- scale(scores_raw, center = FALSE, scale = TRUE)
    colnames(scores) <- paste0("PC", 1:n_components)
    if (!is.null(rownames(X))) {
      rownames(scores) <- rownames(X)
    }
    scores_label <- "Standardized Component Scores (variance = 1)"
  } else {
    scores <- scores_raw
    scores_label <- "Raw Component Scores (variance = eigenvalue)"
  }

  # Compute variance explained
  total_variance <- sum(eigen_decomp$values)
  prop_variance <- eigenvalues / total_variance
  cum_variance <- cumsum(prop_variance)

  # Compute communalities
  communalities <- rowSums(loadings^2)
  names(communalities) <- var_names

  # Create output tables
  cat("===============================================\n")
  cat("PCA RESULTS (Stata/SPSS Style Output)\n")
  cat("===============================================\n\n")

  # Eigenvalues table
  cat("EIGENVALUES AND VARIANCE EXPLAINED\n")
  cat("-----------------------------------\n")
  eigen_table <- data.frame(
    Component = paste0("PC", 1:n_components),
    Eigenvalue = round(eigenvalues, 4),
    Variance_Pct = round(100 * prop_variance, 2),
    Cumulative_Pct = round(100 * cum_variance, 2)
  )
  print(eigen_table, row.names = FALSE)
  cat("\n")

  # Kaiser criterion
  n_kaiser <- sum(eigenvalues > 1)
  cat(paste("Components with eigenvalue > 1 (Kaiser criterion):", n_kaiser, "\n\n"))

  # Loadings matrix
  cat("COMPONENT LOADINGS (Variable-Component Correlations)\n")
  cat("----------------------------------------------------\n")
  loadings_df <- as.data.frame(loadings)
  n_show <- min(10, p)
  print(round(loadings_df[1:n_show, , drop = FALSE], 3))
  if (p > 10) cat("... [", p - 10, " more variables]\n", sep = "")
  cat("\n")

  # Communalities
  cat("COMMUNALITIES (Variance Explained per Variable)\n")
  cat("-----------------------------------------------\n")
  comm_df <- data.frame(
    Variable = var_names[1:n_show],
    Communality = round(communalities[1:n_show], 3)
  )
  print(comm_df, row.names = FALSE)
  if (p > 10) cat("... [", p - 10, " more variables]\n", sep = "")
  cat("\n")

  # Score statistics
  cat(paste0(scores_label, "\n"))
  cat("-----------------------------------------------\n")
  score_summary <- data.frame(
    Component = paste0("PC", 1:min(5, n_components)),
    Mean = round(colMeans(scores[, 1:min(5, n_components), drop = FALSE]), 6),
    StdDev = round(apply(scores[, 1:min(5, n_components), drop = FALSE], 2, sd), 4),
    Min = round(apply(scores[, 1:min(5, n_components), drop = FALSE], 2, min), 3),
    Max = round(apply(scores[, 1:min(5, n_components), drop = FALSE], 2, max), 3)
  )
  print(score_summary, row.names = FALSE)
  cat("\n")

  # Verification
  cat("VERIFICATION\n")
  cat("------------\n")
  actual_cor <- cor(X_std[, 1], scores_raw[, 1])
  expected_cor <- loadings[1, 1]
  cat("Correlation between", var_names[1], "and PC1:\n")
  cat("  Computed from data: ", round(actual_cor, 4), "\n")
  cat("  From loadings matrix: ", round(expected_cor, 4), "\n")
  cat("  Match: ", ifelse(abs(actual_cor - expected_cor) < 0.001, "YES", "NO"), "\n\n")

  # Return results with all names properly set
  result <- list(
    eigenvalues = eigenvalues,
    variance_explained = data.frame(
      component = 1:n_components,
      eigenvalue = eigenvalues,
      prop_variance = prop_variance,
      cum_variance = cum_variance
    ),
    loadings = loadings,
    eigenvectors = eigenvectors,
    scores = scores,
    scores_raw = scores_raw,
    communalities = communalities,
    correlation_matrix = R,
    standardized_data = X_std
  )

  class(result) <- "bpca"
  invisible(result)
}
