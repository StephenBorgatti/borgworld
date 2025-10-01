#' Principal Component Analysis
#'
#' Performs PCA on numeric variables in a data frame
#'
#' @param data A data frame containing numeric variables for PCA
#' @param id Character string naming the column to use as labels (optional)
#' @param scale Logical, whether to scale variables to unit variance (default = TRUE)
#'
#' @return An object of class "bpca" containing:
#'   - scores: The principal component scores (same as x)
#'   - loadings: The loadings matrix (same as rotation)
#'   - correlation: Correlation matrix of the numeric variables
#'   - var_explained: Percentage of variance explained by each PC
#'   - sdev: Standard deviations of the principal components
#'   - numeric_data: The numeric data used for PCA
#'   - All standard prcomp components (x, rotation, center, scale, sdev)
#'
#' @importFrom stats prcomp cor complete.cases
#' @export
#'
#' @examples
#' # Without ID variable
#' p <- bpca(iris[,1:4])
#' summary(p)
#'
#' # With ID variable - exclude it from PCA but use for labels
#' data_with_id <- cbind(id = paste0("S", 1:150), iris[,1:4])
#' p <- bpca(data_with_id, id = "id")
#'
#' # Use with bpcabiplot for visualization
#' # bpcabiplot(p)
#'
bpca <- function(data, id = NULL, scale = TRUE) {

  # Store original data
  original_data <- data

  # Handle id variable
  labels <- NULL
  if (!is.null(id)) {
    if (!(id %in% names(data))) {
      stop(paste("Column", id, "not found in data"))
    }
    # Extract labels before removing from data
    labels <- data[[id]]
    # Remove id column from data for PCA
    data <- data[, names(data) != id, drop = FALSE]
  }

  # Check that remaining columns are numeric
  numeric_cols <- sapply(data, is.numeric)
  if (!all(numeric_cols)) {
    non_numeric <- names(data)[!numeric_cols]
    warning(paste("Removing non-numeric columns:", paste(non_numeric, collapse = ", ")))
    data <- data[, numeric_cols, drop = FALSE]
  }

  # Check for sufficient numeric columns
  if (ncol(data) < 2) {
    stop("Need at least 2 numeric variables for PCA")
  }

  # Remove any rows with missing values
  complete_rows <- complete.cases(data)
  if (sum(!complete_rows) > 0) {
    warning(paste("Removing", sum(!complete_rows), "rows with missing values"))
    data <- data[complete_rows, , drop = FALSE]
    if (!is.null(labels)) {
      labels <- labels[complete_rows]
    }
  }

  # Perform PCA
  pca_result <- prcomp(data, scale = scale, center = TRUE)

  # Add labels to the scores
  if (!is.null(labels)) {
    rownames(pca_result$x) <- labels
  }

  # Calculate variance explained
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

  # Calculate correlation matrix
  cor_matrix <- cor(data)

  # Create result object with all components
  result <- pca_result
  result$var_explained <- var_explained
  result$id_column <- id
  result$labels <- labels
  result$numeric_data <- data  # Store the data actually used for PCA
  result$correlation <- cor_matrix  # Store correlation matrix
  result$scores <- pca_result$x  # Explicit alias for scores
  result$loadings <- pca_result$rotation  # Explicit alias for loadings

  class(result) <- c("bpca", class(result))

  return(result)
}

#' Print method for bpca objects
#'
#' @param x A bpca object
#' @param ... Additional arguments (ignored)
#'
#' @export
#' @method print bpca
print.bpca <- function(x, ...) {
  cat("\nPrincipal Component Analysis Results\n")
  cat("=====================================\n")

  cat("\nData dimensions:", nrow(x$scores), "observations,", ncol(x$loadings), "variables\n")

  if (!is.null(x$id_column)) {
    cat("ID column used for labels:", x$id_column, "(excluded from analysis)\n")
  }

  cat("\nStandard deviations (1, .., p=", length(x$sdev), "):\n", sep = "")
  print(round(x$sdev, 4))

  cat("\nProportion of Variance:\n")
  print(round(x$var_explained[1:min(5, length(x$var_explained))], 2))

  if (length(x$var_explained) > 5) {
    cat("(showing first 5 components)\n")
  }

  cat("\nCumulative Proportion:\n")
  cum_var <- cumsum(x$var_explained)
  print(round(cum_var[1:min(5, length(cum_var))], 2))

  # Show correlation matrix
  cat("\nCorrelation matrix:\n")
  print(round(x$correlation, 3))

  # Show first few PC scores
  cat("\nFirst few PC scores:\n")
  n_show <- min(6, nrow(x$scores))
  pc_show <- min(4, ncol(x$scores))
  print(round(x$scores[1:n_show, 1:pc_show], 3))
  if (nrow(x$scores) > 6 || ncol(x$scores) > 4) {
    cat("(Showing first", n_show, "observations and", pc_show, "PCs)\n")
  }

  invisible(x)
}

#' Summary method for bpca objects
#'
#' @param object A bpca object
#' @param ncomp Number of components to display (default = 7)
#' @param ... Additional arguments (ignored)
#'
#' @export
#' @method summary bpca
summary.bpca <- function(object, ncomp = 7, ...) {
  cat("\nPCA Summary\n")
  cat("===========\n")

  cat("\nCall: bpca with", nrow(object$scores), "observations and",
      ncol(object$loadings), "variables\n")

  cat("\nImportance of components:\n")

  # Create importance table
  importance <- rbind(
    `Standard deviation` = object$sdev,
    `Proportion of Variance` = object$var_explained / 100,
    `Cumulative Proportion` = cumsum(object$var_explained) / 100
  )

  colnames(importance) <- paste0("PC", 1:ncol(importance))

  # Print rounded to 4 decimal places
  n_show <- min(ncomp, ncol(importance))
  print(round(importance[, 1:n_show], 4))

  if (ncol(importance) > n_show) {
    cat("...\n(Showing first", n_show, "components)\n")
  }

  # Show loadings for first few components
  cat("\nLoadings:\n")
  n_comp_show <- min(3, ncol(object$loadings))
  print(round(object$loadings[, 1:n_comp_show], 3))
  if (ncol(object$loadings) > 3) {
    cat("(Showing first 3 components)\n")
  }

  # Determine number of components to retain
  cum_var <- cumsum(object$var_explained)
  n_80 <- which(cum_var >= 80)[1]
  n_90 <- which(cum_var >= 90)[1]

  cat("\nComponents needed for variance explained:\n")
  cat("  80% variance:", n_80, "components\n")
  cat("  90% variance:", n_90, "components\n")

  invisible(object)
}
