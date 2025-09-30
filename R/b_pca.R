#' Principal Component Analysis with Biplot
#'
#' Performs PCA and creates a biplot with optional ID labels
#'
#' @param data A data frame containing numeric variables for PCA
#' @param id Character string naming the column to use as labels (optional)
#' @param scale Logical, whether to scale variables to unit variance (default = TRUE)
#' @param plot Logical, whether to produce a biplot (default = TRUE)
#' @param choices Numeric vector of length 2 indicating which PCs to plot (default = c(1,2))
#'
#' @return An object of class "bpca" containing the prcomp results
#'
#' @importFrom stats prcomp cor
#' @importFrom graphics plot arrows text par
#' @export
#'
#' @examples
#' # Without ID variable
#' bpca(iris[,1:4])
#'
#' # With ID variable - exclude it from PCA but use for labels
#' data_with_id <- cbind(id = paste0("S", 1:150), iris[,1:4])
#' bpca(data_with_id, id = "id")
#'
bpca <- function(data, id = NULL, scale = TRUE, plot = TRUE, choices = c(1, 2)) {

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

  # Create biplot if requested
  if (plot && length(choices) == 2) {
    # Set up plot parameters
    par(mar = c(5, 4, 4, 4))

    # Get scores for chosen PCs
    scores <- pca_result$x[, choices]

    # Get loadings for chosen PCs
    loadings <- pca_result$rotation[, choices]

    # Scale factors for display
    score_scale <- 1
    loading_scale <- max(abs(scores)) / max(abs(loadings)) * 0.8

    # Create axis labels
    xlab <- sprintf("PC%d (%.1f%%)", choices[1], var_explained[choices[1]])
    ylab <- sprintf("PC%d (%.1f%%)", choices[2], var_explained[choices[2]])

    # Plot scores
    plot(scores[,1] * score_scale, scores[,2] * score_scale,
         xlab = xlab, ylab = ylab,
         main = "PCA Biplot",
         pch = 16, col = "steelblue",
         xlim = range(scores[,1]) * 1.2,
         ylim = range(scores[,2]) * 1.2)

    # Add grid
    grid(col = "lightgray", lty = "dotted")
    abline(h = 0, v = 0, col = "gray", lty = 2)

    # Add observation labels if available
    if (!is.null(labels)) {
      text(scores[,1] * score_scale, scores[,2] * score_scale,
           labels = labels,
           pos = 3, cex = 0.7, col = "darkblue")
    }

    # Add loading vectors
    for (i in 1:nrow(loadings)) {
      arrows(0, 0,
             loadings[i,1] * loading_scale,
             loadings[i,2] * loading_scale,
             col = "red", length = 0.1, lwd = 2)

      # Position text at end of arrow with some padding
      text_x <- loadings[i,1] * loading_scale * 1.1
      text_y <- loadings[i,2] * loading_scale * 1.1

      text(text_x, text_y,
           rownames(loadings)[i],
           col = "red", cex = 0.9, font = 2)
    }

    # Add legend
    legend("topright",
           legend = c("Observations", "Variables"),
           col = c("steelblue", "red"),
           pch = c(16, NA),
           lty = c(NA, 1),
           lwd = c(NA, 2),
           cex = 0.8)
  }

  # Create result object
  result <- pca_result
  result$var_explained <- var_explained
  result$id_column <- id
  result$labels <- labels
  result$numeric_data <- data  # Store the data actually used for PCA

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

  cat("\nData dimensions:", nrow(x$x), "observations,", ncol(x$rotation), "variables\n")

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

  # Show correlation matrix of numeric variables only
  cat("\nCorrelation matrix (numeric variables only):\n")
  cor_mat <- cor(x$numeric_data)
  print(round(cor_mat, 3))

  invisible(x)
}

#' Summary method for bpca objects
#'
#' @param object A bpca object
#' @param ... Additional arguments (ignored)
#'
#' @export
#' @method summary bpca
summary.bpca <- function(object, ...) {
  cat("\nImportance of components:\n")

  # Create importance table
  importance <- rbind(
    `Standard deviation` = object$sdev,
    `Proportion of Variance` = object$var_explained / 100,
    `Cumulative Proportion` = cumsum(object$var_explained) / 100
  )

  colnames(importance) <- paste0("PC", 1:ncol(importance))

  # Print rounded to 4 decimal places
  print(round(importance[, 1:min(7, ncol(importance))], 4))

  if (ncol(importance) > 7) {
    cat("...\n(Showing first 7 components)\n")
  }

  invisible(object)
}
